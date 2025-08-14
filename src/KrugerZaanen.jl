using LinearAlgebra, SpecialPolynomials

# Length is now adjusted in multiples of \sqrt{\frac{m \omega}{\hbar}}
β = 1.0

# Define particle number here of easier modification
N = 6

# Backflow parameters
r0 = 0.05374
a = -15.83

# Define backflow function
η(r) = a / (r^3 / r0^3 + 1)

# Build quantum numbers grid
function K(; KMAX=2)
    k = zeros(N, 2) # kludge static allocation
    n = 1

    for kx ∈ 0:KMAX
        for ky ∈ 0:KMAX
            if kx^2 + ky^2 <= KMAX^2
                k[n, 1] = kx
                k[n, 2] = ky

                n = n + 1
            end
        end
    end

    return k
end

# Random coordinates for the Fermion particles, in the 2D slab from -L/2->L/2 or smaller
function randomcoords()
    -1 / 2 .+ rand(N, 2)
end

function A(r, k)
    A = Matrix{ComplexF64}(undef, N, N)  # Pre-allocate without zeros

    rkj = zeros(Float64, N, 2)
    ηrkj = zeros(Float64, N)
    rj = zeros(Float64, 1, 2)

    for j in 1:N
        if a == 0.0 # No backflow, trivial case. Don't really need to do the full matrix det, but simplifies the code + acts as a crosscheck to do it
            rj = @view r[j, :]
        else

            # Attempt at vectorisation. Twice as fast.
            # Very buggy to write, and difficult to reason with.
            rkj = r[j, :]' .- r[:, :]
            ηrkj = η.(sqrt.(sum(abs2.(rkj), dims=2)))
            ηrkj[j] = 0.0 # avoid self-interaction
            rj = r[j, :]' .+ sum(ηrkj .* rkj, dims=1)

        end

        for i in 1:N

            n1 = floor(Int, k[i, 1])
            n2 = floor(Int, k[i, 2])

            # Use QHO basis functions
            A[i, j] = β / sqrt(2^(n1 + n2) * pi * factorial(n1) * factorial(n2))

            A[i, j] *= Basis(Hermite, n1)(β * rj[1])
            A[i, j] *= Basis(Hermite, n2)(β * rj[2])

            A[i, j] *= exp(-β^2 * rj ⋅ rj / 2)
        end
    end

    return A
end

# a=0.4, start to see nodes appearing around particle posn
# a=0.7, about 50% smooth, 50% fractal
# scan across x and y for the first particle and try and sample the wavefunction
"""
    sampleimg(r, k)

Sample the backflow wavefunction by scanning the position of the first particle
across a 2D grid while keeping other particles fixed.
Function is now threaded 🚀

# Arguments
- `r`: Array of particle positions
- `k`: Array of k-space points
- `S`: Grid size for sampling (S×S points)
- `a`: Backflow strength parameter
    - a=0.0: No backflow (standard Slater determinant)
    - a=0.4: Nodes start appearing around particle positions
    - a=0.7: Approximately 50% smooth, 50% fractal structure

Returns a complex-valued S×S array containing the wavefunction values.
"""
function sampleimg(r, S)
    img = zeros(ComplexF64, S + 1, S + 1)
    xs = collect(enumerate(-1/2:1/S:1/2))
    k = K()

    Threads.@threads for (i, x) in xs
        # Create thread-local copy of coordinates
        r_local = copy(r)
        r_local[1, 1] = x

        for (j, y) in enumerate(-1/2:1/S:1/2)
            r_local[1, 2] = y
            An = A(r_local, k)
            img[i, j] = det(An)
        end
    end

    return img
end

"""
Render the wavefunction as an RGB PNG file.
- Red/blue indicates sign of wavefunction
- Brightness indicates magnitude (with POW controlling contrast)
- Yellow dots show particle positions
Returns: Height × Width × 3 array of Float64 values between 0 and 1
"""
function renderimg(r, img, S, POW)
    # Define colors
    red = RGB(1.0, 0.0, 0.0)    # positive values
    blue = RGB(0.0, 0.0, 1.0)   # negative values
    yellow = RGB(1.0, 1.0, 0.0) # particles

    rgb = Array{RGB}(undef, S + 1, S + 1)
    MAX = maximum(abs.(img))

    # Render wavefunction
    for i in 1:S+1, j in 1:S+1
        val = (abs(img[i, j]) / MAX)^POW
        rgb[i, j] = real(img[i, j]) > 0 ? val * red : val * blue
    end

    # Add particle positions
    if S > 50
        particlecoords = S * (r .+ 1 / 2) .|> ceil .|> Int
        for r in eachrow(particlecoords)
            for dx in 0:1, dy in 0:1
                x, y = r[1] + dx, r[2] + dy
                if 1 ≤ x ≤ S + 1 && 1 ≤ y ≤ S + 1
                    rgb[x, y] = yellow
                end
            end
        end
    end

    return rgb
end

