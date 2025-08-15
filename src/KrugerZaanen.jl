using LinearAlgebra, SpecialPolynomials


# Build quantum numbers grid
function Q(N, n_max)

    # kludge static allocation
    n_numbers = zeros(N, 2)
    n = 1

    for n_x ∈ 0:n_max
        for n_y ∈ 0:n_max
            if n_x + n_y <= n_max
                n_numbers[n, 1] = n_x
                n_numbers[n, 2] = n_y

                n = n + 1
            end
        end
    end

    return n_numbers
end

# Random coordinates for the Fermion particles, in the 2D slab from -1/2->1/2, first one is not used
function randomcoords(; N=3, scale=1.0)
    scale * (-1 / 2 .+ rand(N, 2))
end

function A(r, k, a, r0, β, N)

    # Define backflow function
    η(r) = a / (r^3 / r0^3 + 1)

    # Pre-allocate determinant without zeros
    A = Matrix{ComplexF64}(undef, N, N)

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

# Sample the wavefunction, length is now adjusted in multiples of \sqrt{\frac{m \omega}{\hbar}} via β
function sampleimg(r, S; a=0.0, r0=1.0, β=1.0, n_max=1)

    # Allocate pixels
    img = zeros(ComplexF64, S + 1, S + 1)
    xs = collect(enumerate(-1/2:1/S:1/2))

    # Calculate particle number
    N = size(r, 1)

    # Calculate quantum numbers
    q = Q(N, n_max)

    Threads.@threads for (i, x) in xs

        # Create thread-local copy of coordinates
        r_local = copy(r)
        r_local[1, 1] = x

        for (j, y) in enumerate(-1/2:1/S:1/2)
            r_local[1, 2] = y

            An = A(r_local, q, a, r0, β, N)
            img[i, j] = det(An)
        end

    end

    return img
end

# Render the sampled grid
function renderimg(r, img, S; POW=1.0)

    # Calculate particle number
    N = size(r, 1)

    # Define colors
    red = RGB(1.0, 0.0, 0.0)    # positive values
    blue = RGB(0.0, 0.0, 1.0)   # negative values
    yellow = RGB(1.0, 1.0, 0.0) # particles

    rgb = Array{RGB}(undef, S + 1, S + 1)
    MAX = maximum(abs.(img))

    # Render wavefunction (POW creates an unrealistic sampling approach, just for fun?)
    for i in 1:S+1, j in 1:S+1
        val = (abs(img[i, j]) / MAX)^POW
        rgb[i, j] = real(img[i, j]) > 0 ? val * red : val * blue
    end

    # Add particle positions, except the first one
    particlecoords = S * (r .+ 1 / 2) .|> ceil .|> Int

    for n in 2:N
        for dx in 0:1, dy in 0:1
            x, y = particlecoords[n, 1] + dx, particlecoords[n, 2] + dy
            if 1 ≤ x ≤ S + 1 && 1 ≤ y ≤ S + 1
                rgb[x, y] = yellow
            end
        end
    end

    return rgb
end

