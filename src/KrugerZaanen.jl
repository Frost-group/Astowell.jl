using LinearAlgebra

# 2D backflow wavefunction node visualiser
# Following PRB 78 035104 (2008) - Fermionic quantum criticality
# https://doi.org/10.1103/PhysRevB.78.035104
# Physics big love.

#KMAX=4 # really this should decide N below
#N=49

# back to front to how we normally do it? K running from +- 4.0 currently
L=2*pi # also not sure of the maths here; k space normally +- pi/a

"""Build a two dimensional reciprocal space. The 2008 PRB chooses N=49 as it makes a nice regular KMAX<=4 space in 2D."""
function K(;N=49, KMAX=4)
	k=zeros(N,2) # kludge static allocation
	n=1

	for kx ∈ -KMAX:KMAX
		for ky ∈ -KMAX:KMAX
			if kx^2+ky^2<=KMAX^2
				k[n,1]=kx
				k[n,2]=ky
				n=n+1
			end
		end
	end
	
	return k
end

k=K(N=49)

# Random coordinates for the Fermion particles, in the 2D slab from -L/2->L/2
function randomcoords(;N=49)
    -L/2 .+ rand(N,2)*L
end

η(r; a=1.0,r0=0.2)=a^3/(r^3+r0^3) # Backflow function

function A(r,k;N=49, a=0.0)
    A = Matrix{ComplexF64}(undef, N, N)  # Pre-allocate without zeros

	rkj=zeros(Float64,N,2)
	ηrkj=zeros(Float64,N)
	rj=zeros(Float64,1,2)

	for j in 1:N
		if a==0.0 # No backflow, trivial case. Don't really need to do the full matrix det, but simplifies the code + acts as a crosscheck to do it
			rj=@view r[j,:]
		else

# Attempt at vectorisation. Twice as fast.
# Very buggy to write, and difficult to reason with.
			rkj=r[j,:]' .- r[:,:]
			ηrkj=η.(sqrt.(sum(abs2.(rkj), dims=2)), a=a)
			ηrkj[j]=0.0 # avoid self-interaction
			rj= r[j,:]' .+ sum( ηrkj .* rkj, dims=1) 

 #Original explicit for-loop with catch version
# 			rj=r[j,:]
#			for k in 1:N
#				if k!=j
#					rj+=η(norm(r[j,:]-r[k,:]), a=a)*(r[j,:]-r[k,:]) # HOT LOOP, optimise me
#				end
#			end
			
		end
		
		for i in 1:N
#			A[i,j]=exp(im*k[i,:]⋅rj) # works! but is slow by a factor of 2
           @inbounds phase = k[i,1]*rj[1] + k[i,2]*rj[2] # explicit phase calculation
           @inbounds A[i,j] = exp(im*phase)  # exp seems faster than cospi+sinpi
		end
	end
return A
end

# a=0.4, start to see nodes appearing around particle posn
# a=0.7, about 50% smooth, 50% fractal
# scan across x and y for the first particle and try and sample the wavefunction
"""
    sampleimg(r, k; S=100, a=0.0)

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
function sampleimg(r, k; N=49, S=100, a=0.0)
    img = zeros(ComplexF64, S+1, S+1)
    xs = collect(enumerate(-L/2:L/S:L/2))
    
    Threads.@threads for (i,x) in xs
        # Create thread-local copy of coordinates
        r_local = copy(r)
        r_local[1,1] = x
        
        for (j,y) in enumerate(-L/2:L/S:L/2)
            r_local[1,2] = y
            An = A(r_local, k, N=N, a=a)
            img[i,j] = det(An)
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
function renderimg(r,img; S=100, POW=0.07)
    # Define colors
    red = RGB(1.0, 0.0, 0.0)    # positive values
    blue = RGB(0.0, 0.0, 1.0)   # negative values
    yellow = RGB(1.0, 1.0, 0.0) # particles
    
    rgb = Array{RGB}(undef, S+1, S+1)
    MAX = maximum(abs.(img))
    
    # Render wavefunction
    for i in 1:S+1, j in 1:S+1
        val = (abs(img[i,j])/MAX)^POW
        rgb[i,j] = real(img[i,j]) > 0 ? val * red : val * blue
    end

    # Add particle positions
    if S > 50
        particlecoords = S * (r .+ L/2) / L .|> ceil .|> Int
        for r in eachrow(particlecoords)
            for dx in 0:1, dy in 0:1
                x, y = r[1]+dx, r[2]+dy
                if 1 ≤ x ≤ S+1 && 1 ≤ y ≤ S+1
                    rgb[x,y] = yellow
                end
            end
        end
    end
    
    return rgb
end 

