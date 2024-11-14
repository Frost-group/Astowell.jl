using LinearAlgebra, Images

# 2D backflow wavefunction node visualiser
# Following PRB 78 035104 (2008) - Fermionic quantum criticality
# https://doi.org/10.1103/PhysRevB.78.035104
# Physics big love.

KMAX=4 # really this should decide N below
N=49

# back to front to how we normally do it? K running from +- 4.0 currently
L=2*pi # also not sure of the maths here; k space normally +- pi/a

"""Build a two dimensional reciprocal space. The 2008 PRB chooses N=49 as it makes a nice regular KMAX<=4 space in 2D."""
function K(;N=49)
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
	A=zeros(ComplexF64, N,N)
	rkj=zeros(Float64,N,2)
	ηrkj=zeros(Float64,N)
	rj=zeros(Float64,1,2)

	
	for j in 1:N
		if a==0.0 # No backflow, trivial case. Don't really need to do the full matrix det, but simplifies the code + acts as a crosscheck to do it
			rj=r[j,:]
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
			A[i,j]=exp(im*k[i,:]⋅rj) 
		end
	end
return A
end



# scan across x and y for the first particle and try and sample the wavefunction
function sampleimg(r,k;S=100, a=0.0)
	# S size in X and Y to raster across
	# S=100 ==> 8.8 s with backflow, on my M1
	#       ==> 2.6 s with no backflow
	
	img=zeros(ComplexF64, S+1,S+1)
	for (i,x) in enumerate(-L/2:L/S:L/2) 
		r[1,1]=x
		for (j,y) in enumerate(-L/2:L/S:L/2)
			r[1,2]=y
	
			An=A(r,k, a=a)
			# a=0.4, start to see nodes appearing around particle posn
			# a=0.7, about 50% smooth, 50% fractal
			img[i,j]=det(An)
		end
	end
	
	return img	
end

function renderimg(r,img; S=100, POW=0.07)
	p = RGB{N0f8}(1.0, 0.0, 0.0) 
	n = RGB{N0f8}(0.0, 0.0, 1.0)  
	
	rgb = similar(img, RGB{N0f8})

	MAX=maximum(abs.(img))
	
	# Broadcast the conditional logic and color assignment
  	for i in eachindex(img)
    	rgb[i] = real(img[i]) > 0 ? p : n
		rgb[i]*= (abs(img[i])/MAX)^POW # map to 0..1, then take 1/x power to squash up to 1. Otherwise everything is BLACK except for a small region.
  	end

	if (S>50) # only plot electron positions if there is space
	# plot particle positions as dots
	# This would be a lot easier if I was just using a plotting package
	# convert coordinates of particles into x,y within image
		particlecoords = S * (r.+ L/2) / L .|> ceil .|> Int
		for r in eachrow(particlecoords)
#			print(r," ",r[1]," ",r[2],"\n")
			rgb[r[1],r[2]]=RGB{N0f8}(1.0,1.0,0.0) # dot for the particle posn		
			rgb[r[1]+1,r[2]]=RGB{N0f8}(1.0,1.0,0.0) # dot for the particle posn
			rgb[r[1],r[2]+1]=RGB{N0f8}(1.0,1.0,0.0) # dot for the particle posn		
			rgb[r[1]+1,r[2]+1]=RGB{N0f8}(1.0,1.0,0.0) # dot for the particle posn
		end
	end
	
	imresize(rgb, (500,500))
end


