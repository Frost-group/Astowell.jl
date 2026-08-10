using DelimitedFiles
using LinearAlgebra
using Astowell


# Setup particle momenta and positions - currently for two particles, more is too oscillatory to sample
r = [-pi -pi; -pi -pi]
k = [0.0 0.0; 1.0 0.0]

# Allocate 1D line
N = 10^6

data = zeros(N)

for index in 1:N

    # Shift the second particle away from the first, linearly and as much as possible
    r[2, 2] += 2 * pi / N
    r[2, 1] += 2 * pi / N

    # Calculate the difference between trajectories
    data[index] = real(det(A(r, k, a=1.0, N=2)) - det(A(r, k, a=0.0, N=2)))

end

# Save final step bakcflow difference and resulting 1D manifold 
write_ppm("examples/difference.ppm", renderimg(r, sampleimg(r, k, S=1024, a=1.0, N=2) - sampleimg(r, k, S=1024, a=0.0, N=2), S=1024, POW=1.0))
writedlm("examples/manifold.csv", data, ',')

# Calculate Lyapunov exponent - use long time limit, since we are far away from the regions with backflow impact
exponent = 0.0

for index in 2:(N-1)

    # Calculate the relative derivative - I'm pretty sure we should normalizing here like I do
    derivative = (data[index+1] - data[index-1]) / data[index] / 2

    global exponent

    exponent += log(abs(derivative))

end

println(exponent / N)