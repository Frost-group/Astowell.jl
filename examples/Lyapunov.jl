using DelimitedFiles
using LinearAlgebra
using Astowell


# Setup particle momenta and positions - currently for two particles, more is quite noisy
r = [0.0 0.0; 0.0 0.0]
k = [0.0 0.0; 1.0 0.0]

# Allocate 1D line
N = 10^6

data = zeros(N)

for index in 1:N

    # Shift the second particle away from the origin, linearly
    r[2, 2] += pi / N
    r[2, 1] += pi / N

    data[index] = real(det(A(r, k, a=1.0, N=2)))

end

# Save 1D manifold
writedlm("line.csv", data, ',')

# Calculate Lyapunov exponent
exponent = 0.0

for index in 2:(N-1)

    # Calculate the relative derivative - I'm pretty sure we should normalizing here like I do
    derivative = (data[index+1] - data[index-1]) / data[index] / 2

    global exponent

    exponent += log(abs(derivative))

end

println(exponent / N)