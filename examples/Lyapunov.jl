using DelimitedFiles
using LinearAlgebra
using Astowell


# Setup particle momenta and positions - currently for two particles, more is too oscillatory to sample
r = [-pi -pi; -pi -pi]
k = [0.0 0.0; 1.0 0.0]

# Allocate 1D line
N = 10^3

# Set range of whole box size
η = 1e-3

# Allocate data
data = zeros(N - 1, 2)

# Avoid overflow
r[2, 2] += η * 2 * pi / N
r[2, 1] += η * 2 * pi / N

for index in 1:(N-1)

    Dₙ₋₁ = real(det(A(r, k, a=1.0, N=2)))

    # Shift the second particle away from the first, linearly and as much as possible
    r[2, 2] += η * 2 * pi / N
    r[2, 1] += η * 2 * pi / N

    # Calculate the distance between two particles
    xₙ = sqrt((r[2, 2] - r[1, 2])^2 + (r[2, 1] - r[1, 1])^2) / (2 * pi)

    Dₙ = real(det(A(r, k, a=1.0, N=2)))

    # Shift again
    r[2, 2] += η * 2 * pi / N
    r[2, 1] += η * 2 * pi / N

    Dₙ₊₁ = real(det(A(r, k, a=1.0, N=2)))

    # Calculate derivative
    derivative = (Dₙ₊₁ - Dₙ₋₁) / 2

    # Adjust for a relative change
    derivative /= Dₙ

    # Adjust for the length of the box to be a unit
    derivative *= N / η

    # Save the logarithm of the derivative for exponent calculation
    data[index, 2] = log(abs(derivative))

    # Save coordinate value
    data[index, 1] = xₙ

    # Shift back positions
    r[2, 2] -= η * 2 * pi / N
    r[2, 1] -= η * 2 * pi / N

end

# Save final step backflow and resulting 1D manifold 
write_ppm("examples/backflow.ppm", renderimg(r, sampleimg(r, k, S=1024, a=1.0, N=2), S=1024, POW=1.0))
writedlm("examples/manifold.csv", data, ',')