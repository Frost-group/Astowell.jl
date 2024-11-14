using LinearAlgebra
# Test with no backflow
r_test = randomcoords(N=49)
k_test = K(N=49)

An_no_backflow = A(r_test, k_test, a=0.0)
@test size(An_no_backflow) == (49,49)
@test eltype(An_no_backflow) == ComplexF64
# Test determinant is non-zero complex number
d = det(An_no_backflow)
@test !iszero(d)
@test d isa ComplexF64

# Test with backflow
An_backflow = A(r_test, k_test, a=0.5) 
@test size(An_backflow) == (49,49)
@test eltype(An_backflow) == ComplexF64

# Test determinant is non-zero complex number
d = det(An_backflow)
@test !iszero(d)
@test d isa ComplexF64

