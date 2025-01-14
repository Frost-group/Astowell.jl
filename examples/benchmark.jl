using BenchmarkTools
using Astowell
using LinearAlgebra

N=49

# Setup common test data
const r_test = randomcoords(N=N)
const k_test = K(N=N)

println("🚀 Kruger-Zaanen Benchmarks\n")

println("Backflow Matrix Construction: A(randomcoords($N), K($N))")
println("--------------------------------")
println("Without backflow: (a=0.0)")
display(@benchmark A($r_test, $k_test, a=0.0))
println("With backflow: (a=0.5)")
display(@benchmark A($r_test, $k_test, a=0.5))
println()

#println("Image Generation:")
#for S in [50, 100, 200]
#    println("S = $S:")
#    display(@benchmark sampleimg($r_test, $k_test, S=$S, a=0.5))
#end

