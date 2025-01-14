using Astowell
using Test

@testset "Astowell.jl" begin
    @testset "KrugerZaanen" begin include("KrugerZaanen.jl") end
    @testset "Images" begin include("images_tests.jl") end
end
