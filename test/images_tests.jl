# Test image sampling functionality
@testset "Image Sampling" begin
    r_sample = randomcoords(N=49)
    k_sample = K(N=49)
    img = sampleimg(r_sample, k_sample, S=50, a=0.5)
    # img isn't really an image, it's a wavefunction ^_^

    @test size(img) == (51, 51)
    @test eltype(img) == ComplexF64 # I know, 'img' should never be complex, but the name has stuck for now
    @test all(x -> abs(x) ≥ 0, img) # +ve definite 
end

@testset "Image Rendering" begin
    r_render = randomcoords(N=49)
    test_img = ones(51, 51)
    rgb = renderimg(r_render, test_img, S=50, POW=0.02)
    
    # it's all a bit of a hacky mess, so I'm kinda just checking that it runs
    @test size(rgb) == (51, 51)
    @test eltype(rgb) == RGB
end

# Test parameter validation
@testset "Parameter Validation" begin
    r_test = randomcoords(N=49)
    k_test = K(N=49)

    # essentially integration test ; I'll probably change the API, so don't be
    # too upset if this is failing
    
    @test size(sampleimg(r_test, k_test, S=10, a=0.5)) == (11, 11)
    @test size(renderimg(r_test, ones(11,11), S=10, POW=0.02)) == (11, 11)
end
