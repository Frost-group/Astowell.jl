module Astowell

using Random
Random.seed!(0xDEADBEEF) # reproducible runs

include("SimplePNG.jl")
include("KrugerZaanen.jl")

export write_ppm, RGB
export sampleimg, sampleslice, renderimg, K, randomcoords, η, A

end
