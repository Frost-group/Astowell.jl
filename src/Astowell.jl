module Astowell

using Random

Random.seed!(0xDEADBEEF) # reproducible runs

include("KrugerZaanen.jl")
include("SimplePPM.jl")

export sampleimg, renderimg, K, randomcoords, η, A
export write_ppm, RGB

end
