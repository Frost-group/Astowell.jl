module Astowell

using Random

include("KrugerZaanen.jl")
include("SimplePPM.jl")

export sampleimg, renderimg, K, randomcoords, η, A
export write_ppm, RGB

end
