"""
Minimal PPM writer with basic RGB type.
No dependencies, just writes 16-bit RGB PPM files (P6 format).
"""
struct RGB
    r::Float64
    g::Float64
    b::Float64
end

Base.:*(c::RGB, f::Number) = RGB(c.r * f, c.g * f, c.b * f)
Base.:*(f::Number, c::RGB) = c * f

"""Write RGB array to 16-bit PPM file"""
function write_ppm(filename, img::Array{RGB,2})
    height, width = size(img)
    
    open(filename, "w") do io
        # PPM header (P6 = binary RGB format)
        write(io, "P6\n$width $height\n65535\n")
        
        # RGB data (16-bit values in binary)
        for y in 1:height, x in 1:width
            # Convert 0-1 float values to 16-bit integers; bswap to Big Endian for PPM
            r = round(UInt16, clamp(img[y,x].r, 0, 1) * 65535) |> bswap
            g = round(UInt16, clamp(img[y,x].g, 0, 1) * 65535) |> bswap
            b = round(UInt16, clamp(img[y,x].b, 0, 1) * 65535) |> bswap
            
            write(io, r, g, b)
        end
    end
end
