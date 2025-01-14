using Astowell
using ProgressMeter
using Printf

# Setup parameters
N = 49
S = 512  # image resolution
frames = 250
a_range = range(0.0, 1.0, frames)  # backflow strength range

# Initialize system
k = K(N=N)
r = randomcoords(N=N)

# Create output directory if it doesn't exist
mkpath("frames")

println("🎬 Generating Kruger-Zaanen backflow movie frames")
println("Resolution: $(S)×$(S), Frames: $frames")

println("🧵 Using $(Threads.nthreads()) threads")

# Progress bar for the whole process
p = Progress(frames, desc="Rendering frames: ", showspeed=true)

for (frame, a) in enumerate(a_range)
    # Generate wavefunction
    img = sampleimg(r, k, S=S, a=a)
    
    # Render to RGB
    rgb = renderimg(r, img, S=S, POW=0.02)
    
    # Save frame
    filename = @sprintf("frames/backflow_%03d.ppm", frame)
    write_ppm(filename, rgb)
    
    next!(p)
end

println("\n✨ Done! Convert to video with:")
println("# Imagemagick 'convert' PPM to PNG with text overlay")
println("""for f in frames/*.ppm; do
    num=\$(basename \$f .ppm)
    frame=\${num##*_}
    a=\$(printf "%.2f" \$(echo "scale=2; (\$frame - 1) * (1.0/$(frames-1))" | bc))
    convert \$f -pointsize 36 -fill white -annotate +50+$(S-50) "a = \$a" \${f%.*}.png
done""")
println("ffmpeg -framerate 10 -i frames/backflow_%03d.png -c:v libx264 backflow.mp4")
println("rm frames/*.ppm")
