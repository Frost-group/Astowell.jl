using Astowell

N=5
k=K(N=N, KMAX=1)

r=randomcoords(N=N)

img=sampleimg(r,k, N=N, S=250, a=0.5) # a is backflow strength
rgb=renderimg(r,img,S=250, POW=0.02)
write_ppm("2particle.ppm", rgb) 

# nb: this currently moves particles 2 and 3 to the x=0 slice; should prob move that
# functionality here
slice = sampleslice(r,k, N=N, S=100, a=0.0)
@show slice

using Gnuplot


for a in 0:0.25:1.0
    slice=sampleslice(r,k, N=N, S=500, a=a)
# slice is array of complex numbers; result of det() of sampled wavefn
    @gp :- real.(slice) "w l title $(a)"
end

Gnuplot.save("slice.png", term="pngcairo size 550,350 fontscale 0.8")


