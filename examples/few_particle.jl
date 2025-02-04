using Astowell

N=5

k=K(N=N, KMAX=1)

r=randomcoords(N=N)

img=sampleimg(r,k, N=N, S=250, a=0.5) # a is backflow strength

rgb=renderimg(r,img,S=250, POW=0.02)

write_ppm("2particle.ppm", rgb) 
