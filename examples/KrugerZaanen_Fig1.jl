using Astowell

k=K(N=49)
r=randomcoords(N=49)

img=sampleimg(r,k, S=250, a=0.0) # a is backflow strength

rgb=renderimg(r,img,S=250, POW=0.02)

write_ppm("KrugerZaanen_Fig1.ppm", rgb) 
