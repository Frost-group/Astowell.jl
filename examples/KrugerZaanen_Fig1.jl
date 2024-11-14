using Astowell

k=K(N=49)
r=randomcoords(N=49)

img=sampleimg(r,k, S=250, a=0.5) # a is backflow strength

renderimg(r,img,S=250, POW=0.02)

# How to save this sensibly?
