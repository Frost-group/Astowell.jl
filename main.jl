using Astowell

coordinates = randomcoords(N=6, scale=0.125)

image = renderimg(coordinates, sampleimg(coordinates, 1024, a=-15.83, r0=0.05374, n_max=2, β=8.0), 1024)

write_ppm("figure.ppm", image)
