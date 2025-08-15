using Astowell

coordinates = randomcoords()

image = renderimg(coordinates, sampleimg(coordinates, 1024, a=-15.83, r0=0.05374), 1024)

write_ppm("figure.ppm", image)
