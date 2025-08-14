using Astowell

coordinates = randomcoords()

image = renderimg(coordinates, sampleimg(coordinates, 128), 128, 0.01)

write_ppm("figure.ppm", image)
