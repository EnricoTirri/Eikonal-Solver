# we need to generate a 16x16 grind in a text file, specifing the edges between the nodes
# ex line 0.0,0.0 1.0,0.0 as a edge between the nodes (0.0,0.0) and (1.0,0.0)

file = open("test1.txt", "w")
# write nodes
file.write("256\n")
# write edges = edges *4 - with *4 *2
file.write(str(16 * 16 * 4 - 16 * 4 * 2) + "\n")

for i in range(16):
    for j in range(16):
        if i < 15:
            file.write(str(i) + ".0," + str(j) + ".0 " + str(i + 1) + ".0," + str(j) + ".0\n")
        if j < 15:
            file.write(str(i) + ".0," + str(j) + ".0 " + str(i) + ".0," + str(j + 1) + ".0\n")

file.close()
