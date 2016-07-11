infile = open("fort.2", "r")
hist = []
hist.append(0.0)
hist.append(0.0)

for line in infile:
    temp = line.split()
    if int(temp[0]) == 1:
        if float(temp[1]) < float(temp[2]):
            hist[0] = hist[0] + 1.0
        else:
            hist[1] = hist[1] + 1.0


outfile = open("Contribution.txt", "w")
outfile.write( " " + "1" + " " + str(hist[0]) + "\n" )
outfile.write( " " + "2" + " " + str(hist[1]) + "\n" )
