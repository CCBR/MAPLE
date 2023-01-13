#DAC.py

import sys


infile = open(sys.argv[1], 'rt')
limit = int(sys.argv[2])
Max_dist = int(sys.argv[3])
length = int(sys.argv[4])
outfile = open(sys.argv[5], 'w+')

nestedlist =  []
DAC = {}

# Initialize output list
for i in range(1,Max_dist+1):
    DAC[i] = 0


for line in infile:
    array = line.strip().split('\t')
    chrom, center, counts = array[0], int(array[1]), float(array[2])

    # Set unrealistic DYAD counts to 0
    if counts > limit:
        counts = 0

    nestedlist.append([chrom,center,counts])

#testout = open('nested_DAC.csv','w+')
#for item in nestedlist:
#    testout.write(str(item)+'\n')


# Next chunks takes an unique position and a set a consecutive position
# It then finds the distance between the unique position and each position in the set
# The product of counts is added the DAC[distance]
# It then takes the next unique position and repeats the process until exhaustion of positions

i = 0
while i < len(nestedlist):
    currentitem = nestedlist[i]
    currentlist = nestedlist[i+1:i+Max_dist+1]
    i += 1

    for item in currentlist:
        if item[0] == currentitem[0]: #are they on the same chromosome?

            dist = float(item[1]) - float(currentitem[1])
            if dist > Max_dist:
                break
            else:
                if dist > Max_dist or dist < 0:
                    print('Dist Error')
                    print(dist)
                    print(currentitem, item)

                prod = float(currentitem[2]) * float(item[2])

                try:
                    DAC[dist] += prod
                except:
                    print(currentitem, item)

# Create output

outfile.write('Dist,DAC\n')
for dist in DAC:
    try:
        outfile.write(str(dist)+','+str(DAC[dist]/(length - dist + 1))+'\n')
    except: # interval with no DYAD
        outfile.write(str(dist)+',0\n')
