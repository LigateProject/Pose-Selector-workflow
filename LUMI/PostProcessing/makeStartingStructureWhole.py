inputFile = open("input.gro", "r")
outputFile = open("intermediary.gro", "w")

lines = inputFile.readlines()

# .gro is a fixed file format
#  1- 5: residue number
#  6-10: residue name
# 11-15: atom name
# 16-20: atom number
# 21-28: x-coordinate
# 29-36: y-coordinate
# 37-44: z-coordinate

# relevant cell dimensions
dims = [float(lines[-1].split()[0]), float(lines[-1].split()[1]), float(lines[-1].split()[2])]

residuesToBeConsidered = ["ARG", "HIS", "LYS", "ASP", "GLU", "ASN", "GLN", "SER", "THR", "CYS", "GLY", "PRO", "ALA", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "VAL", "ACE", "LIG"]

count = 0

centreOfGeometry=[[0.0, 0.0, 0.0]]
centreOfGeometryIons=[]

# make individual chains whole again

# the name of the system and the number of atoms will not change => just copy
outputFile.write(lines[0])
outputFile.write(lines[1])

for i in range(2, len(lines)-1):
    # read in coordinates
    data=[float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][37:44])]

    # non-protein and non-ligand atoms must be complexated ions
    # all ions are considered as individual chains
    if lines[i][5:10].strip() not in residuesToBeConsidered:
        centreOfGeometryIons.append([0.0, 0.0, 0.0])
        for j in range(len(data)):
            centreOfGeometryIons[-1][j] += data[j]
        continue

    # inside a chain, we can just use the previous atom as reference and make sure that all distances are smaller than half of the relevant cell dimension
    if count > 0:
        for j in range(len(data)):
            diff = data[j] - previous[j]
            data[j] -= round(diff/dims[j]) * dims[j]
    previous = data
    count += 1

    # for inter-chain distance
    for j in range(len(data)):
        centreOfGeometry[-1][j] += data[j]

    # calculate centre of geometry and reset count at the end of the chain
    if "OC2" in lines[i] and not "LIG" in lines[i]:
        for j in range(len(centreOfGeometry[-1])):
            centreOfGeometry[-1][j] /= count
        centreOfGeometry.append([0.0, 0.0, 0.0])
        count = 0

    # print new coordinates
    outputFile.write("%s%8.3f%8.3f%8.3f\n" % (lines[i][0:20], data[0], data[1], data[2]))

# calculate final centre of geometry for non-peptide ligand
if count > 0:
    for j in range(len(centreOfGeometry[-1])):
        centreOfGeometry[-1][j] /= count
    count = 0
    centreOfGeometry.append([0.0, 0.0, 0.0])

# the definition of the unit cell will not change, either => just copy, too
outputFile.write(lines[-1])

inputFile.close()
outputFile.close()

# place chains close to each other
# start with the first chain of the protein and expect the distances between the centres of geometry to always be less than half of the cell size
shift = []
# start by placing all protein chains next to each other
for i in range(1, len(centreOfGeometry) - 2):
    shift.append([])
    for j in range(len(centreOfGeometry[i])):
        diff = centreOfGeometry[i][j] - centreOfGeometry[i-1][j]
        shift[-1].append(-round(diff/dims[j]) * dims[j])
        # if you move a chain, update its centre of geometry such that the subsequent calculations yield the correct shifts
        centreOfGeometry[i][j] += shift[-1][j]
    if (shift[-1][0] != 0.0 or shift[-1][1] != 0.0 or shift[-1][2] != 0.0):
        print("Need to move protein chain %d to obtain a compact representation" % (i+1))
# calculate centre of geometry of the whole protein
centreOfGeometryProtein = [0.0, 0.0, 0.0]
for i in range(0, len(centreOfGeometry) - 2):
    for j in range(len(centreOfGeometry[i])):
        centreOfGeometryProtein[j] += centreOfGeometry[i][j]
for j in range(len(centreOfGeometryProtein)):
    centreOfGeometryProtein[j] /= (len(centreOfGeometry) - 2)
# calculate shift for the complexated ions using the centre of geometry of the protein
for i in range(len(centreOfGeometryIons)):
    shift.append([])
    for j in range(len(centreOfGeometryIons[i])):
        diff = centreOfGeometryIons[i][j] - centreOfGeometryProtein[j]
        shift[-1].append(-round(diff/dims[j]) * dims[j])
    if (shift[-1][0] != 0.0 or shift[-1][1] != 0.0 or shift[-1][2] != 0.0):
        print("Need to move ion %d to obtain a compact representation" % (i+1))
# calculate shift for the ligand using the centre of geometry of the protein
shift.append([])
for j in range(len(centreOfGeometry[-2])):
    diff = centreOfGeometry[-2][j] - centreOfGeometryProtein[j]
    shift[-1].append(-round(diff/dims[j]) * dims[j])
if (shift[-1][0] != 0.0 or shift[-1][1] != 0.0 or shift[-1][2] != 0.0):
    print("Need to move ligand to obtain a compact representation")

inputFile = open("intermediary.gro", "r")
outputFile = open("output.gro", "w")

chainCount = -1

# the name of the system and the number of atoms will not change => just copy
outputFile.write(lines[0])
outputFile.write(lines[1])

for i in range(2, len(lines)-1):
    # read in coordinates
    data=[float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][37:44])]

    if chainCount > -1:
       for j in range(len(data)):
           data[j] += shift[chainCount][j]

    # increase chainCount at the end of the chain or if you read in an ion
    if ("OC2" in lines[i] and not "LIG" in lines[i]) or lines[i][5:10].strip() not in residuesToBeConsidered:
        chainCount += 1

    # print new coordinates
    outputFile.write("%s%8.3f%8.3f%8.3f\n" % (lines[i][0:20], data[0], data[1], data[2]))

# the definition of the unit cell will not change, either => just copy, too
outputFile.write(lines[-1])

inputFile.close()
outputFile.close()
