import math

def calculateAverage(listInput):
    average = 0.0
    for i in range(len(listInput)):
        average += listInput[i]
    average /= len(listInput)
    return average

def calculateStandardDeviation(listInput, average):
    standardDeviation = 0.0
    for i in range(len(listInput)):
        standardDeviation += (listInput[i] - average) * (listInput[i] - average)
    standardDeviation /= len(listInput)
    standardDeviation = math.sqrt(standardDeviation)
    return standardDeviation

def GaussianErrorPropagation(errorInput):
    error = 0.0
    for i in range(len(errorInput)):
        error += errorInput[i] * errorInput[i]
    error = math.sqrt(error)
    error /= len(errorInput)
    return error

f = open("freeEnergySummary.txt", "r")

freeEnergy = []
standardDeviation = []

lines = f.readlines()

for line in lines:
    data = line.split()
    if data[0] == "run_":
        freeEnergy.append(float(data[3]))
        standardDeviation.append(float(data[5]))

maxCount = int(len(freeEnergy)/2)

freeEnergyFull = freeEnergy[0:maxCount]
standardDeviationFull = standardDeviation[0:maxCount]
freeEnergyInteractionEntropy = freeEnergy[maxCount:2*maxCount]
standardDeviationInteractionEntropy = standardDeviation[maxCount:2*maxCount]

m = calculateAverage(freeEnergyFull)
se = calculateStandardDeviation(freeEnergyFull, m)
e = GaussianErrorPropagation(standardDeviationFull)
print("%10.2f %10.2f %10.2f" % (m, se, e))

m = calculateAverage(freeEnergyInteractionEntropy)
se = calculateStandardDeviation(freeEnergyInteractionEntropy, m)
e = GaussianErrorPropagation(standardDeviationInteractionEntropy)
print("%10.2f %10.2f %10.2f" % (m, se, e))

f.close()
