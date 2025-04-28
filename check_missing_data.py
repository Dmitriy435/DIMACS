import csv
import os

q = 5
numCus = 1


folder = "data/" + str(numCus) + 'cusp_q' + str(q)

repcombos = []
with open("data/" + str(numCus) + 'cusp_q' + str(q) + "/representation_combinations.csv", 'r') as csvfile:
    for row in csv.reader(csvfile):
        repcombos.append(row)

V1 = q+1 if numCus < 3 else q-1
V2 = q+1 if numCus < 2 else q-1
V3 = q+1 if numCus < 1 else q-1
numBasisCombos = V1 * V2 * V3

missingData = {}


for i in range(0, len(repcombos)):
    if str(i).zfill(3) in os.listdir(folder):
        missingvals = []
        for j in range(0, numBasisCombos):
            if (str(j).zfill(4) + ".csv") not in os.listdir(folder + "/" + str(i).zfill(3)):
                #print("Missing basis combo " + str(j) + " in representation combo " + str(i))
                missingvals.append(j)
        missingData[i] = missingvals
    else:
        #print("Missing all of representation combo " + str(i))
        missingData[i] = "All"

#print(missingData)

with open(folder + '/missing_data.csv', 'w') as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in missingData.items():
       writer.writerow([key, value])