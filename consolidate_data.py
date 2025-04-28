# from tri_forms_parallelizable, gather all the csv files and combine into one big file (efficiently)

import csv
import numpy as np
from numpy.dtypes import StringDType
import pandas as pd 
from pathlib import Path
import itertools
import os

def cartesian_product(*iterables):
    return list(itertools.product(*iterables))



# Set values:
q = 7
numCus = 1



folder = "data/" + str(numCus) + "cusp_q" + str(q)


numRepCombos = 0
if not Path(folder + "/representation_combinations.csv").exists():
    print("There is no data!!! (or at least no file of combinations of representations)")
else:
    with open(folder + "/representation_combinations.csv", "r") as file:
        csv_reader = csv.reader(file)
        numRepCombos = sum(1 for row in csv_reader)

basisCombos = 0
V1 = q+1 if numCus < 3 else q-1
V2 = q+1 if numCus < 2 else q-1
V3 = q+1 if numCus < 1 else q-1
basisCombos = cartesian_product(range(0, V1), range(0, V2), range(0, V3))
basisCombos = [str(element) for element in basisCombos]
basisCombos.insert(0, "Representation")


# CREATE HEADER
data = np.full((numRepCombos+1, len(basisCombos)), "", dtype=StringDType())
data[0, :] = basisCombos

print(data.dtype)


for subFolder in os.listdir(folder):
    if subFolder[0] == '0' or subFolder[0] == '1':
        t = folder + "/" + subFolder
        for file in os.listdir(t):
            with open(t+"/"+file, 'r') as csvfile:
                r = csv.reader(csvfile)
                for d in r:
                    row = int(subFolder) + 1
                    column = int(file[:-4]) + 1
                    data[row, column] = d[1]
                    data[row, 0] = d[0]


print(data)
# np.savetxt(folder + "/compiled_data.csv", data, delimiter=",") # Gives errors for some reason?
df = pd.DataFrame(data) # this uses up more than necessary memory, idk if this is a problem?
df.to_csv(folder + "/compiled_data.csv", header=False, index=False)