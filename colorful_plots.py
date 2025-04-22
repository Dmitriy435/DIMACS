import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import plotly.graph_objects as go
import numpy as np
import seaborn as sns
import math
import re


q = 8
representation = "pi0 pi2 rho5"

# So far, built only for 1 cusp



def eval_sqrt(e):
    """Evaluates a string expression containing 'sqrt'."""
    expression = ""
    if e[-1] == "?":
        expression = e[:-1]
    else:
        expression = e

    # Find all instances of 'sqrt(x)' in the string
    sqrt_matches = re.findall(r"sqrt\((.*?)\)", expression)

    for match in sqrt_matches:
        # Evaluate the square root and replace in the original string
        expression = expression.replace(f"sqrt({match})", str(math.sqrt(float(match))))

    return eval(expression)





basisVecs = []
data = []

with open('data/1cusp_q' + str(q) + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        if row[0] == "Representation":
            basisVecs = row
            # print(row)
        elif row[0] == representation:
            data = row
            # print(row[5])
            break

arrs = []

for j in range(0, q-1):
    arrTemp = np.zeros(((q+1), (q+1)))
    for i in range(0, len(basisVecs)):
        if i != 0 and basisVecs[i][7] == str(j):
            print("INSIDE!")
            arrTemp[int(basisVecs[i][1])][int(basisVecs[i][4])] = eval_sqrt(data[i])
            print(data[i])
    arrs.append(arrTemp)
    print(arrTemp)


# 2d heat map:
# plt.imshow(arr0, cmap='hot', interpolation='nearest')


# fig, ax = plt.subplots()
# ax.colorbar()
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax = sns.heatmap(arr, linewidth=0.5)
# ax.imshow(arr, cmap='hot', interpolation='nearest')
# ax.invert_yaxis()

# fig = plt.figure(figsize=(10, 10)) 
# ax = fig.add_subplot(111, projection='3d')

# color_map = cm.ScalarMappable(cmap=cm.Greens_r) 
# color_map.set_array(arr)

# img = ax.scatter(arr, marker='s', 
#                  s=200, color='green') 
# plt.colorbar(color_map)



# 3d scatter plot:


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

temp = []
for i in range(0, q+1):
    temp.extend((q+1)*[i])
xvals = (q-1)*temp
# xvals = (q-1)*((q+1)*[0] + (q+1)*[1] + (q+1)*[2] + (q+1)*[3] + (q+1)*[4])
yvals = ((q-1)*(q+1)) * list(range(0,q+1))
zvals = []
for i in range(0, q-1):
    zvals.extend(((q+1)*(q+1))*[i])

values = []
for arr in arrs:
    for i in range(0, q+1):
        # print((arr[i]))
        values.extend(arr[i])

# colorvals0 = [*arr0[0], *arr0[1], *arr0[2], *arr0[3], *arr0[4]]
# colorvals1 = [*arr1[0], *arr1[1], *arr1[2], *arr1[3], *arr1[4]]
# colorvals2 = [*arr2[0], *arr2[1], *arr2[2], *arr2[3], *arr2[4]]
# values = colorvals0 + colorvals1 + colorvals2

ax.scatter(xvals, yvals, zvals, c=values, cmap='Dark2', s=200)
# ax.scatter(xvals, yvals, 1, c=colorvals1, cmap='Dark2', s=200)
# ax.scatter(xvals, yvals, 2, c=colorvals2, cmap='Dark2', s=200)


# fig = go.Figure(data=go.Volume(
#     x=xvals.flatten(),
#     y=yvals.flatten(),
#     z=zvals.flatten(),
#     value=values.flatten(),
#     isomin=-0.1,
#     isomax=0.8,
#     opacity=0.1, # needs to be small to see through all surfaces
#     surface_count=21, # needs to be a large number for good volume rendering
#     ))
# fig.show()



plt.show()
fig.savefig('pictures/plot_q' + str(q)+'.png') 