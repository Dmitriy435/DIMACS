import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import plotly.graph_objects as go
import numpy as np
import seaborn as sns
import math
import re


def eval_sqrt(expression):
    """Evaluates a string expression containing 'sqrt'."""

    # Find all instances of 'sqrt(x)' in the string
    sqrt_matches = re.findall(r"sqrt\((.*?)\)", expression)

    for match in sqrt_matches:
        # Evaluate the square root and replace in the original string
        expression = expression.replace(f"sqrt({match})", str(math.sqrt(float(match))))

    return eval(expression)



basisVecs = []
data = []

with open('1cusp_final_q4.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        if row[0] == "Representation":
            basisVecs = row
            # print(row)
        elif row[0] == "pi0 pi1 rho1":
            data = row
            # print(row[5])
            break

arr0 = np.zeros((5, 5))
for i in range(0, len(basisVecs)):
    if i != 0 and basisVecs[i][7] == '0':
        arr0[int(basisVecs[i][1])][int(basisVecs[i][4])] = eval_sqrt(data[i])
        print(data[i])
print(arr0)

arr1 = np.zeros((5, 5))
for i in range(0, len(basisVecs)):
    if i != 0 and basisVecs[i][7] == '1':
        arr1[int(basisVecs[i][1])][int(basisVecs[i][4])] = eval_sqrt(data[i])
        print(data[i])
print(arr1)

arr2 = np.zeros((5, 5))
for i in range(0, len(basisVecs)):
    if i != 0 and basisVecs[i][7] == '2':
        arr2[int(basisVecs[i][1])][int(basisVecs[i][4])] = eval_sqrt(data[i])
        print(data[i])
print(arr2)


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

xvals = 3*(5*[0] + 5*[1] + 5*[2] + 5*[3] + 5*[4])
yvals = 15 * list(range(0,5))
zvals = 25*[0] + 25*[1] + 25*[2]
colorvals0 = [*arr0[0], *arr0[1], *arr0[2], *arr0[3], *arr0[4]]
colorvals1 = [*arr1[0], *arr1[1], *arr1[2], *arr1[3], *arr1[4]]
colorvals2 = [*arr2[0], *arr2[1], *arr2[2], *arr2[3], *arr2[4]]
values = colorvals0 + colorvals1 + colorvals2

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
fig.savefig('test.png') 