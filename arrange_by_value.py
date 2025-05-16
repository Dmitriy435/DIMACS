import csv

q = 4
numCus = 1

file = "data/" + str(numCus) + 'cusp_q' + str(q) + ".csv"
#file = "data/" + str(numCus) + 'cusp_q' + str(q) + "/compiled_data.csv"

d = {}

with open(file, 'r') as csvfile:
    header = []
    for row in csv.reader(csvfile):
        if row[0][0] == 'R':
            header = row
        else:
            rep = row[0]
            for i in range(1, len(row)):
                entry = row[i]
                if entry != '':
                    if entry in d.keys():
                        d[entry].append((rep, header[i]))
                    else:
                        d[entry] = [(rep, header[i])]

print(d.keys())

with open("data/" + str(numCus) + 'cusp_q' + str(q) + "_unique_values.csv",'w') as f:
    w = csv.writer(f)
    w.writerows(d.items())