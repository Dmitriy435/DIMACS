
# Preliminaries

# Field and general linear group
q = 4
K = GF(q)
G = GL(2, K)

# Borel subgroup
a = K.zeta()
MS = MatrixSpace(K, 2)
gens = [MS([[1,1],[0,1]]),MS([[a,0],[0,1]]), MS([[1,0],[0,a]])]
B = G.subgroup(gens)

# Creating vector space for induced rep
V = VectorSpace(CC, q+1)

# The chosen coset representations 
# First q spots are of the form [[1, 0][-i, 1]] (i counts up from 0 to q-1)
# Last spot (index q) is w
cosetReps = []
repToIndex = {}
index = 0
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetReps.append(rep)
    repToIndex[rep] = index
    #print(rep)
    index = index + 1
rep = G(MS([[0,1],[1,0]]))
cosetReps.append(rep)
repToIndex[rep] = q
#print(cosetReps[q])










# Takes in character of B and returns if the induced representation is irreducible or not

def isInducedIrreducible(chi):
    induced_chi = chi.induct(G)
    return induced_chi.is_irreducible()
ct = B.character_table()
#print(ct)
badChars = []
goodChars = []
for i in range(0, len(B.conjugacy_classes_representatives())):
    if ct[i][0] == 1:
        if isInducedIrreducible(B.character(ct[i])):
            goodChars.append(B.character(ct[i]))
        else:
            badChars.append(B.character(ct[i]))
chi = goodChars[0]
print(chi.values())

#print(len(goodChars))
#print(len(badChars))

# NOTICE - the induced reps of the goodChars are each counted twice, as two diff chars map to the same induced rep







# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentative(g):
    if g.inverse().matrix()[0][0] == 0:
        return G(MS([[0,1],[1,0]])).inverse()
    else:
        x = (K.one() / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        return G(MS([[1,0],[x,1]])).inverse()
# Constant time


# Gives the G action result of group element g from G onto vector v from V
# Uses globally defined character chi
# Uses globally defined vector space V
def gAction(g, vec):
    newVec = V([0] * (q+1))
    for i in range(0, q+1):
        newRep = toRepresentative(cosetReps[i] * g.inverse())
        #print(newRep)
        b = cosetReps[i] * g.inverse() * newRep.inverse()
        
        newIndex = repToIndex[newRep]

        newVec[newIndex] = newVec[newIndex] + chi(b.inverse()) * vec[i]
        #print(chi(b.inverse()) * vec[i])

    return newVec
# Linear time

# Prints the matrix nicely, rounding the actual digits and displays zero as zero
def printMatrix(M):
    decimalPlaces = 4
    print("Printing matrix:")
    for row in M:
        for elem in row:
            temp = 0
            if elem != 0:
                temp = round(elem.real(), decimalPlaces) + round(elem.imag(), decimalPlaces) * I
            print("{:<{mx}}".format(temp, mx=15), end="\t")
        print("")
# Should automate length of mx



# PLAN:
# Take some g
# Compute it as matrix
# find eigenspaces
# take second g
# repeat process until only one eigenspace left

g = G.random_element()
#g = G([[1, 0], [0, 1]])

V = VectorSpace(QQbar, q+1)
H = Hom(V, V)

print(H)
#print(H.domain())
#print(H.domain().basis())

'''
for g in G:
    img = [gAction(g, V([0]*j+[1]+[0]*(q-j))) for j in range(q+1)]
    f = H(img)
    M = f.matrix()
    printMatrix(M)
'''




def eigenVectors(g):
    img = [gAction(g, V([0]*j+[1]+[0]*(q-j))) for j in range(q+1)]
    f = H(img)
    M = f.matrix()
    eigenVecs = M.eigenvectors_left()
    return eigenVecs



memorySet = set()
g = G.random_element()
temp = eigenVectors(g)
for item in temp:
    memorySet.update(item[1])
memory = V.subspace([V(x) for x in memorySet])
print(memory)

count = 0
for g in G:
    temp = eigenVectors(g)
    tempSet = set()
    for item in temp:
        tempSet.update(item[1])
    tempVS = V.subspace([V(x) for x in tempSet])
    if count < 6:
        print("This is tempSet:")
        print(tempSet)
        print(tempVS)
    memory = memory.intersection(tempVS)
    print(memory)
    count = count + 1

print(memory)


# THIS DOESNT WORK AT ALL, NEED TO INTERSECT EIGENSPACES!!!