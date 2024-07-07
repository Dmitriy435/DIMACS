# Preliminaries

# Field and general linear group
q = 3
K = GF(q)
G = GL(2, K)

# Borel subgroup
a = K.zeta()
MS = MatrixSpace(K, 2)
gens = [MS([[1,1],[0,1]]),MS([[a,0],[0,1]]), MS([[1,0],[0,a]])]
B = G.subgroup(gens)

# U subgroup
gens = [MS([[1,x],[0,1]]) for x in K]
U = G.subgroup(gens)

# Creating vector space for induced rep
dim = (q+1) * (q-1)^2
V = VectorSpace(QQbar, dim)
H = Hom(V, V)


# The chosen coset representations of U\G
cosetRepsB = []
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetRepsB.append(rep)
rep = G(MS([[0,1],[1,0]]))
cosetRepsB.append(rep)

cosetRepsD = []
for x in K:
    if x != 0:
        for y in K:
            if y != 0:
                rep = G(MS([[1 / x, 0],[0, 1 / y]]))
                cosetRepsD.append(rep)

repToIndex = {}
index = 0

cosetReps = []
for gRep in cosetRepsB:
    for dRep in cosetRepsD:
        rep = dRep * gRep
        cosetReps.append(rep)
        repToIndex[rep] = index
        index = index + 1
print(cosetReps)
print(len(cosetReps))









ct = U.character_table()
#print(ct)

#for chi in ct:
#    print(U.character(chi).induct(G).values())
# Notice this induced rep same for all non-unit chars!!!

chi = U.character(ct[1])
print(chi.values())






# NEED TO CHANGE THIS!!!

# Takes in element g from G and returns the corresponding representative in U\G
def toRepresentative(g):
    gRep = -1
    if g.inverse().matrix()[0][0] == 0:
        gRep = G(MS([[0,1],[1,0]]))
    else:
        x = (K.one() / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        gRep = G(MS([[1,0],[x,1]]))

    if gRep == -1:
        print("AHHHHHHHH")

    b = gRep.inverse() * g.inverse()
    d = G(MS([[b.matrix()[0][0], 0], [0, b.matrix()[1][1]]]))

    u = gRep * d
    return u.inverse()
# Constant time


# Gives the G action result of group element g from G onto vector v from V
# Uses globally defined character chi
# Uses globally defined vector space V
def gAction(g, vec):
    newVec = V([0] * dim)
    for i in range(0, dim):
        if vec[i] == 0:
            continue
        newRep = toRepresentative(cosetReps[i] * g.inverse())
        #print(newRep)
        u = cosetReps[i] * g.inverse() * newRep.inverse()
        
        newIndex = repToIndex[newRep]

        newVec[newIndex] = newVec[newIndex] + chi(u.inverse()) * vec[i]
        #print(chi(u.inverse()) * vec[i])

    return newVec
# Linear time

# Prints the matrix nicely, rounding the actual digits and displays zero as zero
def printMatrix(M):
    decimalPlaces = 4
    print("Printing matrix:")
    for row in M:
        print("[", end="")
        for elem in row:
            temp = 0
            if elem != 0:
                temp = round(elem.real(), decimalPlaces) + round(elem.imag(), decimalPlaces) * I
            print("{:<{mx}}".format(temp, mx=10), end="\t")
        print("]")
# Should automate length of mx

# Finds eigenSpaces of particular g given the character chi for which this representations is induced
def eigenVectors(g):
    img = [gAction(g, basisVec) for basisVec in V.basis()]
    f = H(img)
    M = f.matrix()
    eigenSpaces = M.eigenvectors_left()
    return eigenSpaces

# Goes through eigenvectors, chooses q-1, 1, and q+1 of them to form a set of spaces
def findInvariantSubspacesH(g):

    # List of eigenvectors:
    temp = eigenVectors(g)
    listVecs = []
    for t in temp:
        listVecs = listVecs + t[1]
    #print(len(listVecs))

    sol = set()

    # q-1:
    s = Subsets(listVecs, q-1)
    #print(s.cardinality())
    for gens in s:
        v = V.subspace(gens)
        sol.add(v)

    # q:
    s = Subsets(listVecs, q)
    #print(s.cardinality())
    for gens in s:
        v = V.subspace(gens)
        sol.add(v)

    # q+1:
    s = Subsets(listVecs, q+1)
    #print(s.cardinality())
    for gens in s:
        v = V.subspace(gens)
        sol.add(v)

    return sol
# Very not efficient

# Iterates through all g, finds all invariant subspaces, and finds intersections
def findInvariantSubspaces():
    memorySet = set()
    g = G.random_element()
    s = findInvariantSubspacesH(g)
    memorySet.update(s)

    for elem in G:
        print("Looking at element:")
        print(elem)
        print("Currently memorySet has " + str(len(memorySet)) + " elements")
        if len(memorySet) == q - 1 + (q^2 - q + q^2 - 3*q + 2) / 2:
            break
        
        s = findInvariantSubspacesH(elem)

        tempSet = set()
        for ogSpace in memorySet:
            for newSpace in s:
                t = ogSpace.intersection(newSpace)

                if t.dimension() == q-1 or t.dimension() == q or t.dimension() == q+1:
                    tempSet.add(t)
        #print(tempSet)
        memorySet = tempSet
        
    print(memorySet)
# Stupidly ineffient
# THIS DOES NOT WORK!!!




'''
g = G(MS([[1, 2], [1, 1]]))
print(g)

print(toRepresentative(g))
print(g * toRepresentative(g).inverse())

vec = V.random_element()

print(vec)
print(gAction(g, vec))
'''






g = G(MS([[1, 2], [1, 1]]))
print(g)

#print(eigenVectors(g))

findInvariantSubspaces()