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
#print(cosetReps)
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



g = G(MS([[1, 2], [1, 1]]))
print(g)

print(toRepresentative(g))
print(g * toRepresentative(g).inverse())

vec = V.random_element()

print(vec)
print(gAction(g, vec))