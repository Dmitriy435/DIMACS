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


# L* inclusion into G
# NOTICE DEPENDS ON CHOICE OF NON SQUARE D!!!

# Write function to determine d - use Legendre symbols!!!!!!!!!!!!! - this only works for certain cases when q=p
d = -1 #q=3
#d = 3 #q=5

gens = []
for x in K:
    for y in K:
        gens.append(MS([[x,d*y],[y,x]]))
gens.remove(MS([[0,0],[0,0]]))
L = G.subgroup(gens)


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
#print(cosetReps)
#print(len(cosetReps))

'''
gens = []
for u in U:
    for g in cosetReps:
        gens.append(u * g)
print(type(gens[0]))
G2 = MatrixGroup(gens)

print(type(list(G2)[0]))
print(type(list(G)[0]))
print(list(G2) == list(G))
'''






# Base character of U

ct = U.character_table()
#print(ct)

#for chi in ct:
#    print(U.character(chi).induct(G).values())
# Notice this induced rep same for all non-unit chars!!!

chi = U.character(ct[1])
print(chi.values())




# Characters of L* - correspond to cuspidal reps

ct2 = L.character_table()

#print(L.order())
#print(ct2)





# Finds non decomposable characters of L*


# Takes a in L, and gives back conjugate
def conjugateInL(a):
    m = copy(a.matrix())
    m[0, 1] = -m[0, 1]
    m[1, 0] = -m[1, 0]
    return G(m)

charsOfL = []
for i in range(0, len(L.conjugacy_classes_representatives())):
    if ct2[i][0] == 1:
        charsOfL.append(L.character(ct2[i]))

nondecomposableChars = []
for char in charsOfL:
    decomposable = True
    for x in L:
        if char(x) != char(conjugateInL(x)):
            decomposable = False
            break
    if not decomposable:
        nondecomposableChars.append(char)

#print("These are nondecomp chars: ")
#for x in nondecomposableChars:
#    print(x.values())
















# REALLY check that this is correct, could be messing everything up so far!!!!


# Takes in element g from G and returns the corresponding representative in U\G
def toRepresentative(g):
    gRep = -1 # This is g_sigma inverse!
    if g.inverse().matrix()[0][0] == 0:
        gRep = G(MS([[0,1],[1,0]]))
    else:
        x = (K.one() / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        gRep = G(MS([[1,0],[x,1]]))

    if gRep == -1:
        print("AHHHHHHHH")

    b = g * gRep
    d = G(MS([[b.inverse().matrix()[0][0], 0], [0, b.inverse().matrix()[1][1]]])).inverse()

    return d * gRep.inverse()
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

        #print(u)
        
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














'''
g = G(MS([[1, 2], [1, 1]]))
print(g)

print(toRepresentative(g))
print(g * toRepresentative(g).inverse())

vec = V.random_element()

print(vec)
print(gAction(g, vec))
'''





nu = nondecomposableChars[0]
print(nu.values())
print("")

vec = V([0] * dim)
vec[0] = 1
#vec[4] = 1
#vec[8] = 1
#vec[12] = 1


'''
g = G.random_element()
print(g)
vec2 = gAction(g, vec)
print(vec2)
vec3 = gAction(g, vec2)
print(vec3)
vec4 = gAction(g, vec3)
print(vec4)

vec = vec + vec2 + vec3 + vec4
'''

'''
for delta in K:
    if delta != 0:
        g = G(MS([[delta, 0], [0, delta]]))
        vec2 = gAction(g, vec)
        print(vec2)
'''





for delta in K:
    if delta != 0:
        g = G(MS([[delta, 0], [0, delta]]))
        vec2 = gAction(g, vec)
        vec = vec + (1 / nu(g)) * vec2
vec = vec / (q-1)


g = G(MS([[2, 0], [0, 2]]))
print("These should be the same: ")
print(gAction(g, vec))
print(nu(g) * vec)
print(gAction(g, vec) == nu(g) * vec)




print("Now testing the subgroup generated by this vector: ")

gens = [gAction(g, vec) for g in G]
#for x in gens:
#    print(x)
#    print("")
W = V.subspace(gens)
print(W)





#for i in range(0, dim):
#    if v2[i] != 0:
#        print(cosetReps[i])
