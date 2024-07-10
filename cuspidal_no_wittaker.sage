# Preliminaries

# Field and general linear group
q = 3
K = GF(q)
#Ltest = GF(q^2).over(K)
L = GF(q^2)


G = GL(2, K)
MS = MatrixSpace(K, 2)


# U subgroup
gens = [MS([[1,x],[0,1]]) for x in K]
U = G.subgroup(gens)


# Creating vector space
dim = q-1
V = VectorSpace(QQbar, dim)
#H = Hom(V, V)
# This vector space of functions K^x to C, basis is fns that take value of 1 on one element and zero elsewehere



# Fix basis:
basisReps = K.list()
basisReps.remove(0)
# Use .index() to find index





# Getting psi 
ct = U.character_table()
charU = U.character(ct[1])
print("This is the character of K we are using: ")
print(charU.values())
print("")

def psi(x):
    return charU(G(MS([[1, x], [0, 1]])))



# Getting nondecomp characters of L^x
H = Hom(K, L)
inclusionMap = H[0]

Lx = GL(1, L)
MSforL = MatrixSpace(L, 1)

ct2 = Lx.character_table()
#print(ct2)

def conjugateOfL(l):
    return l^q

charsOfL = []
for i in range(0, len(Lx.conjugacy_classes_representatives())):
    if ct2[i][0] == 1:
        charsOfL.append(Lx.character(ct2[i]))

nondecomposableChars = []
for char in charsOfL:
    decomposable = True
    for x in Lx:
        if char(x) != char(conjugateOfL(x)):
            decomposable = False
            break
    if not decomposable:
        nondecomposableChars.append(char)


# NEED TO REMOVE CONJUAGETS!
# Currently double counting the cuspidal reps



print("This is the non-decomposable character of L^x we are using: ")
print(nondecomposableChars[0].values())
print("")


# For now, fix nondecomp character:
def nu(l):
    m = nondecomposableChars[0]
    v = Lx(MSforL([l]))
    return m(v)










# Group action as described pg 40 in P-S
def gAction(g, vec):
    if g.matrix()[1, 0] == 0:
        # Easier implementation, g is in B
        newVec = V([0] * dim)

        a = g.matrix()[0, 0]
        b = g.matrix()[0, 1]
        d = g.matrix()[1, 1]

        for i in range(0, dim):
            if vec[i] == 0:
                continue
    
            oldRep = basisReps[i]

            newRep = d * (1 / a) * oldRep
            newIndex = basisReps.index(newRep)
            coefficient = nu(d) * psi(b * (1 / a) * oldRep)

            newVec[newIndex] = newVec[newIndex] + coefficient * vec[i]
        
        return newVec

    else:
        # Harder longer formula
        newVec = V([0] * dim)

        for i in range(0, dim):
            if vec[i] == 0:
                continue
            
            oldRep = basisReps[i]
            for j in range(0, dim):
                y = basisReps[j]
                newVec[j] = coeff(y, oldRep, g) * vec[i]
        
        return newVec
# Shouldn't be too bad efficiency wise?


# Helper function for gAction - calculates the coefficients when g is not in B
def coeff(y, x, g):
    a = g.matrix()[0, 0]
    b = g.matrix()[0, 1]
    c = g.matrix()[1, 0]
    d = g.matrix()[1, 1]
    temp = 0

    comparison = y * (1 / x) * (a*d - b*c)
    for u in L:
        if u != 0 and u * conjugateOfL(u) == comparison:
            temp = temp + nu(u) * psi(- (x / c) * (u + conjugateOfL(u)))
    
    return temp * psi((a * y + d * x) / c) / q
# Pretty fast



g = G.random_element()
#g = G(MS([[1, 2], [1, 1]]))
print(g)
vec = V([0]*dim)
vec[0] = 1

print(gAction(g, vec))
