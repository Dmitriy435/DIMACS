
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

# Creating vector space for induced rep
V = VectorSpace(CC, q+1)






# Takes in character of B and returns if the induced representation is irreducible or not
def isInducedIrreducible(chi):
    induced_chi = chi.induct(G)
    return induced_chi.is_irreducible()


ct = B.character_table()
#print(ct)
chars = []
IIChars = []
for i in range(0, len(B.conjugacy_classes_representatives())):
    if ct[i][0] == 1:
        chars.append(B.character(ct[i]))
        if isInducedIrreducible(B.character(ct[i])):
            IIChars.append(B.character(ct[i]))
chi = IIChars[0]



# The chosen coset representations 
# First q spots counting up from 0 to q-1
# Last spot (index q) is w
cosetReps = []
for i in range(0, q):
    cosetReps.append(G(MS([[1,0],[i,1]])).inverse())
    #print(cosetReps[i])
cosetReps.append(G(MS([[0,1],[1,0]])).inverse())
#print(cosetReps[q])




# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentative(g):
    if g.inverse().matrix()[0][0] == 0:
        return G(MS([[0,1],[1,0]])).inverse()
    else:
        x = (1 / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        return G(MS([[1,0],[x,1]])).inverse()


# Gives the G action result of group element g from G onto vector v from V
# Uses globally defined character chi
# Uses globally defined vector space V
def gAction(g, vec):
    newVec = V([0] * (q+1))
    for i in range(0, q+1):
        newRep = toRepresentative(cosetReps[i] * g.inverse())
        b = cosetReps[i] * g.inverse() * newRep.inverse()
        newIndex = -1

        for j in range(0, q+1):
            if newRep == cosetReps[j]:
                newIndex = j
                break
        
        newVec[newIndex] = chi(b.inverse()) * vec[i]

    return newVec


# Takes naive inner product (built in dot product) and computes the averaged product (invariant under g action on both vectors)
def innerProduct(vec1, vec2):
    sol = 0
    for elem in G.list():
        temp = (gAction(elem, vec1)).dot_product(gAction(elem, vec2))
        sol = sol + temp
    sol = sol / G.order()
    return round(sol.real(), 5) + round(sol.imag(), 5) * I
# Highly inefficient - q^6 runtime, figure out how to somehow iterate only through cosets???


# Generates the matrix of coefficients (idk actual name of this) given element g
def bigMatrix(g):
    M = matrix(CC, q+1, q+1, 0)
    for i in range(0, q+1):
        for j in range(0, q+1):
            vec1 = V([0]*(q+1))
            vec2 = V([0]*(q+1))
            vec1[i] = 1
            vec2[j] = 1

            M[i, j] = innerProduct(gAction(g, vec1), vec2)
            print("Filled in entry in row " + str(i) + " and column " + str(j))
    return M
# Runs in q^8 time :(
# WARNING - even with q=5, this function takes multiple minutes to run










#g = G.random_element()
#vec = V([1, 2, 3, 4 + 5*I, -2, 10])
#print(gAction(g, vec))

#vec1 = V([1, 0, 3, 4])
#vec2 = V([2, 2, 1, 0])
#vec1 = V.random_element()
#vec2 = V.random_element()
#print(vec1)
#print(vec2)
#x = innerProduct(vec1, vec2)
#print(x)
#y = innerProduct(gAction(g, vec1), gAction(g, vec2))
#print(y)

#print(g)
#print(g.inverse())
#print(gAction(g, vec1))


g = G.random_element()
M = bigMatrix(g)
print(M)
