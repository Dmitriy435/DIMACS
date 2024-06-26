
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
print(chi.values())




# The chosen coset representations 
# First q spots counting up from 0 to q-1 (not necessarily?)
# Last spot (index q) is w
cosetReps = []
repToIndex = {}
index = 0
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetReps.append(rep)
    repToIndex[rep] = index
    print(rep)
    index = index + 1
rep = G(MS([[0,1],[1,0]]))
cosetReps.append(rep)
repToIndex[rep] = q
print(cosetReps[q])




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


# Takes naive inner product (built in dot product) and computes the averaged product (invariant under g action on both vectors)
def innerProduct(vec1, vec2):
    sol = 0
    for elem in G: # Roughly q^4 elems
        temp = (gAction(elem, vec1)).dot_product(conjugate(gAction(elem, vec2)))
        #print(temp)
        sol = sol + temp
    sol = sol / G.order()
    return round(sol.real(), 5) + round(sol.imag(), 5) * I
# Highly inefficient - q^5 runtime, figure out how to somehow iterate only through cosets???


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
# Runs in q^7 time :(
# WARNING - even with q=5, this function takes multiple minutes to run


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







#vec1 = V([1, 0, 0, 0, 0])

#t1 = gAction(elem, vec1)
#print(t1)
#temp = (t1).dot_product(conjugate(t1))
#print(temp)
#vec2 = V([2, 2, 1, 0, 0])

#vec1 = V.random_element()
#vec2 = V.random_element()
#print(vec1)
#print(vec2)
#x = innerProduct(vec1, vec1)
#print(x)
#y = innerProduct(gAction(g, vec1), gAction(g, vec2))
#print(y)


g = G.random_element()
M = bigMatrix(g)
printMatrix(M)