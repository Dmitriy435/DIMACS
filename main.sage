import time
from multiprocessing import Pool
# Preliminaries

# Field and general linear group
q = 5
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
print(len(goodChars))
print(len(badChars))









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
        if vec[i] == 0:
            continue
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
    #sol = 0
    #for elem in G: # Roughly q^4 elems
        #temp = (gAction(elem, vec1)).dot_product(conjugate(gAction(elem, vec2)))
        #print(temp)
        #sol = sol + temp

    sol = sum([(gAction(elem, vec1)).dot_product(conjugate(gAction(elem, vec2))) for elem in G])
    sol = sol / G.order()
    return round(sol.real(), 5) + round(sol.imag(), 5) * I
# Highly inefficient - q^5 runtime, figure out how to somehow iterate only through cosets???


# Generates the matrix of coefficients (idk actual name of this) given element g
def bigMatrix(g):
    start_time = time.time()

    basisMatrix = matrix(CC, q+1, q+1, 0)
    for i in range(0, q+1):
        for j in range(i, q+1):
            vec1 = V([0]*(q+1))
            vec2 = V([0]*(q+1))
            vec1[i] = 1
            vec2[j] = 1

            x = innerProduct(vec1, vec2)
            basisMatrix[i, j] = x
            basisMatrix[j, i] = conjugate(x)
        print("Row " + str(i) + " of basis matrix")
    print(basisMatrix)

    M = matrix(CC, q+1, q+1, 0)
    for i in range(0, q+1):
        for j in range(0, q+1):
            vec1 = V([0]*(q+1))
            vec2 = V([0]*(q+1))
            vec1[i] = 1
            vec2[j] = 1
            
            realvec1 = gAction(g, vec1)

            M[i, j] = realvec1.dot_product(basisMatrix.column(j))
            #M[i, j] = realvec1.dot_product(vec2)
            #M[i, j] = innerProduct(gAction(g, vec1), vec2)
            print("Filled in entry in row " + str(i) + " and column " + str(j))
    end_time = time.time()
    print('bigMatrix time: %f'%(end_time - start_time))
    return M
# Runs in q^7 time :(
# WARNING - even with q=5, this function takes multiple minutes to run


'''
def oneInnerProduct(vec1, vec2):
    return(innerProduct(gAction(g, vec1), vec2))

def alt_bigMatrix(g):
    start_time = time.time()
    vecs = [V([0]*j+[1]+[0]*(q-j)) for j in range(q+1)]
    veclist = [[vec1, vec2] for vec1 in vecs for vec2 in vecs]
    with Pool() as pool:
        M = pool.starmap(oneInnerProduct, veclist)
    end_time = time.time()
    sol = matrix(CC, q+1, q+1, M)
    print('alt_bigMatrix time: %f'%(end_time - start_time))
    return(sol)
'''


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

vec1 = V.random_element()
vec2 = V.random_element()
print(vec1)
print(vec2)
x = innerProduct(vec1, vec2)
print(x)
#y = innerProduct(gAction(g, vec1), gAction(g, vec2))
y = vec1.dot_product(conjugate(vec2))
print(y)


g = G.random_element()
print(g)
#g = G([[1, 0], [0, 1]])
M = bigMatrix(g)
printMatrix(M)


#g = G.random_element()
#M = bigMatrix(g)
#printMatrix(M)
#M = alt_bigMatrix(g)
#print(M)


#vec1 = V([0]*(q+1))
#vec1[0] = 1

#for elem in G:
    #print(elem)
    #print(gAction(elem, vec1))