
# Preliminaries

q = 3
K = GF(q)

G = GL(2, K)

a = K.zeta()
MS = MatrixSpace(K, 2)
gens = [MS([[1,1],[0,1]]),MS([[a,0],[0,1]]), MS([[1,0],[0,a]])]
B = G.subgroup(gens)


ct = B.character_table()
print(ct[4])
chi = B.character(ct[5])




# Creating vector space for induced rep

V = VectorSpace(CC, q+1)



# First q spots counting up from 0 to q-1
# Last spot (index q) is w
cosetReps = []
for i in range(0, q):
    cosetReps.append(G(MS([[1,0],[i,1]])).inverse())
cosetReps.append(G(MS([[0,1],[1,0]])).inverse())




# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentative(g):
    if g.inverse().matrix()[0][0] == 0:
        return G(MS([[0,1],[1,0]])).inverse()
    else:
        x = (1 / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        return G(MS([[1,0],[x,1]])).inverse()


# Gives the G action result of group element g from G onto vector v from V - VERIFY THIS, MAKE SURE IT IS CORRECT
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


# Takes naive inner product (built in dot product) and computes the averaged dot product (invariant under g action on both vectors)
def innerProduct(vec1, vec2):
    sol = 0
    for elem in G.list():
        temp = (gAction(elem, vec1)).dot_product(gAction(elem, vec2))
        #print(elem)
        #print(temp)
        sol = sol + temp
    sol = sol / G.order()
    return round(sol.real(), 5) + round(sol.imag(), 5) * I
# Highly inefficient - q^6 runtime, figure out how to somehow iterate only through cosets???




#g = G.random_element()
#vec = V([1, 2, 3, 4 + 5*I, -2, 10])
#print(gAction(g, vec))

#vec1 = V([1, 0, 3, 4])
#vec2 = V([2, 2, 1, 0])
vec1 = V.random_element()
vec2 = V.random_element()
print(vec1)
print(vec2)
x = innerProduct(vec1, vec2)
print(x)

















# Playing with characters - clean up

#ct = B.character_table()
#print(ct)
#chi = B.character(ct[1])

#b = B.random_element()
# b
#print(chi(b))
#print(b)

# Iterate through character table for the 1-d reps

# induced_chi = chi.induct(G)
# induced_chi
# induced_chi.values()
# induced_chi(G.random_element())
# induced_chi.is_irreducible()
# chi.values() == chi.restrict(B).values()