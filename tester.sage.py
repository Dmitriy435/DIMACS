

# This file was *autogenerated* from the file tester.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_5 = Integer(5)
# Preliminaries

q = _sage_const_3 
K = GF(q)

G = GL(_sage_const_2 , K)

a = K.zeta()
MS = MatrixSpace(K, _sage_const_2 )
gens = [MS([[_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ]]),MS([[a,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ]]), MS([[_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,a]])]
B = G.subgroup(gens)


ct = B.character_table()
print(ct[_sage_const_4 ])
chi = B.character(ct[_sage_const_5 ])




# Creating vector space for induced rep

V = VectorSpace(CC, q+_sage_const_1 )



# First q spots counting up from 0 to q-1
# Last spot (index q) is w
cosetReps = []
for i in range(_sage_const_0 , q):
    cosetReps.append(G(MS([[_sage_const_1 ,_sage_const_0 ],[i,_sage_const_1 ]])).inverse())
cosetReps.append(G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]])).inverse())




# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentative(g):
    if g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ] == _sage_const_0 :
        return G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]])).inverse()
    else:
        x = (_sage_const_1  / g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ]) * g.inverse().matrix()[_sage_const_1 ][_sage_const_0 ]
        return G(MS([[_sage_const_1 ,_sage_const_0 ],[x,_sage_const_1 ]])).inverse()


# Gives the G action result of group element g from G onto vector v from V - VERIFY THIS, MAKE SURE IT IS CORRECT
def Gaction(g, vec):
    newVec = V([_sage_const_0 ] * (q+_sage_const_1 ))
    for i in range(_sage_const_0 , q+_sage_const_1 ):
        newRep = toRepresentative(cosetReps[i] * g.inverse())
        b = cosetReps[i] * g.inverse() * newRep.inverse()
        newIndex = -_sage_const_1 

        for j in range(_sage_const_0 , q+_sage_const_1 ):
            if newRep == cosetReps[j]:
                newIndex = j
                break
        
        newVec[newIndex] = chi(b.inverse()) * vec[i]

    return newVec


# Takes naive inner product (built in dot product) and computes the averaged dot product (invariant under g action on both vectors)
def innerProduct(vec1, vec2):
    sol = _sage_const_0 
    for elem in G.list():
        temp = (Gaction(elem, vec1)).dot_product(Gaction(elem, vec2))
        #print(elem)
        #print(temp)
        sol = sol + temp
    sol = sol / G.order()
    return round(sol.real(), _sage_const_5 ) + round(sol.imag(), _sage_const_5 ) * I
# Highly inefficient - q^6 runtime




#g = G.random_element()
#vec = V([1, 2, 3, 4 + 5*I, -2, 10])
#print(Gaction(g, vec))

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
