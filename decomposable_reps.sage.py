

# This file was *autogenerated* from the file decomposable_reps.sage
from sage.all_cmdline import *   # import sage library

_sage_const_5 = Integer(5); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_15 = Integer(15)
import numpy as np

# Preliminaries

# Field and general linear group
q = _sage_const_5 
K = GF(q)
G = GL(_sage_const_2 , K)

# Borel subgroup
a = K.zeta() #Generator of K*
MS = MatrixSpace(K, _sage_const_2 )
gens = [MS([[_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ]]),MS([[a,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ]]), MS([[_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,a]])]
B = G.subgroup(gens)

# Creating vector space for induced rep
V = VectorSpace(QQbar, q+_sage_const_1 )
H = Hom(V, V)

# The chosen coset representations 
# First q spots are of the form [[1, 0][-i, 1]] (i counts up from 0 to q-1)
# Last spot (index q) is w
cosetReps = []
repToIndex = {}
index = _sage_const_0 
for x in K:
    rep = G(MS([[_sage_const_1 ,_sage_const_0 ],[-x,_sage_const_1 ]]))
    cosetReps.append(rep)
    repToIndex[rep] = index
    #print(rep)
    index = index + _sage_const_1 
rep = G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]]))
cosetReps.append(rep)
repToIndex[rep] = q
#print(cosetReps[q])










# Takes in character of B and returns if the induced representation is irreducible or not
def isInducedIrreducible(chi):
    induced_chi = chi.induct(G)
    return induced_chi.is_irreducible()

# Producing the lists of good characters and bad characters
ct = B.character_table()
badChars = []
goodChars = []
for i in range(_sage_const_0 , len(B.conjugacy_classes_representatives())):
    if ct[i][_sage_const_0 ] == _sage_const_1 :
        if isInducedIrreducible(B.character(ct[i])):
            goodChars.append(B.character(ct[i]))
        else:
            badChars.append(B.character(ct[i]))
print(len(goodChars))
print(len(badChars))

# NOTICE - the induced reps of the goodChars are each counted twice, as two diff chars map to the same induced rep







# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentative(g):
    if g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ] == _sage_const_0 :
        return G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]])).inverse()
    else:
        x = (K.one() / g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ]) * g.inverse().matrix()[_sage_const_1 ][_sage_const_0 ]
        return G(MS([[_sage_const_1 ,_sage_const_0 ],[x,_sage_const_1 ]])).inverse()
# Constant time


# Gives the G action result of group element g from G onto vector v from V, with rep induced by chi
# Uses globally defined vector space V
def gAction(g, vec, chi):
    newVec = V([_sage_const_0 ] * (q+_sage_const_1 ))
    for i in range(_sage_const_0 , q+_sage_const_1 ):
        if vec[i] == _sage_const_0 :
            continue
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
    decimalPlaces = _sage_const_4 
    print("Printing matrix:")
    for row in M:
        for elem in row:
            temp = _sage_const_0 
            if elem != _sage_const_0 :
                temp = round(elem.real(), decimalPlaces) + round(elem.imag(), decimalPlaces) * I
            print("{:<{mx}}".format(temp, mx=_sage_const_15 ), end="\t")
        print("")
# Should automate length of mx





# Finds eigenSpaces of particular g given the character chi for which this representations is induced
def eigenSpaces(g, chi):
    img = [gAction(g, basisVec, chi) for basisVec in V.basis()]
    f = H(img)
    M = f.matrix()
    eigenSpaces = M.eigenspaces_left()
    return eigenSpaces


# Given character of B chi, finds the 1d (if exists) G-invariant subspace of the induced representation
def findGsubspace(chi):
    memorySet = set()
    g = G.random_element()
    spaces = eigenSpaces(g, chi)
    for t in spaces:
        #print(t[1])
        memorySet.add(t[_sage_const_1 ])
    #print(memorySet)

    for g in G:
        if len(memorySet) == _sage_const_1  and list(memorySet)[_sage_const_0 ].dimension() == _sage_const_0 :
            break

        # ONLY FOR DEALING WITH BAD CHARACTERS!!! MAY CAUSE ERRORS WHEN TESTING THIS ON THE GOOD CHARACTERS!!!
        if len(memorySet) == _sage_const_2 :
            l = []
            for x in memorySet:
                l.append(x.dimension())
            if l == [_sage_const_0 , _sage_const_1 ] or l == [_sage_const_1 , _sage_const_0 ]:
                break
        # Although reduces accuracy, the speed is improved 100 fold

        spaces = eigenSpaces(g, chi)
        gSet = set()
        for t in spaces:
            gSet.add(t[_sage_const_1 ])

        tempSet = set()
        for ogSpace in memorySet:
            for newSpace in gSet:
                t = ogSpace.intersection(newSpace)
                tempSet.add(t)
        #print(tempSet)
        memorySet = tempSet

    if len(memorySet) == _sage_const_1 :
        print("This was a good character! No G-invariant subspace!")
    else:
        print("This is the 1d G-invariant subspace:")
        for item in memorySet:
            if item.dimension()==_sage_const_1 :
                print(item)
                return item
# Runs pretty slowly when bad character - any way to speed this up?
# Could start checking if only 1d subspace left, then just simply check if this remains to be eigenvector for remainding elems

'''
for chi in goodChars:
    findGsubspace(chi)
for chi in badChars:
    findGsubspace(chi)
'''




chi = badChars[_sage_const_1 ]
print(chi.values())
W = findGsubspace(chi)

'''
V2 = V / W
print(V2)
liftMap = V2.lift_map()
quotientMap = V2.quotient_map()

print("")


print(liftMap.image()) # This is the space we want to be working over!
vec = liftMap(V2([5, 1, -2]))
print(vec)
g = G([[1, 2], [2, 2]])
print(g)
print(gAction(g, vec, chi))
print(liftMap(quotientMap(gAction(g, vec, chi)))) # How we compute gAction on this space!
'''


'''
vec = W.basis()[0]
print(vec)

s = set()
for g in G:
    x = tuple(gAction(g, vec, chi))
    s.add(x)
print(s)
'''


print("")


# FOUND THE Q DIMENSIONAL SUBSPACE!!!!
# JUST THE ORTHO COMPLEMENT!!!!

U = W.complement()
print(U)
print(U.basis())
vec = U.basis()[_sage_const_0 ]

gens = []
for g in G:
    v = gAction(g, vec, chi)
    gens.append(v)

test = V.subspace(gens)
print(test)



g = G([[_sage_const_1 , _sage_const_2 ], [_sage_const_1 , _sage_const_1 ]])
vec = U.random_element()
vec2 = W.basis()[_sage_const_0 ]

print(vec)
print(gAction(g, vec, chi))
print(vec.dot_product(conjugate(gAction(g, vec, chi))))
print(vec2)
print(gAction(g, vec2, chi))

