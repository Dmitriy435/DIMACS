

# This file was *autogenerated* from the file trilinear_form.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0)# Combine are the representations here and calculate the trilinear form

# define inner product on cuspidal? - seems like normal dot product works


# Field and general linear group
q = _sage_const_3 
K = GF(q)
G = GL(_sage_const_2 , K)
MS = MatrixSpace(K, _sage_const_2 )

L = GF(q**_sage_const_2 )

# Borel subgroup
a = K.zeta()
gens = [MS([[_sage_const_1 ,_sage_const_1 ],[_sage_const_0 ,_sage_const_1 ]]),MS([[a,_sage_const_0 ],[_sage_const_0 ,_sage_const_1 ]]), MS([[_sage_const_1 ,_sage_const_0 ],[_sage_const_0 ,a]])]
B = G.subgroup(gens)

# U subgroup
gens = [MS([[_sage_const_1 ,x],[_sage_const_0 ,_sage_const_1 ]]) for x in K]
U = G.subgroup(gens)


# Vector subspaces:
Vcuspidal = VectorSpace(QQbar, q-_sage_const_1 )
Vinduced = VectorSpace(QQbar, q+_sage_const_1 )
H = Hom(Vinduced, Vinduced)








# Representatives
# Use .index() to find index

# For Cuspidal:
basisRepsCuspidal = K.list()
basisRepsCuspidal.remove(_sage_const_0 )

# Reps of B
cosetRepsB = []
for x in K:
    rep = G(MS([[_sage_const_1 ,_sage_const_0 ],[-x,_sage_const_1 ]]))
    cosetRepsB.append(rep)
rep = G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]]))
cosetRepsB.append(rep)







# Finding all necessary characters

# Takes in character of B and returns if the induced representation is irreducible or not
def isInducedIrreducible(chi):
    induced_chi = chi.induct(G)
    return induced_chi.is_irreducible()

# Producing the lists of good characters and bad characters of B
ct = B.character_table()
badCharsB = []
goodCharsB = []
for i in range(_sage_const_0 , len(B.conjugacy_classes_representatives())):
    if ct[i][_sage_const_0 ] == _sage_const_1 :
        if isInducedIrreducible(B.character(ct[i])):
            goodCharsB.append(B.character(ct[i]))
        else:
            badCharsB.append(B.character(ct[i]))

for i in range(_sage_const_0 , len(goodCharsB)):
    if i >= len(goodCharsB):
        break

    char1 = goodCharsB[i]
    for j in range(i+_sage_const_1 , len(goodCharsB)):
        char2 = goodCharsB[j]

        equalsInverse = True
        for x in K:
            for y in K:
                if x == _sage_const_0  or y == _sage_const_0 :
                    continue
                g = G(MS([[x, _sage_const_0 ], [_sage_const_0 , y]]))
                w = G(MS([[_sage_const_0 , _sage_const_1 ], [_sage_const_1 , _sage_const_0 ]]))
                if char1(g) != char2(w * g * w):
                    equalsInverse = False
                    break
        if equalsInverse:
            goodCharsB.remove(char2)
            break


print("Lengths of the good chars and bad chars of B (double counted good chars removed)")
print(len(goodCharsB))
print(len(badCharsB))
print("")




# Fixing character of U (this choice doesn't matter)
ct = U.character_table()
charU = U.character(ct[_sage_const_1 ])
print("This is the character of K we are using: ")
print(charU.values())
print("")

# Pass in member of K and get the value of the char of U back
def psi(x):
    return charU(G(MS([[_sage_const_1 , x], [_sage_const_0 , _sage_const_1 ]])))

# Getting nondecomp characters of L^x
Lx = GL(_sage_const_1 , L)
MSforL = MatrixSpace(L, _sage_const_1 )
ct2 = Lx.character_table()

# Conjugate of elem in L
def conjugateOfL(l):
    return l**q

nondecomposableChars = []
for i in range(_sage_const_0 , len(Lx.conjugacy_classes_representatives())):
    if ct2[i][_sage_const_0 ] == _sage_const_1 :
        char = Lx.character(ct2[i])
        decomposable = True
        for x in Lx:
            if char(x) != char(conjugateOfL(x)):
                decomposable = False
                break
        if not decomposable:
            nondecomposableChars.append(char)

for i in range(_sage_const_0 , len(nondecomposableChars)):
    if i >= len(nondecomposableChars):
        break

    char1 = nondecomposableChars[i]
    for j in range(i+_sage_const_1 , len(nondecomposableChars)):
        char2 = nondecomposableChars[j]
        equalsConjugate = True
        for x in Lx:
            if char1(conjugateOfL(x)) != char2(x):
                equalsConjugate = False
                break
        if equalsConjugate:
            nondecomposableChars.remove(char2)
            break

print("This is how many nondecomp chars of L we have (double counted reps removed)")
print(len(nondecomposableChars))
print("")



# Char of L - takes in elem in L and an elem in nondecomposableChars
def nu(l, nondecompChar):
    v = Lx(MSforL([l]))
    return nondecompChar(v)













# For the induced irreps:

# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentativeInduced(g):
    if g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ] == _sage_const_0 :
        return G(MS([[_sage_const_0 ,_sage_const_1 ],[_sage_const_1 ,_sage_const_0 ]])).inverse()
    else:
        x = (K.one() / g.inverse().matrix()[_sage_const_0 ][_sage_const_0 ]) * g.inverse().matrix()[_sage_const_1 ][_sage_const_0 ]
        return G(MS([[_sage_const_1 ,_sage_const_0 ],[x,_sage_const_1 ]])).inverse()
# Constant time


# Gives the G action result of group element g from G onto vector v from V, with rep induced by chi
def gActionInduced(g, vec, chi):
    newVec = Vinduced([_sage_const_0 ] * (q+_sage_const_1 ))
    for i in range(_sage_const_0 , q+_sage_const_1 ):
        if vec[i] == _sage_const_0 :
            continue
        newRep = toRepresentativeInduced(cosetRepsB[i] * g.inverse())
        #print(newRep)
        b = cosetRepsB[i] * g.inverse() * newRep.inverse()
        
        newIndex = cosetRepsB.index(newRep)

        newVec[newIndex] = newVec[newIndex] + chi(b.inverse()) * vec[i]
        #print(chi(b.inverse()) * vec[i])

    return newVec
# Linear time




# Group action as described pg 40 in P-S, given a nondecomp character of L
def gActionCuspidal(g, vec, nondecompChar):
    if g.matrix()[_sage_const_1 , _sage_const_0 ] == _sage_const_0 :
        # Easier implementation, g is in B
        newVec = Vcuspidal([_sage_const_0 ] * (q-_sage_const_1 ))

        a = g.matrix()[_sage_const_0 , _sage_const_0 ]
        b = g.matrix()[_sage_const_0 , _sage_const_1 ]
        d = g.matrix()[_sage_const_1 , _sage_const_1 ]

        for i in range(_sage_const_0 , q-_sage_const_1 ):
            if vec[i] == _sage_const_0 :
                continue
    
            oldRep = basisRepsCuspidal[i]

            newRep = d * (_sage_const_1  / a) * oldRep
            newIndex = basisRepsCuspidal.index(newRep)
            coefficient = nu(d, nondecompChar) * psi(b * (_sage_const_1  / a) * oldRep)

            newVec[newIndex] = newVec[newIndex] + coefficient * vec[i]
        
        return newVec

    else:
        # Harder longer formula
        newVec = Vcuspidal([_sage_const_0 ] * (q-_sage_const_1 ))

        for i in range(_sage_const_0 , q-_sage_const_1 ):
            if vec[i] == _sage_const_0 :
                continue
            
            oldRep = basisRepsCuspidal[i]
            for j in range(_sage_const_0 , q-_sage_const_1 ):
                y = basisRepsCuspidal[j]
                newVec[j] = coeff(y, oldRep, g, nondecompChar) * vec[i]
        
        return newVec
# Fast


# Helper function for gAction - calculates the coefficients when g is not in B
def coeff(y, x, g, nondecompChar):
    a = g.matrix()[_sage_const_0 , _sage_const_0 ]
    b = g.matrix()[_sage_const_0 , _sage_const_1 ]
    c = g.matrix()[_sage_const_1 , _sage_const_0 ]
    d = g.matrix()[_sage_const_1 , _sage_const_1 ]
    temp = _sage_const_0 

    comparison = y * (_sage_const_1  / x) * (a*d - b*c)
    for u in L:
        if u != _sage_const_0  and u * conjugateOfL(u) == comparison:
            temp = temp + nu(u, nondecompChar) * psi(- (x / c) * (u + conjugateOfL(u)))
    
    return temp * psi((a * y + d * x) / c) / q
# Fast as well







# Finds eigenSpaces of particular g given the character chi for which this representations is induced
def eigenSpaces(g, chi):
    img = [gActionInduced(g, basisVec, chi) for basisVec in Vinduced.basis()]
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









# Matrix coeff of cuspidal
def matrixCoeffCuspidal(g, vec1, vec2, nondecompChar):
    v = gActionCuspidal(g, vec1, nondecompChar)
    return v.dot_product(conjugate(vec2))

# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    v = gActionInduced(g, vec1, chi)
    return v.dot_product(conjugate(vec2))
# Do actual sum manually



# Gives the triple product - have to specify the representations and the vectors inside the function
def tripleProduct(g):
    cusp = nondecomposableChars[_sage_const_0 ]
    induced1 = goodCharsB[_sage_const_0 ]
    induced2 = badCharsB[_sage_const_0 ]

    Vtiny = findGsubspace(induced2)
    VtinyComplement = Vtiny.complement()

    one = matrixCoeffCuspidal(g, Vcuspidal.random_element(), Vcuspidal.random_element(), cusp)
    two = matrixCoeffInduced(g, Vinduced.random_element(), Vinduced.random_element(), induced1)
    three = matrixCoeffInduced(g, Vtiny.random_element(), Vtiny.random_element(), induced2)

    return one * two * three

'''
print("")
g = G.random_element()
print(g)

print(tripleProduct(g))
'''














# Trying to implement basically the tensor product of reps



cusp = nondecomposableChars[_sage_const_0 ]
induced1 = goodCharsB[_sage_const_0 ]
induced2 = goodCharsB[_sage_const_0 ]

#Vtiny = findGsubspace(induced2)
#Vinduced2 = Vtiny.complement()




print("")

from sage.modules.tensor_operations import VectorCollection, TensorOperation

temp1 = VectorCollection(Vcuspidal.basis(), QQbar, q-_sage_const_1 )
temp2 = VectorCollection(Vinduced.basis(), QQbar, q+_sage_const_1 )
temp3 = VectorCollection(Vinduced.basis(), QQbar, q+_sage_const_1 )


Vtest = TensorOperation([temp1, temp2, temp3])
homomorphismsTripleProduct = Hom(Vtest, Vtest)
print(Vtest)
#print(Vtest.basis())



# WARNING! For now, the representations used are hard coded into these functions

def gActionTripleProduct(g, vec1, vec2, vec3):
    v1 = gActionCuspidal(g, vec1, cusp)
    v2 = gActionInduced(g, vec2, induced1)
    v3 = gActionInduced(g, vec3, induced2)


    newVec = Vtest([_sage_const_0 ]*(q-_sage_const_1 )*(q+_sage_const_1 )*(q+_sage_const_1 ))

    for i, j, k in cartesian_product((range(q-_sage_const_1 ), range(q+_sage_const_1 ), range(q+_sage_const_1 ))):
        #print(i)
        #print(j)
        #print(temp1.vectors()[i])
        #print(temp2.vectors()[j])
        #print("")
        newVec[Vtest.index_map(i, j, k)] = v1[i] * v2[j] * v3[k]



    return newVec

vec1 = Vcuspidal.random_element()
vec2 = Vinduced.random_element()
vec3 = Vinduced.random_element()
g = G.random_element()

#print(gActionTripleProduct(g, vec1, vec2, vec3))



def tripleProductInvariantSpaces(g):
    img = []

    for vec1 in Vcuspidal.basis():
        #print("")
        for vec2 in Vinduced.basis():
            for vec3 in Vinduced.basis():
                img.append(gActionTripleProduct(g, vec1, vec2, vec3))
    
    f = homomorphismsTripleProduct(img)
    M = f.matrix()
    eigenSpaces = M.eigenspaces_left()
    return eigenSpaces


def findTripleProductInvariantSubspace():
    memorySet = set()
    g = G.random_element()
    spaces = tripleProductInvariantSpaces(g)
    for t in spaces:
        #print(t[1])
        memorySet.add(t[_sage_const_1 ])
    print(len(memorySet))

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

        spaces = tripleProductInvariantSpaces(g)
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
        print(len(memorySet))
        print(g)

    if len(memorySet) == _sage_const_1 :
        print("No G-invariant subspace, uh-oh")
    else:
        print("This is the 1d G-invariant subspace:")
        for item in memorySet:
            if item.dimension()==_sage_const_1 :
                print(item)
                return item


findTripleProductInvariantSubspace()

