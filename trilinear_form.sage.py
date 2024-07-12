

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
Vsmall = VectorSpace(CC, q-_sage_const_1 )
Vbig = VectorSpace(QQbar, q+_sage_const_1 )
H = Hom(Vbig, Vbig)








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
print("Lengths of the good chars and bad chars of B (good chars double counted)")
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
print("This is how many nondecomp chars of L we have (also double counted)")
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
    newVec = Vbig([_sage_const_0 ] * (q+_sage_const_1 ))
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
        newVec = Vsmall([_sage_const_0 ] * (q-_sage_const_1 ))

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
        newVec = Vsmall([_sage_const_0 ] * (q-_sage_const_1 ))

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
    img = [gActionInduced(g, basisVec, chi) for basisVec in Vbig.basis()]
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




chi = badCharsB[_sage_const_1 ]
findGsubspace(chi)

g = G.random_element()

vec1 = Vbig([_sage_const_0 ] * (q+_sage_const_1 ))
vec1[_sage_const_0 ] = _sage_const_1 
vec2 = Vbig([_sage_const_0 ] * (q+_sage_const_1 ))
vec2[_sage_const_1 ] = _sage_const_1 
vec3 = Vbig([_sage_const_0 ] * (q+_sage_const_1 ))
vec3[_sage_const_2 ] = _sage_const_1 

print(gActionInduced(g, vec1, chi))
print(gActionInduced(g, vec2, chi))
print(gActionInduced(g, vec3, chi))
