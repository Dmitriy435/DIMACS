

# This file was *autogenerated* from the file trilinear_forms_cuspidal.sage
from sage.all_cmdline import *   # import sage library

_sage_const_4 = Integer(4); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_5 = Integer(5)# Combine are the representations here and calculate the trilinear form

# define inner product on cuspidal? - seems like normal dot product works


# Field and general linear group
q = _sage_const_4 

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

# Determining the characters of K^x that make up each char of B
'''
test1 = G(MS([[a, 0], [0, 1]]))
test2 = G(MS([[1, 0], [0, a]]))
print("Our generator of K^x is " + str(a))
print("")
for i in range(len(goodCharsB)):
    print("This is character " + str(i) + ", and the following are the two values of chi_1 and chi_2 on the generator:")
    print(goodCharsB[i](test1))
    print(goodCharsB[i](test2))
    print("")
'''






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
        for x in K:
            if x!= _sage_const_0  and conjugate(char1(Lx(MSforL([x])))) != char2(Lx(MSforL([x]))):
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



# Find 2 q+1 and 1 cuspidal whose central char is trivial !!!
C = Combinations(list(range(len(goodCharsB))) * _sage_const_2 , _sage_const_2 )
print("")
count = _sage_const_0 
for chars in C:
    chi1 = goodCharsB[chars[_sage_const_0 ]]
    chi2 = goodCharsB[chars[_sage_const_1 ]]
    for i in range(_sage_const_0 , len(nondecomposableChars)):
        cusp = nondecomposableChars[i]

        centralCharTrivial = True
        for x in K:
            if x != _sage_const_0 :
                d = G(MS([[x, _sage_const_0 ], [_sage_const_0 , x]]))
                prod = chi1(d) * chi2(d) * nu(x, cusp)
                if prod != _sage_const_1 :
                    centralCharTrivial = False
                    break
        if centralCharTrivial:
            print("Triple of two good characters and one cuspidal rep with trivial central character:")
            print(chars, end=", ")
            print(i)
            count = count + _sage_const_1 
print("There are " + str(count) + " combos")









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





# CHANGE CHARACTERS USED HERE!!!

induced1 = goodCharsB[_sage_const_1 ]
induced2 = goodCharsB[_sage_const_2 ]
cusp = nondecomposableChars[_sage_const_5 ]


#vec1 = Vinduced.basis()[0]
#vec2 = Vinduced.basis()[1]
#vec3 = Vcuspidal.basis()[0]


# Evaluates function vec at value g
def evalInduced(g, vec, char):
    rep = toRepresentativeInduced(g)
    index = cosetRepsB.index(rep)
    if vec[index] == _sage_const_0 :
        return _sage_const_0 
    b = g * rep.inverse()
    return vec[index] * char(b)



# Gives the triple product - have to specify the representations inside the function
def tripleProduct(g, v1, v2, v3):
    one = matrixCoeffInduced(g, v1, v1, induced1)
    if one == _sage_const_0 :
        return _sage_const_0 
    two = matrixCoeffInduced(g, v2, v2, induced2)
    if two == _sage_const_0 :
        return _sage_const_0 
    three = matrixCoeffCuspidal(g, v3, v3, cusp)
    sol = one * two * three
    #if sol != 0:
    #    print(g.inverse())
    #print(one * two * three)
    return one * two * three

# Gives trilinear form by avergaing matrix coefficients over the whole group
def trilinearForm(v1, v2, v3):
    l = [tripleProduct(g, v1, v2, v3) for g in G]
    l = [i for i in l if i != _sage_const_0 ]
    print(len(l))
    s = sum(l)
    return s / G.order()

# Given vector of induced rep and element g, returns the Whittaker fn of that element
def whittaker(g, vec, char):
    s = _sage_const_0 
    for z in K:
        inp = G(MS([[_sage_const_0 , _sage_const_1 ], [_sage_const_1 , _sage_const_0 ]])) * G(MS([[_sage_const_1 , z], [_sage_const_0 , _sage_const_1 ]])) * g 
        s = s + evalInduced(inp, vec, char) * (_sage_const_1  / psi(z))
    return s

def whittaker2(g, vec, char):
    s = _sage_const_0 
    for z in K:
        inp = G(MS([[_sage_const_0 , _sage_const_1 ], [_sage_const_1 , _sage_const_0 ]])) * G(MS([[_sage_const_1 , z], [_sage_const_0 , _sage_const_1 ]])) * g 
        s = s + evalInduced(inp, vec, char) * psi(z)
    return s

#Computes RS trilinear form
def RStrilinearForm(v1, v2, v3):
    s = _sage_const_0 
    for g in G:
        val1 = evalInduced(g, v1, induced1)
        if val1 == _sage_const_0 :
            continue
        val2 = whittaker(g, v2, induced2)
        if val2 == _sage_const_0 :
            continue
        s = s + val1 * val2 * whittaker2(g, v3, induced3)

    s = s / G.order()
    return norm(s)
print("")




# Iterates over all basis vectors and computes the trilinear forms of them

for i, j, k in cartesian_product((range(q+_sage_const_1 ), range(q+_sage_const_1 ), range(q-_sage_const_1 ))):
    v1 = Vinduced.basis()[i]
    v2 = Vinduced.basis()[j]
    v3 = Vcuspidal.basis()[k]
    if v1 == v2 or v1 == v3 or v2 == v3:
        print("Not all different basis vecs")
    calc1 = trilinearForm(v1, v2, v3)
    #calc2 = RStrilinearForm(v1, v2, v3)
    #if calc2 < 0.0000001:
    #    calc2 = 0
    print(calc1)
    #print(calc2)
    if calc1 != _sage_const_0 :
        #print(calc2 / calc1)
        print(v1)
        print(v2)
        print(v3)
    print("")

