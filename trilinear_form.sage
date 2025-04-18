# THIS FILE IS OBSOLETE - DO NOT USE
from mpmath import *

#################

q = 4

#################




# Field and general linear group
K = GF(q)
G = GL(2, K)
MS = MatrixSpace(K, 2)

L = GF(q^2)

# Borel subgroup
a = K.zeta()
print("This is the generator of K^*:")
print(a)
print()
gens = [MS([[1,1],[0,1]]),MS([[a,0],[0,1]]), MS([[1,0],[0,a]])]
B = G.subgroup(gens)

# U subgroup
gens = [MS([[1,x],[0,1]]) for x in K]
U = G.subgroup(gens)


# Vector subspaces:
Vcuspidal = VectorSpace(QQbar, q-1)
Vinduced = VectorSpace(QQbar, q+1)
H = Hom(Vinduced, Vinduced)








# Representatives
# Use .index() to find index

# For Cuspidal:
basisRepsCuspidal = K.list()
basisRepsCuspidal.remove(0)

#print(basisRepsCuspidal[0])

# Reps of B
cosetRepsB = []
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetRepsB.append(rep)
    print(rep)
rep = G(MS([[0,1],[1,0]]))
print(rep)
cosetRepsB.append(rep)


# Reps of U
cosetRepsD = []
for x in K:
    if x != 0:
        for y in K:
            if y != 0:
                rep = G(MS([[1 / x, 0],[0, 1 / y]]))
                cosetRepsD.append(rep)

cosetRepsU = []
for gRep in cosetRepsB:
    for dRep in cosetRepsD:
        rep = dRep * gRep
        cosetRepsU.append(rep)







# Finding all necessary characters

# Takes in character of B and returns if the induced representation is irreducible or not
def isInducedIrreducible(chi):
    induced_chi = chi.induct(G)
    return induced_chi.is_irreducible()

# Producing the lists of good characters and bad characters of B
ct = B.character_table()
badCharsB = []
goodCharsB = []
for i in range(0, len(B.conjugacy_classes_representatives())):
    if ct[i][0] == 1:
        if isInducedIrreducible(B.character(ct[i])):
            goodCharsB.append(B.character(ct[i]))
        else:
            badCharsB.append(B.character(ct[i]))

for i in range(0, len(goodCharsB)):
    if i >= len(goodCharsB):
        break

    char1 = goodCharsB[i]
    for j in range(i+1, len(goodCharsB)):
        char2 = goodCharsB[j]

        equalsInverse = True
        for x in K:
            for y in K:
                if x == 0 or y == 0:
                    continue
                g = G(MS([[x, 0], [0, y]]))
                w = G(MS([[0, 1], [1, 0]]))
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
C = Combinations(list(range(len(goodCharsB))) * 3, 3)
print("")

for chars in C:
    chi1 = goodCharsB[chars[0]]
    chi2 = goodCharsB[chars[1]]
    chi3 = goodCharsB[chars[2]]

    centralCharTrivial = True
    for x in K:
        if x != 0:
            d = G(MS([[x, 0], [0, x]]))
            prod = chi1(d) * chi2(d) * chi3(d)
            if prod != 1:
                centralCharTrivial = False
                break
    if centralCharTrivial:
        print("Triple of good characters with trivial central character:")
        print(chars)
        

print("")
'''



# Fixing character of U (this choice doesn't matter)
ct = U.character_table()
charU = U.character(ct[1])
print("This is the character of K we are using: ")
print(charU.values())
print("")

# Pass in member of K and get the value of the char of U back
def psi(x):
    return charU(G(MS([[1, x], [0, 1]])))

# Getting nondecomp characters of L^x
Lx = GL(1, L)
MSforL = MatrixSpace(L, 1)
ct2 = Lx.character_table()

# Conjugate of elem in L
def conjugateOfL(l):
    return l^q

nondecomposableChars = []
for i in range(0, len(Lx.conjugacy_classes_representatives())):
    if ct2[i][0] == 1:
        char = Lx.character(ct2[i])
        decomposable = True
        for x in Lx:
            if char(x) != char(conjugateOfL(x)):
                decomposable = False
                break
        if not decomposable:
            nondecomposableChars.append(char)

for i in range(0, len(nondecomposableChars)):
    if i >= len(nondecomposableChars):
        break

    char1 = nondecomposableChars[i]
    for j in range(i+1, len(nondecomposableChars)):
        char2 = nondecomposableChars[j]
        equalsConjugate = True
        for x in K:
            if x!= 0 and conjugate(char1(Lx(MSforL([x])))) != char2(Lx(MSforL([x]))):
                equalsConjugate = False
                break
        if equalsConjugate:
            nondecomposableChars.remove(char2)
            break

print("This is how many nondecomp chars of L we have (double counted reps removed)")
print(len(nondecomposableChars))
print("")

#Giving what these chars actually are:
aL = L.zeta()


# See "https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field_element.html#sage.rings.number_field.number_field_element.NumberFieldElement.coordinates_in_terms_of_powers"


# Char of L - takes in elem in L and an elem in nondecomposableChars
def nu(l, nondecompChar):
    v = Lx(MSforL([l]))
    return nondecompChar(v)


print(aL)
print("This is cuspidal char " + str(2) + "'s value on generator of L^*")
temp = (nu(aL, nondecomposableChars[2])).complex_embedding()
print(temp)
print(arg(temp)*15 / 2 / pi)


print()

determinantG = 0
i = 0
for x in K:
    print("Value of psi on " + str(x))
    print(psi(x))
    i = i + 1
    if i == 3:
        determinantG = x
print()

print(determinantG)
print(determinantG * determinantG * determinantG)
#print(determinantG + determinantG)
print()

print("Finding u:")
print()
for u in L:
    if u * conjugateOfL(u) == determinantG:
        print(u)
        print(u + conjugateOfL(u))
        print(psi(u + conjugateOfL(u)))
        print(arg(nu(u, nondecomposableChars[2]).complex_embedding())* 15 / 2 / pi)
        print()

kappa = e^(2 * pi / 15)
t = -kappa^(-2) + kappa^(-5) - kappa^7 + kappa^4 + kappa
t = -t / 4
print(t)




# For the induced irreps:

# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentativeInduced(g):
    if g.inverse().matrix()[0][0] == 0:
        return G(MS([[0,1],[1,0]])).inverse()
    else:
        x = (K.one() / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        return G(MS([[1,0],[x,1]])).inverse()
# Constant time


# Gives the G action result of group element g from G onto vector v from V, with rep induced by chi
def gActionInduced(g, vec, chi):
    newVec = Vinduced([0] * (q+1))
    for i in range(0, q+1):
        if vec[i] == 0:
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
    if g.matrix()[1, 0] == 0:
        # Easier implementation, g is in B
        newVec = Vcuspidal([0] * (q-1))

        a = g.matrix()[0, 0]
        b = g.matrix()[0, 1]
        d = g.matrix()[1, 1]

        for i in range(0, q-1):
            if vec[i] == 0:
                continue
    
            oldRep = basisRepsCuspidal[i]

            newRep = d * (1 / a) * oldRep
            newIndex = basisRepsCuspidal.index(newRep)
            coefficient = nu(d, nondecompChar) * psi(b * (1 / a) * oldRep)

            newVec[newIndex] = newVec[newIndex] + coefficient * vec[i]
        
        return newVec

    else:
        # Harder longer formula
        newVec = Vcuspidal([0] * (q-1))

        for i in range(0, q-1):
            if vec[i] == 0:
                continue
            
            oldRep = basisRepsCuspidal[i]
            for j in range(0, q-1):
                y = basisRepsCuspidal[j]
                newVec[j] = newVec[j] + coeff(y, oldRep, g, nondecompChar) * vec[i]
        
        return newVec
# Fast


# Helper function for gActionCuspidal - calculates the coefficients when g is not in B
def coeff(y, x, g, nondecompChar):
    a = g.matrix()[0, 0]
    b = g.matrix()[0, 1]
    c = g.matrix()[1, 0]
    d = g.matrix()[1, 1]
    temp = 0

    comparison = y * (1 / x) * (a*d - b*c)
    for u in L:
        if u != 0 and u * conjugateOfL(u) == comparison:
            temp = temp + nu(u, nondecompChar) * psi(- (x / c) * (u + conjugateOfL(u)))
    
    return temp * psi((a * y + d * x) / c) / q
# Fast as well




def findGsubspace(chi):
    sol = Vinduced([1]*(q+1))
    sol[q] = chi(G(MS([[-1, 0], [0, 1]])))
    return Vinduced.subspace([sol])











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

induced1 = goodCharsB[0]
induced2 = goodCharsB[1]
induced3 = goodCharsB[1]

#cusp = nondecomposableChars[0]


vec1 = Vinduced.basis()[0]
vec2 = Vinduced.basis()[1]
vec3 = Vinduced.basis()[2]



# Evaluates function vec at value g
def evalInduced(g, vec, char):
    rep = toRepresentativeInduced(g)
    index = cosetRepsB.index(rep)
    if vec[index] == 0:
        return 0
    b = g * rep.inverse()
    return vec[index] * char(b)


# Gives the triple product - have to specify the representations inside the function
def tripleProduct(g, v1, v2, v3):
    one = matrixCoeffInduced(g, v1, v1, induced1)
    if one == 0:
        return 0
    two = matrixCoeffInduced(g, v2, v2, induced2)
    if two == 0:
        return 0
    three = matrixCoeffInduced(g, v3, v3, induced3)
    sol = one * two * three
    #if sol != 0:
    #    print(g.inverse())
    #print(one * two * three)
    return one * two * three

# Gives trilinear form by avergaing matrix coefficients over the whole group
def trilinearForm(v1, v2, v3):
    l = [tripleProduct(g, v1, v2, v3) for g in G]
    l = [i for i in l if i != 0]
    print(len(l))
    s = sum(l)
    return QQ(s / G.order())

# Given vector of induced rep and element g, returns the Whittaker fn of that element
def whittaker(g, vec, char):
    s = 0
    for z in K:
        inp = G(MS([[0, 1], [1, 0]])) * G(MS([[1, z], [0, 1]])) * g 
        s = s + evalInduced(inp, vec, char) * (1 / psi(z))
    return s

def whittaker2(g, vec, char):
    s = 0
    for z in K:
        inp = G(MS([[0, 1], [1, 0]])) * G(MS([[1, z], [0, 1]])) * g 
        s = s + evalInduced(inp, vec, char) * psi(z)
    return s

#Computes RS trilinear form
def RStrilinearForm(v1, v2, v3):
    s = 0
    for g in G:
        val1 = evalInduced(g, v1, induced1)
        if val1 == 0:
            continue
        val2 = whittaker(g, v2, induced2)
        if val2 == 0:
            continue
        s = s + val1 * val2 * whittaker2(g, v3, induced3)

    s = s / G.order()
    return QQ(norm(s))
print("")


# Testing the normalizing factors of the trilinear forms

'''
s1 = 0
for g in cosetRepsB:
    temp = evalInduced(g, vec1, induced1)
    s1 = s1 + temp * conjugate(temp)
print(s1)

s1 = 0
for g in cosetRepsB:
    temp = evalInduced(g, vec2, induced2)
    s1 = s1 + temp * conjugate(temp)
print(s1)

s1 = 0
for g in cosetRepsB:
    temp = evalInduced(g, vec3, induced3)
    s1 = s1 + temp * conjugate(temp)
print(s1)


s = 0
for g in cosetRepsB:
    temp = whittaker(g, vec2, induced2)
    s = s + temp * conjugate(temp)
print(s)

s2 = 0
for g in cosetRepsB:
    temp = whittaker2(g, vec3, induced3)
    s2 = s2 + temp * conjugate(temp)
print(s2)

print(trilinearForm(vec1, vec2, vec3))
print(RStrilinearForm(vec1, vec2, vec3))
'''



# Iterates over all basis vectors and computes the trilinear forms of them
'''
for i, j, k in cartesian_product((range(q+1), range(q+1), range(q+1))):
    v1 = Vinduced.basis()[i]
    v2 = Vinduced.basis()[j]
    v3 = Vinduced.basis()[k]
    if v1 == v2 or v1 == v3 or v2 == v3:
        print("Not all different basis vecs")
    calc1 = trilinearForm(v1, v2, v3)
    calc2 = RStrilinearForm(v1, v2, v3)
    if calc2 < 0.0000001:
        calc2 = 0
    print(calc1)
    print(calc2)
    if calc1 != 0:
        print(calc2 / calc1)
        print(v1)
        print(v2)
        print(v3)
    print("")
'''