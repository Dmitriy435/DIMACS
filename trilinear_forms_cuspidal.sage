import csv
from mpmath import *

################# CHANGE HERE!

q = 5


# Set how many induced and cuspidal we want - must add up to 3


numCus = 1 #This is at least 1 for this file
numInd = 3 - numCus

#################





# Field and general linear group
K = GF(q)
G = GL(2, K)
MS = MatrixSpace(K, 2)

L = GF(q^2)

# Borel subgroup
a = K.zeta()
gens = [MS([[1,1],[0,1]]),MS([[a,0],[0,1]]), MS([[1,0],[0,a]])]
B = G.subgroup(gens)

# U subgroup
gens = [MS([[1,x],[0,1]]) for x in K]
U = G.subgroup(gens)


# Vector subspaces:
Vcuspidal = VectorSpace(QQbar, q-1)
Vinduced = VectorSpace(QQbar, q+1)

# Causes floating point errors, needed for the identification package
#Vcuspidal = VectorSpace(ComplexField(170), q-1)
#Vinduced = VectorSpace(ComplexField(170), q+1)

H = Hom(Vinduced, Vinduced)




# Representatives
# Use .index() to find index

# For Cuspidal:
basisRepsCuspidal = K.list()
basisRepsCuspidal.remove(0)


# Reps of B
cosetRepsB = []
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetRepsB.append(rep)
rep = G(MS([[0,1],[1,0]]))
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



# Char of L - takes in elem in L and an elem in nondecomposableChars
def nu(l, nondecompChar):
    v = Lx(MSforL([l]))
    return nondecompChar(v)

'''
def j(u):
    sum = 0
    for l in Lx:
        if (l^(q+1)).matrix()[0,0] == u:
            sum = sum + psi((l).matrix()[0, 0] + (l^q).matrix()[0, 0]) * nu(l.matrix()[0,0], nondecomposableChars[0])
    return - sum / q

x = 1
y = 3

s = 0
for v in K:
    if v != 0:
        s = s + j(x * v) * j(y * v) * psi(v) / nu(v, nondecomposableChars[0])
print("The left hand side is " + str(s))

s2 = nu(-1, nondecomposableChars[0]) * psi(-x-y) * j(x*y)
print("The right hand side is " + str(s2))
'''




# Determining the characters of K^x that make up each char of B
def decomposeChars():
    test1 = G(MS([[a, 0], [0, 1]]))
    test2 = G(MS([[1, 0], [0, a]]))
    print("Our generator of K^x is " + str(a))
    print("")
    if numInd != 0:
        for i in range(len(goodCharsB)):
            print("This is good character " + str(i) + ", and the following are the two values of chi_1 and chi_2 on the generator:")
            print(goodCharsB[i](test1))
            print(goodCharsB[i](test2))
            print("")
        print("")
    if numCus != 0:
        for i in range(len(nondecomposableChars)):
            print("This is cuspidal character " + str(i) + ", and the following is the value of it on the generator:")
            print(nu(a, nondecomposableChars[i]))
            print("")
        print("")
    print("")


# Prints all valid combinations as specificed by start of how many cusp
def validCombinations():
    combos = []
    if numCus == 1:
        C = Combinations(list(range(len(goodCharsB))) * 2, 2)
        count = 0
        for chars in C:
            chi1 = goodCharsB[chars[0]]
            chi2 = goodCharsB[chars[1]]
            for i in range(0, len(nondecomposableChars)):
                cusp = nondecomposableChars[i]

                centralCharTrivial = True
                for x in K:
                    if x != 0:
                        d = G(MS([[x, 0], [0, x]]))
                        prod = chi1(d) * chi2(d) * nu(x, cusp)
                        if prod != 1:
                            centralCharTrivial = False
                            break
                if centralCharTrivial:
                    print("Triple of two good characters and one cuspidal rep with trivial central character:")
                    print(chars, end=", ")
                    print(i)
                    count = count + 1
                    combos.append((chars[0], chars[1], i))
        print("There are " + str(count) + " combos")
    if numCus == 2:
        C = Combinations(list(range(len(nondecomposableChars))) * 2, 2)
        count = 0
        for chars in C:
            chi1 = nondecomposableChars[chars[0]]
            chi2 = nondecomposableChars[chars[1]]

            for i in range(0, len(goodCharsB)):
                chi3 = goodCharsB[i]

                centralCharTrivial = True
                for x in K:
                    if x != 0:
                        d = G(MS([[x, 0], [0, x]]))
                        prod = nu(x, chi1) * nu(x, chi2) * chi3(d)
                        if prod != 1:
                            centralCharTrivial = False
                            break
                if centralCharTrivial:
                    print("Triple of one good character and two cuspidal reps with trivial central character:")
                    print(i, end=", ")
                    print(chars)
                    count = count + 1
                    combos.append((i, chars[0], chars[1]))
        print("There are " + str(count) + " combos")
    if numCus == 3:
        C = Combinations(list(range(len(nondecomposableChars))) * 3, 3)
        count = 0
        for chars in C:
            chi1 = nondecomposableChars[chars[0]]
            chi2 = nondecomposableChars[chars[1]]
            chi3 = nondecomposableChars[chars[2]]

            centralCharTrivial = True
            for x in K:
                if x != 0:
                    prod = nu(x, chi1) * nu(x, chi2) * nu(x, chi3)
                    if prod != 1:
                        centralCharTrivial = False
                        break
            if centralCharTrivial:
                print("Triple of three cuspidal reps with trivial central character:")
                print(chars)
                count = count + 1
                combos.append((chars[0], chars[1], chars[2]))
        print("There are " + str(count) + " combos")
    print("\n")
    return combos









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
                if i == j: # TEMPORARY DELETE AFTERWARDS!!!!! MAYBE MORE EFFICIENT THIS WAY??? Only for purposes of dot product with itself
                    y = basisRepsCuspidal[j]
                    newVec[j] = newVec[j] + coeff(y, oldRep, g, nondecompChar) * vec[i]     
        return newVec
# Fast


# Helper function for gAction - calculates the coefficients when g is not in B
def coeff(y, x, g, nondecompChar):
    a = g.matrix()[0, 0]
    b = g.matrix()[0, 1]
    c = g.matrix()[1, 0]
    d = g.matrix()[1, 1]
    temp = 0

    comparison = y * (1 / x) * (a*d - b*c)
    for u in L:
        if u * conjugateOfL(u) == comparison:
            temp = temp + nu(u, nondecompChar) * psi(- (x / c) * (u + conjugateOfL(u)))
    return - temp * psi((a * y + d * x) / c) / q # Negative sign has been added
# Fast as well








# Evaluates function vec at value g
def evalInduced(g, vec, char):
    rep = toRepresentativeInduced(g)
    index = cosetRepsB.index(rep)
    if vec[index] == 0:
        return 0
    b = g * rep.inverse()
    return vec[index] * char(b)
# Relatively fast


# Matrix coeff of cuspidal
def matrixCoeffCuspidal(g, vec1, vec2, nondecompChar):
    v = gActionCuspidal(g, vec1, nondecompChar)
    return v.dot_product(conjugate(vec2))
# Fast


# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    v = gActionInduced(g, vec1, chi)
    return v.dot_product(conjugate(vec2))
# Fast


# Gives the triple product - have to specify the representations inside the function
def tripleProduct(g, v1, v2, v3):
    if char1 in goodCharsB:
        one = matrixCoeffInduced(g, v1, v1, char1)
    else:
        one = matrixCoeffCuspidal(g, v1, v1, char1)
    if one == 0:
        return 0

    if char2 in goodCharsB:
        two = matrixCoeffInduced(g, v2, v2, char2)
    else:
        two = matrixCoeffCuspidal(g, v2, v2, char2)
    if two == 0:
        return 0
    
    if char3 in goodCharsB:
        three = matrixCoeffInduced(g, v3, v3, char3)
    else:
        three = matrixCoeffCuspidal(g, v3, v3, char3)

    return one * two * three

# Gives trilinear form by avergaing matrix coefficients over the whole group
def trilinearForm(v1, v2, v3):
    l = [tripleProduct(g, v1, v2, v3) for g in G]
    l = [i for i in l if i != 0]
    print(len(l))
    s = sum(l)
    sol = (s / G.order())
    if sol.imag() < 0.000001:
        sol = sol.real()
    try:
        return QQ(sol)
    except:
        return sol
# We don't know how to find Whittaker fns of cuspidal reps, only matrix coeffs


# Iterates over all basis vectors and computes the trilinear forms of them
def triformBasisVecs():
    sol = []
    for i, j, k in cartesian_product((range(V1.dimension()), range(V2.dimension()), range(V3.dimension()))):
        v1 = V1.basis()[i]
        v2 = V2.basis()[j]
        v3 = V3.basis()[k]
        if v1 == v2 or v1 == v3 or v2 == v3:
            print("Not all different basis vecs")
        calc1 = trilinearForm(v1, v2, v3)
        if calc1 < 0.0000001:
            calc1 = 0
        s = calc1
        #if q == 4:
        #    s = identify(calc1, ['sqrt(5)'])
        #if q == 5:
        #    s = identify(calc1, ['sqrt(3)', 'sqrt(2)', 'sqrt(6)'])

        print(s)
        sol.append(s)
        #if (s == "0"):
        #    print(s)
        #    sol.append(s)
        #else:
        #    print(s[1:-1])
        #    sol.append(s[1:-1])
        #if calc1 != 0:
        #    print(v1)
        #    print(v2)
        #    print(v3)
        print("")
    return sol
# How to do Whittaker?

def createHeader():
    header = ["Representation"]
    dim1, dim2, dim3 = q+1, q+1, q+1

    if numCus == 1:
        dim3 = q-1
    elif numCus == 2:
        dim2, dim3 = q-1, q-1
    elif numCus == 3:
        dim1, dim2, dim3 = q-1, q-1, q-1
    
    for i in range(dim1):
            for j in range(dim2):
                for k in range(dim3):
                    header.append((i, j, k))
    print(header)
    return header



V1 = Vinduced if numCus < 3 else Vcuspidal
V2 = Vinduced if numCus < 2 else Vcuspidal
V3 = Vinduced if numCus < 1 else Vcuspidal

char1, char2, char3 = 0, 0, 0



################# CHANGE HERE!

#decomposeChars()
#validCombinations()

'''
header = createHeader()
print(len(header))
print("")


reps = [
        (0, 1, 0)
        ]

data = [header]
char1, char2, char3 = 0, 0, 0

for rep in reps:
    char1 = goodCharsB[rep[0]]
    char2 = goodCharsB[rep[1]]
    char3 = nondecomposableChars[rep[2]]

    v1 = V1.basis()[4]
    v2 = V2.basis()[3]
    v3 = V3.basis()[2]

    calc = trilinearForm(v1, v2, v3)
    print(identify(calc, ['sqrt(5)']))
    
    #s = triformBasisVecs()
    #s.insert(0, "rho" + str(rep[0]) + " rho" + str(rep[1]) + " rho" + str(rep[2]))
    #data.append(s)
    #print("Completed rep " + str(rep))


'''
'''
header = createHeader()
combos = validCombinations()
data = [header]
countreps = 0
for combo in combos:
    countreps = countreps + 1
    print(combo)
    char1 = goodCharsB[combo[0]] if numCus < 3 else nondecomposableChars[combo[0]]
    char2 = goodCharsB[combo[1]] if numCus < 2 else nondecomposableChars[combo[1]]
    char3 = goodCharsB[combo[2]] if numCus < 1 else nondecomposableChars[combo[2]]

    s = triformBasisVecs()
    if numCus == 3:
        s.insert(0, "rho" + str(combo[0]) + " rho" + str(combo[1]) + " rho" + str(combo[2]))
    elif numCus == 2:
        s.insert(0, "pi" + str(combo[0]) + " rho" + str(combo[1]) + " rho" + str(combo[2]))
    elif numCus == 1:
        s.insert(0, "pi" + str(combo[0]) + " pi" + str(combo[1]) + " rho" + str(combo[2]))
    data.append(s)
    print("Completed rep " + str(combo))
    print("This was " + str(countreps) + " out of " + str(len(combos)) + "\n")


with open(str(numCus) + 'cusp_q' + str(q) + '.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
'''


#################

# https://mpmath.org/doc/current/identification.html#identify