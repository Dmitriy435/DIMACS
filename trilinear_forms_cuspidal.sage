################# CHANGE HERE!

q = 5


# Set how many induced and cuspidal we want - must add up to 3

numInd = 1 #For now only working with good chars!
numCus = 2 #This is at least 1 for this file

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
        print("There are " + str(count) + " combos")
    print("\n")









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


# Helper function for gAction - calculates the coefficients when g is not in B
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








# Evaluates function vec at value g
def evalInduced(g, vec, char):
    rep = toRepresentativeInduced(g)
    index = cosetRepsB.index(rep)
    if vec[index] == 0:
        return 0
    b = g * rep.inverse()
    return vec[index] * char(b)


# Matrix coeff of cuspidal
def matrixCoeffCuspidal(g, vec1, vec2, nondecompChar):
    v = gActionCuspidal(g, vec1, nondecompChar)
    return v.dot_product(conjugate(vec2))


# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    v = gActionInduced(g, vec1, chi)
    return v.dot_product(conjugate(vec2))


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
    for i, j, k in cartesian_product((range(V1.dimension()), range(V2.dimension()), range(V3.dimension()))):
        v1 = V1.basis()[i]
        v2 = V2.basis()[j]
        v3 = V3.basis()[k]
        if v1 == v2 or v1 == v3 or v2 == v3:
            print("Not all different basis vecs")
        calc1 = trilinearForm(v1, v2, v3)
        if calc1 < 0.0000001:
            calc1 = 0
        print(calc1)
        if calc1 != 0:
            print(v1)
            print(v2)
            print(v3)
        print("")
# How to do Whittaker?






################# CHANGE HERE!

decomposeChars()
validCombinations()


char1 = goodCharsB[0]
char2 = nondecomposableChars[7]
char3 = nondecomposableChars[8]


V1 = Vinduced if char1 in goodCharsB else Vcuspidal
V2 = Vinduced if char2 in goodCharsB else Vcuspidal
V3 = Vinduced if char3 in goodCharsB else Vcuspidal



triformBasisVecs()

#################