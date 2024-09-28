################# CHANGE HERE!

q = 5


# Set how many q and q+1 dimensional we want - must add up to 3
numQp1 = 3
numQ = 0

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
Vinduced = VectorSpace(QQbar, q+1)
H = Hom(Vinduced, Vinduced)








# Representatives
# Use .index() to find index

# Reps of B
cosetRepsB = []
for x in K:
    rep = G(MS([[1,0],[-x,1]]))
    cosetRepsB.append(rep)
rep = G(MS([[0,1],[1,0]]))
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
def decomposeChars():
    test1 = G(MS([[a, 0], [0, 1]]))
    test2 = G(MS([[1, 0], [0, a]]))
    print("Our generator of K^x is " + str(a))
    print("")
    if numQp1 != 0:
        for i in range(len(goodCharsB)):
            print("This is good character " + str(i) + ", and the following are the two values of chi_1 and chi_2 on the generator:")
            print(goodCharsB[i](test1))
            print(goodCharsB[i](test2))
            print("")
        print("")

    if numQ != 0:
        for i in range(len(badCharsB)):
            print("This is bad character " + str(i) + ", and the following is the values of chi_1 on the generator:")
            print(badCharsB[i](test1))
            print("")
        print("")


combos = []

# Finding combos with trivial central character:
def validCombinations():
    # 3 q+1
    if numQ == 0:
        C = Combinations(list(range(len(goodCharsB))) * 3, 3)
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
                combos.append(chars)
        print("")
    # 2 q+1 dim, 1 q dim
    if numQ == 1:
        C = Combinations(list(range(len(goodCharsB))) * 2, 2)
        count = 0
        for chars in C:
            for i in range(0, len(badCharsB)):
                chi1 = goodCharsB[chars[0]]
                chi2 = goodCharsB[chars[1]]
                chi3 = badCharsB[i]

                centralCharTrivial = True
                for x in K:
                    if x != 0:
                        d = G(MS([[x, 0], [0, x]]))
                        prod = chi1(d) * chi2(d) * chi3(d)
                        if prod != 1:
                            centralCharTrivial = False
                            break
                if centralCharTrivial:
                    print("Triple of 2 good characters and 1 bad one with trivial central character:")
                    print(chars, end=", ")
                    print(i)
                    combos.append((chars[0], chars[1], i))
                    count = count+1
        print("There are " + str(count) + " valid combinations")


    # 1 q+1 dim, 2 q dim
    if numQ == 2:
        C = Combinations(list(range(len(badCharsB))) * 2, 2)
        count = 0
        for chars in C:
            for i in range(0, len(goodCharsB)):
                chi1 = badCharsB[chars[0]]
                chi2 = badCharsB[chars[1]]
                chi3 = goodCharsB[i]

                centralCharTrivial = True
                for x in K:
                    if x != 0:
                        d = G(MS([[x, 0], [0, x]]))
                        prod = chi1(d) * chi2(d) * chi3(d)
                        if prod != 1:
                            centralCharTrivial = False
                            break
                if centralCharTrivial:
                    print("Triple of 1 good character and 2 bad ones with trivial central character:")
                    print(i, end=", ")
                    print(chars)
                    combos.append((i, chars[0], chars[1]))
                    count = count+1
        print("There are " + str(count)+" valid combinations")



    # 3 q dim
    if numQ == 3:
        C = Combinations(list(range(len(badCharsB))) * 3, 3)
        count = 0
        for chars in C:
            chi1 = badCharsB[chars[0]]
            chi2 = badCharsB[chars[1]]
            chi3 = badCharsB[chars[2]]

            centralCharTrivial = True
            for x in K:
                if x != 0:
                    d = G(MS([[x, 0], [0, x]]))
                    prod = chi1(d) * chi2(d) * chi3(d)
                    if prod != 1:
                        centralCharTrivial = False
                        break
            if centralCharTrivial:
                print("Triple of 3 bad characters with trivial central character:")
                print(chars)
                combos.append(chars)
                count = count+1
        print("There are " + str(count)+" valid combinations")









            

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


# Finds the 1d subsapce given a bad character
def findGsubspace(chi):
    sol = Vinduced([1]*(q+1))
    sol[q] = chi(G(MS([[-1, 0], [0, 1]])))
    return Vinduced.subspace([sol])
# Constant time




print(rep)


# Evaluates function vec at value g
def evalInduced(g, vec, char):
    rep = toRepresentativeInduced(g)
    index = cosetRepsB.index(rep)
    if vec[index] == 0:
        return 0
    b = g * rep.inverse()
    return vec[index] * char(b)


# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    v = gActionInduced(g, vec1, chi)
    return v.dot_product(conjugate(vec2))


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
        val3 = whittaker2(g, v3, induced3)
        if val3 == 0:
            continue
        #print(g)
        #print(val1 * val2 * val3)
        s = s + val1 * val2 * val3

    s = s / G.order()
    #print(s)
    return QQ(norm(s))

    '''
    if almostSol == 0:
        return 0

    normalizingFactor1 = 0
    for uRep in cosetRepsU:
        normalizingFactor1 = normalizingFactor1 + norm(whittaker(uRep, v2, induced2))
    normalizingFactor1 = normalizingFactor1 / (q-1)^2

    normalizingFactor2 = 0
    for uRep in cosetRepsU:
        normalizingFactor2 = normalizingFactor2 + norm(whittaker2(uRep, v3, induced3))
    normalizingFactor2 = normalizingFactor2 / (q-1)^2

    print(normalizingFactor1)
    print(normalizingFactor2)

    #return QQ(almostSol / normalizingFactor1 / normalizingFactor2)
    return almostSol
    '''
print("")



# Helper function
def triformBasisVecsH(i, j, k):
    v1 = V1.basis()[i]
    v2 = V2.basis()[j]
    v3 = V3.basis()[k]
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

# Iterates over all basis vectors and computes the trilinear forms of them
def triformBasisVecs():
    for i, j, k in cartesian_product((range(V1.dimension()), range(V2.dimension()), range(V3.dimension()))):
        triformBasisVecsH(i, j, k)
# This is a function now so below can be easily commented in/out to toggle


# Hits the varieties that matter
def triformBasisVecsShort():
    if numQ == 1:
        triformBasisVecsH(0, 0, 0)
        triformBasisVecsH(1, 0, 0)
        triformBasisVecsH(0, 1, 0)
        triformBasisVecsH(0, 0, 1)
        triformBasisVecsH(0, q, 0)
        triformBasisVecsH(q, 0, 0)
    if numQ == 0:
        triformBasisVecsH(0, 1, 2)
        triformBasisVecsH(0, 0, 0)
        triformBasisVecsH(1, 0, 0)
        triformBasisVecsH(0, 1, 0)
        triformBasisVecsH(0, 0, 1)




################# CHANGE HERE!

#decomposeChars()
#validCombinations()


induced1 = goodCharsB[0]
induced2 = goodCharsB[2]
induced3 = goodCharsB[3]


V1 = Vinduced if induced1 in goodCharsB else findGsubspace(induced1).complement()
V2 = Vinduced if induced2 in goodCharsB else findGsubspace(induced2).complement()
V3 = Vinduced if induced3 in goodCharsB else findGsubspace(induced3).complement()


vec1 = Vinduced.basis()[q]
vec2 = Vinduced.basis()[0]
vec3 = Vinduced.basis()[1]

#print(RStrilinearForm(vec1, vec2, vec3))



#triformBasisVecs()
triformBasisVecsShort()

#################