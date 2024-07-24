# Combine are the representations here and calculate the trilinear form

# define inner product on cuspidal? - seems like normal dot product works


# Field and general linear group
q = 5
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





# Finding combos with trivial central character:

# 2 q+1 dim, 1 q dim
'''
C = Combinations(list(range(len(goodCharsB))) * 2, 2)
print("")

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
            count = count+1
print("There are " + str(count)+" valid combinations")
'''


# 1 q+1 dim, 2 q dim
'''
C = Combinations(list(range(len(badCharsB))) * 2, 2)
print("")

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
            count = count+1
print("There are " + str(count)+" valid combinations")
'''


# 3 q dim - not fixed yet

C = Combinations(list(range(len(badCharsB))) * 3, 3)
print("")

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
        memorySet.add(t[1])
    #print(memorySet)

    for g in G:
        if len(memorySet) == 1 and list(memorySet)[0].dimension() == 0:
            break

        # ONLY FOR DEALING WITH BAD CHARACTERS!!! MAY CAUSE ERRORS WHEN TESTING THIS ON THE GOOD CHARACTERS!!!
        if len(memorySet) == 2:
            l = []
            for x in memorySet:
                l.append(x.dimension())
            if l == [0, 1] or l == [1, 0]:
                break
        # Although reduces accuracy, the speed is improved 100 fold

        spaces = eigenSpaces(g, chi)
        gSet = set()
        for t in spaces:
            gSet.add(t[1])

        tempSet = set()
        for ogSpace in memorySet:
            for newSpace in gSet:
                t = ogSpace.intersection(newSpace)
                tempSet.add(t)
        #print(tempSet)
        memorySet = tempSet

    if len(memorySet) == 1:
        print("This was a good character! No G-invariant subspace!")
    else:
        print("This is the 1d G-invariant subspace:")
        for item in memorySet:
            if item.dimension()==1:
                print(item)
                return item
# Runs pretty slowly when bad character - any way to speed this up?
# Could start checking if only 1d subspace left, then just simply check if this remains to be eigenvector for remainding elems

# ABOVE FUNCTION IS SO UNNECESSARY, WE HAVE EXPLICIT FORMULA FOR THIS VECTOR!!!







# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    v = gActionInduced(g, vec1, chi)
    return v.dot_product(conjugate(vec2))
# Do actual sum manually





# CHANGE CHARACTERS USED HERE!!!

induced1 = badCharsB[1]
induced2 = badCharsB[2]
induced3 = badCharsB[3]


OneDSpace1 = findGsubspace(induced1)
Vqdim1 = OneDSpace1.complement()

OneDSpace2 = findGsubspace(induced2)
Vqdim2 = OneDSpace2.complement()

OneDSpace3 = findGsubspace(induced3)
Vqdim3 = OneDSpace3.complement()



#vec1 = Vinduced.basis()[0]
#vec2 = Vinduced.basis()[1]
#vec3 = Vinduced.basis()[2]



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
    return s / G.order()

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




# Iterates over all basis vectors and computes the trilinear forms of them

for i, j, k in cartesian_product((range(q), range(q), range(q))):
    v1 = Vqdim1.basis()[i]
    v2 = Vqdim2.basis()[j]
    v3 = Vqdim3.basis()[k]
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