import csv
from mpmath import *
import pickle
from pathlib import Path
import os
import numpy as np
import time

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

V1 = Vinduced if numCus < 3 else Vcuspidal
V2 = Vinduced if numCus < 2 else Vcuspidal
V3 = Vinduced if numCus < 1 else Vcuspidal

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



# shouldnt need to use array here, not significant enough
goodCharsB = []
nondecomposableChars = []

with open('precomputedChars/' + str(numCus) + 'cusp_q' + str(q) + '_goodCharsB', 'rb') as fp:
    goodCharsB = pickle.load(fp)

with open('precomputedChars/' + str(numCus) + 'cusp_q' + str(q) + '_nondecomposableChars','rb') as fp:
    nondecomposableChars = pickle.load(fp)



# Fixing character of U (this choice doesn't matter)
ct = U.character_table()
charU = U.character(ct[1])
#print("This is the character of K we are using: ")
#print(charU.values())
#print("")

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

# Char of L - takes in elem in L and an elem in nondecomposableChars
def nu(l, nondecompChar):
    v = Lx(MSforL([l]))
    return nondecompChar(v)


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
                    #print("Triple of two good characters and one cuspidal rep with trivial central character:")
                    #print(chars, end=", ")
                    #print(i)
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
                    #print("Triple of one good character and two cuspidal reps with trivial central character:")
                    #print(i, end=", ")
                    #print(chars)
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
                #print("Triple of three cuspidal reps with trivial central character:")
                #print(chars)
                count = count + 1
                combos.append((chars[0], chars[1], chars[2]))
        print("There are " + str(count) + " combos")
    print("\n")
    return combos
#


# For the induced irreps:

# Takes in element g from G and returns the corresponding representative in B\G
def toRepresentativeInduced(g):
    if g.inverse().matrix()[0][0] == 0:
        return G(MS([[0,1],[1,0]])).inverse()
    else:
        x = (K.one() / g.inverse().matrix()[0][0]) * g.inverse().matrix()[1][0]
        return G(MS([[1,0],[x,1]])).inverse()
# Could be optimized, lots of inverses?

#sum1, sum2, sum3, sum4 = 0
#count1, count2, count3, count4 = 0
# Gives the G action result of group element g from G onto vector v from V, with rep induced by chi
def gActionInduced(g, vec, chi):
    #global sum1, count1, sum2, count2, sum3, count3, sum4, count4

    newVec = Vinduced([0] * (q+1))
    for i in range(0, q+1):
        if vec[i] == 0:
            continue
        

        #start = time.time()

        newRep = toRepresentativeInduced(cosetRepsB[i] * g.inverse())

        #end = time.time()
        #sum1 = sum1 + end - start
        #count1 = count1 + 1

        #start = time.time()
        b = cosetRepsB[i] * g.inverse() * newRep.inverse()
        #end = time.time()
        #sum2 = sum2 + end - start
        #count2 = count2 + 1
        
        newIndex = cosetRepsB.index(newRep)

        #start = time.time()
        newVec[newIndex] = newVec[newIndex] + chi(b.inverse()) * vec[i]
       # end = time.time()
        #sum3 = sum3 + end - start
        #count3 = count3 + 1

       # start = time.time()
        g.inverse()
       # end = time.time()
        #sum4 = sum4 + end - start
        #count4 = count4 + 1

    return newVec
# CAN POSSIBLY OPTIMIZE THIS TOO???


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
            y = basisRepsCuspidal[i]
            newVec[i] = newVec[i] + coeff(y, oldRep, g, nondecompChar) * vec[i] # TEMPORARY DELETE AFTERWARDS!!!!! MAYBE MORE EFFICIENT THIS WAY??? Only for purposes of dot product with itself
            #for j in range(0, q-1):
            #    y = basisRepsCuspidal[j]
            #    newVec[j] = newVec[j] + coeff(y, oldRep, g, nondecompChar) * vec[i]     
        return newVec
# LOOK INTO OPTIMIZING THIS!!! - this is probably the bottleneck


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

#sumOfTimesCuspidal = 0
#countCuspidal = 0
# Matrix coeff of cuspidal
def matrixCoeffCuspidal(g, vec1, vec2, nondecompChar):
    #global sumOfTimesCuspidal
    #global countCuspidal
   # start = time.time()
    v = gActionCuspidal(g, vec1, nondecompChar)
    #end = time.time()
    #print("Matrix coeff cusp time: ")
    #print(end - start)
    #sumOfTimesCuspidal = sumOfTimesCuspidal + end - start
    #countCuspidal = countCuspidal + 1

    return v.dot_product(conjugate(vec2))
# 


#sumOfTimesInduced = 0
#countInduced = 0
# Matrix coeff of Induced (no matter what irrep in particular)
def matrixCoeffInduced(g, vec1, vec2, chi):
    #global sumOfTimesInduced
    #global countInduced
    #start = time.time()
    v = gActionInduced(g, vec1, chi)
    #end = time.time()
    #sumOfTimesInduced = sumOfTimesInduced + end - start
    #countInduced = countInduced + 1

    return v.dot_product(conjugate(vec2))
# 


# Gives the triple product - have to specify the representations inside the function
def tripleProduct(g, v1, v2, v3):
    if numCus < 3:
        one = matrixCoeffInduced(g, v1, v1, char1)
    else:
        one = matrixCoeffCuspidal(g, v1, v1, char1)
    if one == 0:
        return 0

    if numCus < 2:
        two = matrixCoeffInduced(g, v2, v2, char2)
    else:
        two = matrixCoeffCuspidal(g, v2, v2, char2)
    if two == 0:
        return 0
    

    if numCus < 1:
        three = matrixCoeffInduced(g, v3, v3, char3)
    else:
        three = matrixCoeffCuspidal(g, v3, v3, char3)

    return one * two * three

# Gives trilinear form by avergaing matrix coefficients over the whole group
def trilinearForm(v1, v2, v3):
    l = np.empty(G.order(), dtype=object)
    i = 0
    for g in G:
        l[i] = tripleProduct(g, v1, v2, v3)
        i = i + 1
    s = np.sum(l)

    # Use above for memory improvements!!
    #l = [tripleProduct(g, v1, v2, v3) for g in G]
    #s = sum(l)

    sol = (s / G.order())
    if sol.imag() < 0.000001:
        sol = sol.real()
    try:
        return QQ(sol)
    except:
        return sol
# We don't know how to find Whittaker fns of cuspidal reps, only matrix coeffs

basisCombos = cartesian_product((range(V1.dimension()), range(V2.dimension()), range(V3.dimension())))

# Iterates over all basis vectors and computes the trilinear forms of them
def triformBasisVecs(basisCombo):
    sol = []
    (i, j, k) = basisCombos[basisCombo]

    v1 = V1.basis()[i]
    v2 = V2.basis()[j]
    v3 = V3.basis()[k]
    if v1 == v2 or v1 == v3 or v2 == v3:
        print("Not all different basis vecs")
    calc1 = trilinearForm(v1, v2, v3)
    if calc1 < 0.0000001:
        calc1 = 0
    s = calc1

    print(s)
    sol.append(s)
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





char1, char2, char3 = 0, 0, 0





combos = []

if Path("data/" + str(numCus) + 'cusp_q' + str(q) + "/representation_combinations.csv").exists():
    with open("data/" + str(numCus) + 'cusp_q' + str(q) + "/representation_combinations.csv", 'r') as csvfile:
        for row in csv.reader(csvfile):
            combos.append(row)
else:
    combos = validCombinations()
    with open("data/" + str(numCus) + 'cusp_q' + str(q) + "/representation_combinations.csv", 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(combos)




print("There are " + str(len(combos)) + " representations, input which you want to find: ")
whichcombo = int(input(""))
combo = combos[whichcombo]

print("There are " + str(V1.dimension() * V2.dimension() * V3.dimension()) + " triples of basis vecs, input which you want to find: ")
whichBasis = int(input(""))
print("")

# Need to include safeguards if inputs outside of range???

folder_path = "data/" + str(numCus) + 'cusp_q' + str(q) + "/" + str(whichcombo).zfill(3)
Path(folder_path).mkdir(parents=True, exist_ok=True)


print(combo)
char1 = goodCharsB[int(combo[0])] if numCus < 3 else nondecomposableChars[int(combo[0])]
char2 = goodCharsB[int(combo[1])] if numCus < 2 else nondecomposableChars[int(combo[1])]
char3 = goodCharsB[int(combo[2])] if numCus < 1 else nondecomposableChars[int(combo[2])]


if Path(folder_path + "/" + str(whichBasis).zfill(4) + '.csv').exists():
    print("Already done!")
else:

    s = triformBasisVecs(whichBasis)
    data = []
    if numCus == 3:
        s.insert(0, "rho" + str(combo[0]) + " rho" + str(combo[1]) + " rho" + str(combo[2]))
    elif numCus == 2:
        s.insert(0, "pi" + str(combo[0]) + " rho" + str(combo[1]) + " rho" + str(combo[2]))
    elif numCus == 1:
        s.insert(0, "pi" + str(combo[0]) + " pi" + str(combo[1]) + " rho" + str(combo[2]))
    data.append(s)
    with open(folder_path + "/" + str(whichBasis).zfill(4) + '.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
        print("Completed basis " + str(basisCombos[whichBasis]) + " of rep " + str(combo))


# Optimization testing:
'''
print("Average time for matrixCoeffInduced:")
if (countInduced == 0):
    print(0)
else:
    print(sumOfTimesInduced / countInduced)
print(countInduced)
print("Average time for matrixCoeffCuspidal:")
print(sumOfTimesCuspidal / countCuspidal)
print(countCuspidal)




print("Average time for sum1:")
print(sum1 / count1)
print(count1)
print("Average time for sum2:")
print(sum2 / count2)
print(count2)
print("Average time for sum3:")
print(sum3 / count3)
print(count3)
print("Average time for sum4:")
print(sum4 / count4)
print(count4)
'''
