

# This file was *autogenerated* from the file cuspidal_no_whittaker.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_15 = Integer(15)# Preliminaries

# Field and general linear group
q = _sage_const_3 
K = GF(q)
#Ltest = GF(q^2).over(K)
L = GF(q**_sage_const_2 )


G = GL(_sage_const_2 , K)
MS = MatrixSpace(K, _sage_const_2 )


# U subgroup
gens = [MS([[_sage_const_1 ,x],[_sage_const_0 ,_sage_const_1 ]]) for x in K]
U = G.subgroup(gens)


# Creating vector space
dim = q-_sage_const_1 
V = VectorSpace(CC, dim)
#This vector space of functions K^x to C, basis is fns that take value of 1 on one element and zero elsewehere



# Fix basis:
basisRepsCuspidal = K.list()
basisRepsCuspidal.remove(_sage_const_0 )
# Use .index() to find index





# Getting psi 
ct = U.character_table()
charU = U.character(ct[_sage_const_1 ])
print("This is the character of K we are using: ")
print(charU.values())
print("")

def psi(x):
    return charU(G(MS([[_sage_const_1 , x], [_sage_const_0 , _sage_const_1 ]])))





# Getting nondecomp characters of L^x
Lx = GL(_sage_const_1 , L)
MSforL = MatrixSpace(L, _sage_const_1 )

#H = Hom(K, L)
#inclusionMap = H[0]

ct2 = Lx.character_table()
#print(ct2)

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


# NEED TO REMOVE CONJUAGETS!
# Currently double counting the cuspidal reps



print("This is the non-decomposable character of L^x we are using: ")
print(nondecomposableChars[_sage_const_0 ].values())
print("")


# For now, fix nondecomp character:
def nu(l):
    m = nondecomposableChars[_sage_const_0 ]
    v = Lx(MSforL([l]))
    return m(v)










# Group action as described pg 40 in P-S
def gAction(g, vec):
    if g.matrix()[_sage_const_1 , _sage_const_0 ] == _sage_const_0 :
        # Easier implementation, g is in B
        newVec = V([_sage_const_0 ] * dim)

        a = g.matrix()[_sage_const_0 , _sage_const_0 ]
        b = g.matrix()[_sage_const_0 , _sage_const_1 ]
        d = g.matrix()[_sage_const_1 , _sage_const_1 ]

        for i in range(_sage_const_0 , dim):
            if vec[i] == _sage_const_0 :
                continue
    
            oldRep = basisRepsCuspidal[i]

            newRep = d * (_sage_const_1  / a) * oldRep
            newIndex = basisRepsCuspidal.index(newRep)
            coefficient = nu(d) * psi(b * (_sage_const_1  / a) * oldRep)

            newVec[newIndex] = newVec[newIndex] + coefficient * vec[i]
        
        return newVec

    else:
        # Harder longer formula
        newVec = V([_sage_const_0 ] * dim)

        for i in range(_sage_const_0 , dim):
            if vec[i] == _sage_const_0 :
                continue
            
            oldRep = basisRepsCuspidal[i]
            for j in range(_sage_const_0 , dim):
                y = basisRepsCuspidal[j]
                newVec[j] = coeff(y, oldRep, g) * vec[i]
        
        return newVec
# Fast


# Helper function for gAction - calculates the coefficients when g is not in B
def coeff(y, x, g):
    a = g.matrix()[_sage_const_0 , _sage_const_0 ]
    b = g.matrix()[_sage_const_0 , _sage_const_1 ]
    c = g.matrix()[_sage_const_1 , _sage_const_0 ]
    d = g.matrix()[_sage_const_1 , _sage_const_1 ]
    temp = _sage_const_0 

    comparison = y * (_sage_const_1  / x) * (a*d - b*c)
    for u in L:
        if u != _sage_const_0  and u * conjugateOfL(u) == comparison:
            temp = temp + nu(u) * psi(- (x / c) * (u + conjugateOfL(u)))
    
    return temp * psi((a * y + d * x) / c) / q
# Fast as well



g = G.random_element()
#g = G(MS([[1, 2], [1, 1]]))
print(g)
vec = V([_sage_const_0 ]*dim)
vec[_sage_const_0 ] = _sage_const_1 

print(gAction(g, vec))











# With testing so far, it seems that 

def innerProduct(vec1, vec2):
    #sol = 0
    #for elem in G: # Roughly q^4 elems
        #temp = (gAction(elem, vec1)).dot_product(conjugate(gAction(elem, vec2)))
        #print(temp)
        #sol = sol + temp

    sol = sum([(gAction(elem, vec1)).dot_product(conjugate(gAction(elem, vec2))) for elem in G])
    sol = sol / G.order()
    return round(sol.real(), _sage_const_5 ) + round(sol.imag(), _sage_const_5 ) * I

def bigMatrix(g):
    basisMatrix = matrix(CC, dim, dim, _sage_const_0 )
    for i in range(_sage_const_0 , dim):
        for j in range(i, dim):
            vec1 = V([_sage_const_0 ]*(dim))
            vec2 = V([_sage_const_0 ]*(dim))
            vec1[i] = _sage_const_1 
            vec2[j] = _sage_const_1 

            x = innerProduct(vec1, vec2)
            basisMatrix[i, j] = x
            basisMatrix[j, i] = conjugate(x)
        print("Row " + str(i) + " of basis matrix")
    print(basisMatrix)

    M = matrix(CC, dim, dim, _sage_const_0 )
    for i in range(_sage_const_0 , dim):
        for j in range(_sage_const_0 , dim):
            vec1 = V([_sage_const_0 ]*(dim))
            vec2 = V([_sage_const_0 ]*(dim))
            vec1[i] = _sage_const_1 
            vec2[j] = _sage_const_1 
            
            realvec1 = gAction(g, vec1)

            M[i, j] = realvec1.dot_product(basisMatrix.column(j))
            #M[i, j] = realvec1.dot_product(vec2)
            #M[i, j] = innerProduct(gAction(g, vec1), vec2)
            print("Filled in entry in row " + str(i) + " and column " + str(j))
    return M

# Prints the matrix nicely, rounding the actual digits and displays zero as zero
def printMatrix(M):
    decimalPlaces = _sage_const_4 
    print("Printing matrix:")
    for row in M:
        print("[", end="")
        for elem in row:
            temp = _sage_const_0 
            if elem != _sage_const_0 :
                temp = round(elem.real(), decimalPlaces) + round(elem.imag(), decimalPlaces) * I
            print("{:<{mx}}".format(temp, mx=_sage_const_15 ), end="\t")
        print("]")
# Should automate length of mx


M = bigMatrix(g)
printMatrix(M)

