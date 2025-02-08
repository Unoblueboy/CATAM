# -*- coding: utf-8 -*-
"""
Created on Sun Dec  2 21:44:09 2018

@author: Natha
"""
import numpy as np

# QUESTION 1

def inverses(p):
    inv = []
    for a in range(1,p):
        a_inv = 1
        for a_inv in range(1,p):
            if a*a_inv % p == 1:
                inv.append(a_inv)
                break
        else:
            inv.append(0)
    return inv
    

# QUESTION 3

def transpose(M, i, j ,mod):
    # swap rows i and j
    M = M.copy()
    y = M[i].copy()
    M[i] = M[j].copy()
    M[j] = y
    return M

def multiply(M, i, a, mod):
    # multiply row i by a
    M = M.copy()
    M[i] = a*M[i] % mod
    return M

def subtract(M, i, j, a, mod):
    # subtract a * row j from row i
    M = M.copy()
    M[i] = (M[i] - a*M[j]) % mod
    return M

def gauss_elim(M, mod):
    M = M.copy()
    M = M % mod
    rows = M.shape[0]
    cols = M.shape[1]
    inv = inverses(mod)
    cur_col = 0
    cur_row = 0
    while cur_col < cols and cur_row < rows:
        if M[cur_row][cur_col] == 0:
            for row in range(cur_row, rows):
                if M[row][cur_col]!=0:
                    M=transpose(M, row, cur_col, mod)
                    break
            else:
                cur_col += 1
                continue
        lead_coef = M[cur_row][cur_col]
        M = multiply(M, cur_row, inv[lead_coef-1], mod)
        for row in range(0, rows):
            if row == cur_row:
                continue
            row_lead_coef = M[row][cur_col]
            M = subtract(M, row, cur_row, row_lead_coef, mod)
        cur_row+=1
        cur_col+=1
    return M
 
A1 = np.array([[0,1,7,2,10],[8,0,2,5,1],[2,1,2,5,5],[7,4,5,3,0]])
A2 = np.array([[6,16,11,14,1,4],[7,9,1,1,21,0],[8,2,9,12,17,7],[2,19,2,19,7,12]])

print("GAUSS ELIMINATION")
print("A1, 11")
print(gauss_elim(A1,11))

print("A1, 19")
print(gauss_elim(A1,19))

print("A2, 23")
print(gauss_elim(A2,23))
print()

#QUESTION 4
       
def rank(M, mod):
    M = gauss_elim(M,mod)
    r=0
    rows = M.shape[0]
    cols = M.shape[1]
    cur_col = 0
    cur_row = 0
    while cur_col < cols and cur_row < rows:
        if M[cur_row][cur_col] > 0:
            r+=1
            cur_row+=1
        else:
            cur_col+=1
    return r

def undetermined_x(M, mod):
    M = gauss_elim(M,mod)
    rows = M.shape[0]
    cols = M.shape[1]
    cur_col = 0
    cur_row = 0
    undet_x = []
    while cur_col < cols and cur_row < rows:
        if M[cur_row][cur_col] > 0:
            cur_row+=1
            cur_col+=1
        else:
            undet_x.append(cur_col)
            cur_col+=1
    if cur_col != cols:
        undet_x = undet_x + [x for x in range(cur_col, cols)]
    return undet_x

def kernel_basis(M, mod):
    M = gauss_elim(M ,mod)
    r = rank(M, mod)
    cols = M.shape[1]
    undet_x = undetermined_x(M, mod)
    basis = []
    for i in undet_x:
        x = np.zeros((cols,1), dtype=int)
        x[i][0]=1
        j=cols-1
        cur_row = r-1
        determined_x = np.matmul(M, x) 
        while j>=0:
            if j not in undet_x:
                x[j] = -determined_x[cur_row] % mod
                cur_row -= 1
            j-=1
        basis.append(x)
    if len(basis) == 0:
        basis.append(np.zeros((cols,1), dtype=int))
    return basis
        
B1 = np.array([[4,6,5,2,3],[5,0,3,0,1],[1,5,7,1,0],[5,5,0,3,1],[2,1,2,4,0]])
B2 = np.array([[3,7,19,3,9,6],[10,2,20,15,3,0],[14,1,3,14,11,3],[26,1,21,6,3,5],[0,1,3,19,0,3]])

print("KERNEL BASES")
print("B1, 13")
print(kernel_basis(B1, 13))

print("B1, 17")
print(kernel_basis(B1, 17))

print("B1, 13")
print(kernel_basis(B2, 23))
print()

# QUESTION 5
# dim(U) + dim(Uo) = n

def get_rid_of_empty(U):
    rows = U.shape[0]
    cols = U.shape[1]
    non_empty = []
    for i in range(rows):
        for j in range(cols):
            if U[i][j]!=0:
                non_empty.append(i)
                break
    U=U[non_empty]
    return U

def vector_row_space_sum(U, W, mod):
    # U and W given in the form of a matrix where U has row space U
    # and W has row space W
    U = gauss_elim(U,mod)
    W = gauss_elim(W,mod)
    full_space = np.concatenate((U, W))
    full_space = gauss_elim(full_space,mod)
    full_space=get_rid_of_empty(full_space)
    return full_space

def vector_col_space_sum(U, W, mod):
    row_U = np.transpose(U)
    row_W = np.transpose(W)
    row_full_space = vector_row_space_sum(row_U, row_W, mod)
    full_space = np.transpose(row_full_space)
    return full_space

def row_annihilator(U, mod):
    # U given in the form of a matrix with row space U
    basis = kernel_basis(U, mod)
    matrix_basis = np.concatenate(basis, axis=1)
    return matrix_basis

def col_annihilator(U, mod):
    return np.transpose(row_annihilator(np.transpose(U),mod))
    
# QUESTION 6
U_annihilate = row_annihilator(A1, 19)
U_annihilate_annihilate = col_annihilator(U_annihilate, 19)
print("ANNIHILATION")
print("Uo")
print(U_annihilate)

print("Uoo")
print(U_annihilate_annihilate)
print()

#QUESTION 7

def bases(U, W, mod):
    # U
    U = gauss_elim(U, mod)
    U = get_rid_of_empty(U)

    
    # W
    W = gauss_elim(W, mod)
    W = get_rid_of_empty(W)

    
    # U + W
    UpW = vector_row_space_sum(U, W, mod)

    
    # U & W
    UaW = col_annihilator(vector_col_space_sum(row_annihilator(U, mod),row_annihilator(W, mod), mod), mod)

    
    return U, W, UpW, UaW

A3 = np.array([[1,0,0,0,3,0,0],[0,5,0,1,6,3,0],[0,0,5,0,2,0,0],[2,4,0,0,0,5,1],[4,3,0,0,6,2,6]])

print("BASES OF SUMS AND INTERSECTIONS")
print("A1, B1")
U1, W1, U1pW1, U1aW1 = bases(A1, B1, 11)
print(U1)
print(W1)
print(U1pW1)
print(U1aW1)

B3 = np.transpose(row_annihilator(A3, 19))
print("A3, B3")
U2, W2, U2pW2, U2aW2 = bases(A3, B3, 19)
print(U2)
print(W2)
print(U2pW2)
print(U2aW2)

B4 = np.transpose(row_annihilator(A3, 23))
print("A3, B4")
U3, W3, U3pW3, U3aW3 = bases(A3, B4, 23)
print(U3)
print(W3)
print(U3pW3)
print(U3aW3)


# In the reals we would find U + W = R^n, so dim(U+W) = n
# but in the last case we see in fact that dim(U+W) = 6 = 7-1 so
# yh...