# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 22:45:42 2019

@author: Natha
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
import pandas as pd

def U1(X,T, N=150):
    answer = 1-X
    for n in range(1, N+1):
        term = -2*np.sin(n*np.pi*X)*np.e**-(n**2*np.pi**2*T)/(np.pi*n)
        answer += term
    return answer

def U2(X,T, N=150):
    answer = 1
    for n in range(0, N+1):
        term = -2*np.sin((n+0.5)*np.pi*X)*np.e**-((n+0.5)**2*np.pi**2*T)/(np.pi*(n+0.5))
        answer += term
    return answer

def U3(X,T, K = 1):
    return erfc(0.5*X/(K*T)**0.5) # assuming K = 1

X_Values = [0.125*n for n in range(9)]
T_Values = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0]

grid_1 = [U1(X, 0.25) for X in X_Values]
grid_2 = [U2(X, 0.25) for X in X_Values]
grid_3 = [U3(X, 0.25) for X in X_Values]

def draw_table(table, file_name, column_format):
    template = r'''\documentclass[preview]{{standalone}}
    \usepackage{{booktabs}}
    \begin{{document}}
    {}
    \end{{document}}
    '''
    dataframe = pd.DataFrame.from_dict(table)
    with open(file_name+".tex", 'w') as file:
        file.write(template.format(dataframe.to_latex(
                escape=False,
                column_format = column_format,
                index = False)))

Table_PT1 = {r"$x$":X_Values,
             r"\multicolumn{1}{|p{3cm}|}{\centering fixed-endpoint-temperature problem}":["{0:.6f}".format(y) for y in grid_1],
             r"\multicolumn{1}{|p{3cm}|}{\centering insulated-end problem}":["{0:.6f}".format(y) for y in grid_2],
             r"\multicolumn{1}{|p{3cm}|}{\centering semi-infinite bar problem}":["{0:.6f}".format(y) for y in grid_3]}

draw_table(Table_PT1, "table_PT1", "|c|l|l|l|")

grid_4 = [[U1(X, T) for X in X_Values] for T in T_Values]
grid_5 = [[U2(X, T) for X in X_Values] for T in T_Values]
grid_6 = [[U3(X, T) for X in X_Values] for T in T_Values]


fig = plt.figure(figsize=(13,5))


for i in range(len(grid_4)):
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(X_Values, grid_4[i], label="T="+str(T_Values[i]))
    ax1.legend()
    ax1.set_title("fixed-endpoint-temperature problem")
    ax1.set_xlabel("X")
    ax1.set_ylabel(r"$U_1\left(X,T\right)$")


for i in range(len(grid_5)):
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(X_Values, grid_5[i], label="T="+str(T_Values[i]))
    ax2.legend()
    ax2.set_title("insulated-end problem")
    ax2.set_xlabel("X")
    ax2.set_ylabel(r"$U_2\left(X,T\right)$")


for i in range(len(grid_6)):
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(X_Values, grid_6[i], label="T="+str(T_Values[i]))
    ax3.legend()
    ax3.set_title("semi-infinite bar problem")
    ax3.set_xlabel("X")
    ax3.set_ylabel(r"$f\left(\xi\right)$")

plt.tight_layout()
fig.savefig("PT1_1.pdf")



def U1_X(X,T, N=150):
    answer = -1
    for n in range(1, N+1):
        term = -2*np.cos(n*np.pi*X)*np.e**-(n**2*np.pi**2*T)
        answer += term
    return answer

def U2_X(X,T, N=150):
    answer = 0
    for n in range(1, N+1):
        term = -2*np.cos((n-0.5)*np.pi*X)*np.e**-((n-0.5)**2*np.pi**2*T)
        answer += term
    return answer

def U3_X(X,T, K = 1):
    return -np.e**(-X**2/(4*K*T))/(np.pi*K*T)**0.5 # assuming K = 1

grid_7 = [-U1_X(0, T) for T in T_Values]
grid_8 = [-U2_X(0, T) for T in T_Values]
grid_9 = [-U3_X(0, T) for T in T_Values]

fig = plt.figure(figsize=(13,5))

ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(T_Values, grid_7)
ax1.set_title("fixed-endpoint-temperature problem")
ax1.set_xlabel("t")
ax1.set_ylabel(r"$U_1\left(x,t\right)$")


ax2 = fig.add_subplot(1, 3, 2)
ax2.plot(T_Values, grid_8)
ax2.set_title("insulated-end problem")
ax2.set_xlabel("t")
ax2.set_ylabel(r"$U_2\left(x,t\right)$")


ax3 = fig.add_subplot(1, 3, 3)
ax3.plot(T_Values, grid_9)
ax3.set_title("semi-infinite bar problem")
ax3.set_xlabel("t")
ax3.set_ylabel(r"$f\left(\xi\right)$")


plt.tight_layout()
fig.savefig("PT1_2.pdf")


# Question 3

def U2_num(C, N, maxT = 2):
    dx = 1/N
    dt = C*dx*dx
    
    T_Values = [0]
    X_Values = [j*dx for j in range(0, N+1)]
    
    U = [[0.5] + [0 for _ in range(N)]]
    while T_Values[-1] < maxT:
        U_T = [1]
        T_Values.append(T_Values[-1] + dt)
        U_prev = U[-1]
        for i in range(1,N):
            U_i = U_prev[i] + C*(U_prev[i+1]-2*U_prev[i] + U_prev[i-1])
            U_T.append(U_i)
        U_N = U_prev[N] + C*(2*U_prev[N-1]-2*U_prev[N])
        U_T.append(U_N)
        U.append(U_T)
    return U, X_Values, T_Values

C_Values = [1/12, 1/6, 1/3, 1/2, 2/3]
N_Values = [8,16,32]

for i in range(len(C_Values)):
    fig = plt.figure(figsize=(13,5))
    axs = [fig.add_subplot(1, 3, l+1)  for l in range(len(N_Values))]
    C = C_Values[i]
    for j in range(len(N_Values)):
        ax = axs[j]
        N = N_Values[j]
        ax.set_title("N={}".format(N))
        grid, X_Vals, T_Vals = U2_num(C,N)
        for k in range(5):
            ax.plot(X_Vals, grid[k*len(grid)//5], label="t={:.4f}".format(T_Vals[k*len(grid)//5]))
        ax.plot(X_Vals, grid[-1], label="t={:.4f}".format(T_Vals[-1]))
        ax.legend()
        ax.set_xlabel(r"$X$")
        ax.set_ylabel(r"$U_2\left(X,T\right)$")
    plt.tight_layout()
    plt.show()
    fig.savefig("PT2_{}.pdf".format(i+1))

U2_ana_T_00625 = [U2(n/8, 0.0625) for n in range(9)]
U2_ana_T_0125 = [U2(n/8, 0.125) for n in range(9)]
U2_ana_T_025 = [U2(n/8, 0.25) for n in range(9)]
U2_ana_T_05 = [U2(n/8, 0.5) for n in range(9)]
U2_ana_T_10 = [U2(n/8, 1.0) for n in range(9)]
U2_ana_T_20 = [U2(n/8, 2.0) for n in range(9)]

U2_num_grid, _, T_num_Values = U2_num(1/2, 8)

U2_num_T_00625 = U2_num_grid[T_num_Values.index(0.0625)]
U2_num_T_0125 = U2_num_grid[T_num_Values.index(0.125)]
U2_num_T_025 = U2_num_grid[T_num_Values.index(0.25)]
U2_num_T_05 = U2_num_grid[T_num_Values.index(0.5)]
U2_num_T_10 = U2_num_grid[T_num_Values.index(1.0)]
U2_num_T_20 = U2_num_grid[T_num_Values.index(2.0)]


Table_00625 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_00625,
               "Numerical Solution":U2_num_T_00625,
               "Error":[abs(U2_ana_T_00625[i]-U2_num_T_00625[i]) for i in range(9)]}

Table_0125 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_0125,
               "Numerical Solution":U2_num_T_0125,
               "Error":[abs(U2_ana_T_0125[i]-U2_num_T_0125[i]) for i in range(9)]}

Table_025 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_025,
               "Numerical Solution":U2_num_T_025,
               "Error":[abs(U2_ana_T_025[i]-U2_num_T_025[i]) for i in range(9)]}

Table_05 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_05,
               "Numerical Solution":U2_num_T_05,
               "Error":[abs(U2_ana_T_05[i]-U2_num_T_05[i]) for i in range(9)]}

Table_10 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_10,
               "Numerical Solution":U2_num_T_10,
               "Error":[abs(U2_ana_T_10[i]-U2_num_T_10[i]) for i in range(9)]}

Table_20 = {"$X$":[n/8 for n in range(9)],
               "Analytical Solution":U2_ana_T_20,
               "Numerical Solution":U2_num_T_20,
               "Error":[abs(U2_ana_T_20[i]-U2_num_T_20[i]) for i in range(9)]}


draw_table(Table_0125, "table_0125", "|c||l|l|l|")
draw_table(Table_025, "table_025", "|c||l|l|l|")
draw_table(Table_05, "table_05", "|c||l|l|l|")
draw_table(Table_10, "table_10", "|c||l|l|l|")


fig = plt.figure(figsize=(13,9))

ax = fig.add_subplot(2, 3, 1)
ax.plot([n/8 for n in range(9)], Table_00625["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_00625["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=0.0625")

ax = fig.add_subplot(2, 3, 2)
ax.plot([n/8 for n in range(9)], Table_0125["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_0125["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=0.125")

ax = fig.add_subplot(2, 3, 3)
ax.plot([n/8 for n in range(9)], Table_025["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_025["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=0.25")

ax = fig.add_subplot(2, 3, 4)
ax.plot([n/8 for n in range(9)], Table_05["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_05["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=0.5")

ax = fig.add_subplot(2, 3, 5)
ax.plot([n/8 for n in range(9)], Table_10["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_10["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=1.0")

ax = fig.add_subplot(2, 3, 6)
ax.plot([n/8 for n in range(9)], Table_20["Analytical Solution"],label="Analytical Solution")
ax.plot([n/8 for n in range(9)], Table_20["Numerical Solution"],label="Numerical Solution")
ax.set_xlabel("X")
ax.set_ylabel(r"$U_2\left(x,t\right)$")
ax.legend()
ax.set_title("T=2.0")
plt.tight_layout()
plt.savefig("PT2_ii.pdf")

U2_num_1, X_Values_1, T_Values_1 = U2_num(1/12, 8) # C = 1/12, N = 8
U2_num_2, X_Values_2, T_Values_2 = U2_num(1/12, 16) # C = 1/12, N = 16
U2_num_3, X_Values_3, T_Values_3 = U2_num(1/12, 32) # C = 1/12, N = 32

max_error_1 = max([abs(U2(n/8,2)-U2_num_1[-1][n]) for n in range(9)])
max_error_2 = max([abs(U2(n/16,2)-U2_num_2[-1][n]) for n in range(17)])
max_error_3 = max([abs(U2(n/32,2)-U2_num_3[-1][n]) for n in range(33)])

U2_num_4, X_Values_4, T_Values_4 = U2_num(1/3, 8) # C = 1/3, N = 8
U2_num_5, X_Values_5, T_Values_5 = U2_num(1/3, 16) # C = 1/3, N = 16
U2_num_6, X_Values_6, T_Values_6 = U2_num(1/3, 32) # C = 1/3, N = 32

max_error_4 = max([abs(U2(n/8,2)-U2_num_4[-1][n]) for n in range(9)])
max_error_5 = max([abs(U2(n/16,2)-U2_num_5[-1][n]) for n in range(17)])
max_error_6 = max([abs(U2(n/32,2)-U2_num_6[-1][n]) for n in range(33)])

U2_num_7, X_Values_7, T_Values_7 = U2_num(1/6, 8) # C = 1/6, N = 8
U2_num_8, X_Values_8, T_Values_8 = U2_num(1/6, 16) # C = 1/6, N = 16
U2_num_9, X_Values_9, T_Values_9 = U2_num(1/6, 32) # C = 1/6, N = 32

max_error_7 = max([abs(U2(n/8,2)-U2_num_7[-1][n]) for n in range(9)])
max_error_8 = max([abs(U2(n/16,2)-U2_num_8[-1][n]) for n in range(17)])
max_error_9 = max([abs(U2(n/32,2)-U2_num_9[-1][n]) for n in range(33)])

Table_errors = {"C":[r"$\frac{1}{12}$",r"$\frac{1}{6}$",r"$\frac{1}{3}$"],
                r"$N=8$": [max_error_1, max_error_7, max_error_4],
                r"$N=16$": [max_error_2, max_error_8, max_error_5],
                r"$N=32$": [max_error_3, max_error_9, max_error_6],}

draw_table(Table_errors, "table_errors", "|c||l|l|l|")