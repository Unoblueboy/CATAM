# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 04:14:02 2018

@author: Natha
"""

from numpy import e, log, sin, cos
import matplotlib.pyplot as plt
import pandas as pd

def euler_method(x0, xN, N, func, Y0):
    h = (xN-x0)/N
    Y = [Y0]
    for n in range(1, N+1):
        Yn_1 = Y[-1]
        xn_1 = x0 + (n-1)*h
        Yn = Yn_1 + h * func(xn_1, Yn_1)
        Y.append(Yn)
    return Y

def leapfrog_method(x0, xN, N, func, Y0):
    h = (xN-x0)/N
    Y = [Y0, Y0 + h*func(x0, Y0)]
    for n in range(2, N+1):
        Yn_1 = Y[-1]
        Yn_2 = Y[-2]
        xn_1 = x0 + (n-1)*h
        Y.append(Yn_2+2*h*func(xn_1, Yn_1))
    return Y

def rk4_method(x0, xN, N, func, Y0):
    h = (xN-x0)/N
    Y = [Y0]
    for n in range(1, N+1):
        Yn_1 = Y[-1]
        xn_1 = x0 + (n-1)*h
        k1 = h*func(xn_1, Yn_1)
        k2 = h*func(xn_1+0.5*h, Yn_1+0.5*k1)
        k3 = h*func(xn_1+0.5*h, Yn_1+0.5*k2)
        k4 = h*func(xn_1+h, Yn_1+k3)
        Yn = Yn_1+ 1/6*(k1+2*k2+2*k3+k4)
        Y.append(Yn)
    return Y

def x_axis_gen(x0, xN, n):
    h = (xN-x0)/n
    X = [x0]
    for i in range(1, n+1):
        Xn = X[-1]
        X.append(Xn+h)
    return X
    




f = lambda x, y: -4*y + 3*e**-x
y_ana = lambda x: e**-x - e**(-4*x)


def latex_format(num):
    parts = "{0:.3e}".format(num).split("e")
    if len(parts) == 1:
        part = parts[0]
        part = part.replace(r"inf", r"\infty")
        return "$"+part+"$"
    base = parts[0]
    exp = parts[1]
    latex_str = r"${0} \times 10^{{{1}}}$".format(str(float(base)), str(int(exp)))
    return latex_str

def Tabulate_leapfrog(n, table_length):
    Yl_values = leapfrog_method(0, 10, n, f, 0)
    X_values = x_axis_gen(0, 10, n)
    tick = round(n/table_length)
    
    ln_E_n = []
    
    table = {"$n$":[],
             "$x_n$":[],
             "$Y_n$":[],
             "$y(x_n)$":[],
             "$E_n$":[],
             "$ln(|E_n|)$":[]}
    for i in range(len(X_values)):
        if i%tick!=0:
            continue
        x_n = X_values[i]
        Y_n = Yl_values[i]
        y_n = y_ana(x_n)
        E_n = Y_n - y_n
        ln_mag_E_n = log(abs(E_n))
        table["$n$"].append(i)
        table["$x_n$"].append(x_n)
        table["$Y_n$"].append(latex_format(Y_n))
        table["$y(x_n)$"].append(y_n)
        table["$E_n$"].append(latex_format(E_n))
        table["$ln(|E_n|)$"].append(ln_mag_E_n)
        ln_E_n.append(ln_mag_E_n)
    
    ln_E_N = ln_E_n[-1]
    ln_E_0 = ln_E_n[1]
    x_N = float(table["$x_n$"][-1])
    x_0 = float(table["$x_n$"][1])
    growth_rate = (ln_E_N-ln_E_0)/(x_N-x_0)
    return table, growth_rate


# can see E against ln_mag_E has gradient approx 3

table1, gr1 = Tabulate_leapfrog(25,25)
table2, gr2 = Tabulate_leapfrog(50,10)
table3, gr3 = Tabulate_leapfrog(100,10)
table4, gr4 = Tabulate_leapfrog(200,10)

def draw_table(table, index, file_name, column_format):
    template = r'''\documentclass[preview]{{standalone}}
    \usepackage{{booktabs}}
    \begin{{document}}
    {}
    \end{{document}}
    '''
    dataframe = pd.DataFrame.from_dict(table).set_index(index)
    with open(file_name+".tex", 'w') as file:
        file.write(template.format(dataframe.to_latex(
                escape=False,
                column_format = column_format)))



draw_table(table1, "$n$", "table1", "|c||c|r|c|r|r|")
draw_table(table2, "$n$", "table2", "|c||c|r|c|r|r|")
draw_table(table3, "$n$", "table3", "|c||c|r|c|r|r|")
draw_table(table4, "$n$", "table4", "|c||c|r|c|r|r|")

print(gr1)
print(gr2)
print(gr3)
print(gr4)

# By using larger and larger n we see that the growth rate in fact increases
# but we start with a smaller error

# QUESTION 3

Yr_values = rk4_method(0, 4, 10, f, 0)
Ye_values = euler_method(0, 4, 10, f, 0)
X_values = x_axis_gen(0, 4, 10)
Xa_values = x_axis_gen(0, 4, 1000)
Ya_values = [y_ana(x) for x in Xa_values]

table5 = {"$x$": X_values,
        "Runge-Kutta" : Yr_values,
        "Euler": Ye_values,
        "Analytical": [Ya_values[i*100] for i in range(0,11)]}

draw_table(table5, "$x$", "table5", "|c||c|c|c|")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(X_values, Ye_values)
ax.plot(X_values, Yr_values)
ax.plot(Xa_values, Ya_values)
ax.legend(['Euler', 'Runge-Kutta', 'Analytical'], loc='upper right')

plt.show()

fig.savefig("Q3.pdf")
# QUESTION 4

y_ana_04 = y_ana(0.4)
euler = []
leapfrog = []
rk4 = []
h_vals = []
for k in range(0,16):
    n = 2**k
    h_vals.append(0.4/n)
    Ye_n = euler_method(0,0.4,n,f,0)[n]
    Yl_n = leapfrog_method(0,0.4,n,f,0)[n]
    Yr_n = rk4_method(0,0.4,n,f,0)[n]
    euler.append(Ye_n-y_ana_04)
    leapfrog.append(Yl_n-y_ana_04)
    rk4.append(Yr_n-y_ana_04)

table6 = {"$h$" : h_vals,
         "Euler" : euler,
         "Leapfrog" : leapfrog,
         "Runge-Kutta" : rk4}

draw_table(table6, "$h$", "table6", "|c||c|c|c|")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot([log(abs(x)) for x in h_vals], [log(abs(x)) for x in euler])
ax.plot([log(abs(x)) for x in h_vals], [log(abs(x)) for x in leapfrog])
ax.plot([log(abs(x)) for x in h_vals], [log(abs(x)) for x in rk4])
ax.legend(['Euler', 'Leapfrog', 'Runge-Kutta'], loc='lower right')

plt.show()

fig.savefig("Q4.pdf")

# Question 6

def coupled_rk4_method(t0, tN, N, gamma, delta, omega, a, w):
    f1 = lambda t, y1, y2: y2
    f2 = lambda t, y1, y2: -gamma*y2 - (delta**3)*(y1**2)*y2 - (omega**2)*y1+a*sin(w*t)
    h = (tN-t0)/N
    Y = [[0,0]]
    
    for n in range(1,N+1):
        tn_1 = t0 +(n-1)*h
        Y1n_1= Y[-1][0]
        Y2n_1= Y[-1][1]
        
        k11 = h*f1(tn_1, Y1n_1, Y2n_1)
        k12 = h*f2(tn_1, Y1n_1, Y2n_1)
        
        k21 = h*f1(tn_1+0.5*h, Y1n_1+0.5*k11, Y2n_1+0.5*k12)
        k22 = h*f2(tn_1+0.5*h, Y1n_1+0.5*k11, Y2n_1+0.5*k12)
        
        k31 = h*f1(tn_1+0.5*h, Y1n_1+0.5*k21, Y2n_1+0.5*k22)
        k32 = h*f2(tn_1+0.5*h, Y1n_1+0.5*k21, Y2n_1+0.5*k22)
        
        k41 = h*f1(tn_1+h, Y1n_1+k31, Y2n_1+k32)
        k42 = h*f2(tn_1+h, Y1n_1+k31, Y2n_1+k32)
        
        Y1n = Y1n_1 + 1/6*(k11+2*k21+2*k31+k41)
        Y2n = Y2n_1 + 1/6*(k12+2*k22+2*k32+k42)
        
        Y.append([Y1n, Y2n])
    
    soln = [y[0] for y in Y]
    return soln

#y_ana = lambda t: 1/7*((e**(-t/2))*((3**0.5)*cos((3**0.5)*t/2) + 5*sin((3**0.5)*t/2)) - 2 * sin((3**0.5)*t) - (3**0.5)*cos((3**0.5)*t))

def get_y_ana(gamma, w):
    A = w*gamma/((1-w**2)**2+(gamma*w)**2)
    B = w*(gamma**2 - 2*(1-w**2))/(((1-w**2)**2 +(gamma*w)**2)*(4-gamma**2)**0.5)
    C = (1-w**2)/((1-w**2)**2 +(gamma*w)**2)
    D = -(gamma*w)/((1-w**2)**2 +(gamma*w)**2)
    k = ((4-gamma**2)**0.5)/2
    
    y_ana = lambda t: A*(e**(-gamma*t/2))*cos(k*t) + B*(e**(-gamma*t/2))*sin(k*t) + C*sin(w*t) + D*cos(w*t)
    
    return y_ana

def Tabulate_rk4(t0, tN, N, gamma, w, num_output):
    Yr_values = coupled_rk4_method(t0, tN, N, gamma, 0, 1, 1, w)
    X_values = x_axis_gen(t0, tN, N)
    tick = round(N/num_output)
    
    y_ana = get_y_ana(gamma, w)
    
    table = {"$n$":[],
             "$x_n$":[],
             "$Y_n$":[],
             "$y(x_n)$":[],
             "$E_n$":[]}
    for i in range(N+1):
        if i % tick != 0:
            continue
        x_n = X_values[i]
        Y_n = Yr_values[i]
        y_n = y_ana(x_n)
        E_n = Y_n - y_n
        table["$n$"].append(i)
        table["$x_n$"].append(x_n)
        table["$Y_n$"].append(Y_n)
        table["$y(x_n)$"].append(y_n)
        table["$E_n$"].append(latex_format(E_n))
    return table

table7 = Tabulate_rk4(0, 10, 25, 1, 3**(1/2), 25)
table8 = Tabulate_rk4(0, 10, 50, 1, 3**(1/2), 10)
table9 = Tabulate_rk4(0, 10, 100, 1, 3**(1/2), 10)

draw_table(table7, "$n$", "table7", "|c||c|r|c|r|")
draw_table(table8, "$n$", "table8", "|c||c|r|c|r|")
draw_table(table9, "$n$", "table9", "|c||c|r|c|r|")

# QUESTION 7

X_ana = x_axis_gen(0, 40, 4000)

X_40 = x_axis_gen(0, 40, 100)
Y_1 = coupled_rk4_method(0, 40, 100, 0.25, 0, 1, 1, 1)
y_ana_1 = get_y_ana(0.25, 1)
Y_ana_1 = [y_ana_1(x) for x in X_ana]


Y_2 = coupled_rk4_method(0, 40, 100, 0.5, 0, 1, 1, 1)
y_ana_2 = get_y_ana(0.5 , 1)
Y_ana_2 = [y_ana_2(x) for x in X_ana]


Y_3 = coupled_rk4_method(0, 40, 100, 1, 0, 1, 1, 1)
y_ana_3 = get_y_ana(1   , 1)
Y_ana_3 = [y_ana_3(x) for x in X_ana]


Y_4 = coupled_rk4_method(0, 40, 100, 1.9, 0, 1, 1, 1)
y_ana_4 = get_y_ana(1.9 , 1)
Y_ana_4 = [y_ana_4(x) for x in X_ana]

fig, ax = plt.subplots(4, sharex=True)

ax[0].plot(X_40, Y_1)
ax[0].plot(X_ana, Y_ana_1)
ax[0].legend(['Runge-Kutta', 'Analytic'], loc='upper center', 
              bbox_to_anchor=(0.5, 1.5),
              ncol = 2, fancybox=True, shadow=True)
ax[0].set_title(r"$\gamma=0.25$", loc="right")

ax[1].plot(X_40, Y_2)
ax[1].plot(X_ana, Y_ana_2)
ax[1].set_title(r"$\gamma=0.5$", loc="right")

ax[2].plot(X_40, Y_3)
ax[2].plot(X_ana, Y_ana_3)
ax[2].set_title(r"$\gamma=1$", loc="right")

ax[3].plot(X_40, Y_4)
ax[3].plot(X_ana, Y_ana_4)
ax[3].set_title(r"$\gamma=1.9$", loc="right")

fig.subplots_adjust(hspace=0.6)

plt.show()

fig.savefig("Q7a.pdf")


Y_5 = coupled_rk4_method(0, 40, 100, 0.25, 0, 1, 1, 2)
y_ana_5 = get_y_ana(0.25, 2)
Y_ana_5 = [y_ana_5(x) for x in X_ana]

Y_6 = coupled_rk4_method(0, 40, 100, 0.5, 0, 1, 1, 2)
y_ana_6 = get_y_ana(0.5 , 2)
Y_ana_6 = [y_ana_6(x) for x in X_ana]

Y_7 = coupled_rk4_method(0, 40, 100, 1, 0, 1, 1, 2)
y_ana_7 = get_y_ana(1   , 2)
Y_ana_7 = [y_ana_7(x) for x in X_ana]

Y_8 = coupled_rk4_method(0, 40, 100, 1.9, 0, 1, 1, 2)
y_ana_8 = get_y_ana(1.9 , 2)
Y_ana_8 = [y_ana_8(x) for x in X_ana]

fig, ax = plt.subplots(4, sharex=True)

ax[0].plot(X_40, Y_5)
ax[0].plot(X_ana, Y_ana_5)
ax[0].legend(['Runge-Kutta', 'Analytic'], loc='upper center', 
              bbox_to_anchor=(0.5, 1.5),
              ncol = 2, fancybox=True, shadow=True)
ax[0].set_title(r"$\gamma=0.25$", loc="right")

ax[1].plot(X_40, Y_6)
ax[1].plot(X_ana, Y_ana_6)
ax[1].set_title(r"$\gamma=0.5$", loc="right")

ax[2].plot(X_40, Y_7)
ax[2].plot(X_ana, Y_ana_7)
ax[2].set_title(r"$\gamma=1$", loc="right")

ax[3].plot(X_40, Y_8)
ax[3].plot(X_ana, Y_ana_8)
ax[3].set_title(r"$\gamma=1.9$", loc="right")

fig.subplots_adjust(hspace=0.6)

plt.show()

fig.savefig("Q7b.pdf")

# QUESTION 8

X_60 = x_axis_gen(0, 60, 1000)
Y_9 = coupled_rk4_method(0,60,1000,0,0.25,1,1,1)

Y_10 = coupled_rk4_method(0,60,1000,0,0.5,1,1,1)

Y_11 = coupled_rk4_method(0,60,1000,0,1,1,1,1)


X_60_1 = x_axis_gen(0, 60, 30000)
Y_12 = coupled_rk4_method(0,60,30000,0,20,1,1,1)

fig, ax = plt.subplots(4, sharex=True)

ax[0].plot(X_60, Y_9)
ax[0].set_title(r"$\delta=0.25$", loc="right")

ax[1].plot(X_60, Y_10)
ax[1].set_title(r"$\delta=0.5$", loc="right")

ax[2].plot(X_60, Y_11)
ax[2].set_title(r"$\delta=1$", loc="right")

ax[3].plot(X_60_1, Y_12)
ax[3].set_title(r"$\delta=20$", loc="right")

fig.subplots_adjust(hspace=0.6)

plt.show()

fig.savefig("Q8.pdf")