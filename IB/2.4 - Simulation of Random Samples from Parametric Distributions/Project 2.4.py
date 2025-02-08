# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:29:15 2019

@author: Natha
"""

from random import uniform
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.warnings.filterwarnings('ignore')

# Question 2

def exp_sample(n, t0):
    samples = []
    for _ in range(n):
        samples.append(-np.log(1-uniform(0,1))/t0)
    return samples

def log_like_median(samples):
    def l(m):
        result = 1
        for xi in samples:
            result *= np.log(2)/m*2**(-xi/m)
        return np.log(result)
    return l

samples = exp_sample(6, 1.2)
samples.sort()
l6 = log_like_median(samples)
x_values = np.linspace(0,5,100)
y_values_6 = l6(x_values)

median_MLE = np.log(2)*sum(samples)/6
with open('data.txt', 'w') as file:
    file.write("Question 2\n")
    file.write("The true median is: " + str(np.log(2)/1.2) +"\n")
    file.write("The median MLE is: " + str(median_MLE)+"\n")
    file.write("The samples:" + str(samples)+"\n\n")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(x_values, y_values_6)
plt.show()

fig.savefig("Q2.pdf")

# Question 3
# n=25

samples_25 = exp_sample(25, 1.2)
samples_25.sort()
l25 = log_like_median(samples_25)
x_values = np.linspace(0,5,100)
y_values_25 = l25(x_values)

median_MLE_25 = np.log(2)*sum(samples_25)/25
with open('data.txt', 'a') as file:
    file.write("Question 3\n")
    file.write("n = 25\n")
    file.write("The median MLE is: " + str(median_MLE_25) + "\n")
    file.write("The samples:" + str(samples_25) + "\n\n")

# n=50

samples_50 = exp_sample(50, 1.2)
samples_50.sort()
l50 = log_like_median(samples_50)
x_values = np.linspace(0,5,100)
y_values_50 = l50(x_values)

median_MLE_50 = np.log(2)*sum(samples_50)/50
with open('data.txt', 'a') as file:
    file.write("n = 50\n")
    file.write("The median MLE is: " + str(median_MLE_50) + "\n")
    file.write("The samples:" + str(samples_50) + "\n\n")

# n=100

samples_100 = exp_sample(100, 1.2)
samples_100.sort()
l100 = log_like_median(samples_100)
x_values = np.linspace(0,5,100)
y_values_100 = l100(x_values)

median_MLE_100 = np.log(2)*sum(samples_100)/100
with open('data.txt', 'a') as file:
    file.write("n = 100\n")
    file.write("The median MLE is: " + str(median_MLE_100) + "\n")
    file.write("The samples:" + str(samples_100) + "\n\n")

# The Graph

fig = plt.figure(figsize=(13,5))
ax = fig.add_subplot(1, 3, 1)
ax.plot(x_values, y_values_25)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_25(m)$")
ax.set_title("$n=25$")

ax = fig.add_subplot(1, 3, 2)
ax.plot(x_values, y_values_50)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_50(m)$")
ax.set_title("$n=50$")


ax = fig.add_subplot(1, 3, 3)
ax.plot(x_values, y_values_100)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_100(m)$")
ax.set_title("$n=100$")
plt.tight_layout()

fig.savefig("Q3.pdf")


# Question 7

def gamma_sample(n, t0):
    samples = []
    for _ in range(n):
        samples.append(sum(exp_sample(2, t0)))
    return samples

def log_like_theta(samples):
    def l(t):
        result = 1
        for xi in samples:
            result *= t*t*xi*np.e**(-t*xi)
        return np.log(result)
    return l

samples_dist_10 = gamma_sample(10, 2.2)
samples_dist_10.sort()
l_dist_10 = log_like_theta(samples_dist_10)
x_values = np.linspace(0,5,100)
y_values_dist_10 = l_dist_10(x_values)

theta_MLE_dist_10 = 2/(sum(samples_dist_10)/10)
with open('data.txt', 'a') as file:
    file.write("Question 7\n")
    file.write("n = 10\n")
    file.write("The Theta MLE is: " + str(theta_MLE_dist_10) + "\n")
    file.write("The samples:" + str(samples_dist_10) + "\n\n")

# n = 30

samples_dist_30 = gamma_sample(30, 2.2)
samples_dist_30.sort()
l_dist_30 = log_like_theta(samples_dist_30)
x_values = np.linspace(0,5,100)
y_values_dist_30 = l_dist_30(x_values)

theta_MLE_dist_30 = 2/(sum(samples_dist_30)/30)
with open('data.txt', 'a') as file:
    file.write("n = 30\n")
    file.write("The Theta MLE is: " + str(theta_MLE_dist_30) + "\n")
    file.write("The samples:" + str(samples_dist_30) + "\n\n")



# n = 50

samples_dist_50 = gamma_sample(50, 2.2)
samples_dist_50.sort()
l_dist_50 = log_like_theta(samples_dist_50)
x_values = np.linspace(0,5,100)
y_values_dist_50 = l_dist_50(x_values)

theta_MLE_dist_50 = 2/(sum(samples_dist_50)/50)
with open('data.txt', 'a') as file:
    file.write("n = 50\n")
    file.write("The Theta MLE is: " + str(theta_MLE_dist_50) + "\n")
    file.write("The samples:" + str(samples_dist_50) + "\n\n")

# The Graph

fig = plt.figure(figsize = (13,5))
ax = fig.add_subplot(1, 3, 1)
ax.plot(x_values, y_values_dist_10)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_10(m)$")
ax.set_title("$n=10$")

ax = fig.add_subplot(1, 3, 2)
ax.plot(x_values, y_values_dist_30)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_30(m)$")
ax.set_title("$n=30$")

ax = fig.add_subplot(1, 3, 3)
ax.plot(x_values, y_values_dist_50)
ax.set_xlabel("$m$")
ax.set_ylabel("$l_50(m)$")
ax.set_title("$n=50$")
plt.tight_layout()

fig.savefig("Q7.pdf")


# Question 8

n=10
N=200
set_of_samples = [gamma_sample(n, 2.2) for _ in range(N)]
set_of_MLEs_10 = [2/(sum(x)/n) for x in set_of_samples]

n=50

set_of_samples = [gamma_sample(n, 2.2) for _ in range(N)]
set_of_MLEs_50 = [2/(sum(x)/n) for x in set_of_samples]

fig = plt.figure(figsize=(9,5))

ax = fig.add_subplot(1, 2, 1)
ax.hist(set_of_MLEs_10, bins=5)


ax = fig.add_subplot(1, 2, 2)
ax.hist(set_of_MLEs_50, bins=10)

fig.savefig("Q8.pdf")


# Question 10

def norm_sample(mu, n):
    samples = []
    while len(samples)<n:
        A = uniform(0,1)
        B = uniform(0,1)
        Phi = 2*np.pi*A
        V = -2*np.log(1-B)
        X = mu + V**(0.5)*np.cos(Phi)
        Y = mu + V**(0.5)*np.sin(Phi)
        samples.append(X)
        samples.append(Y)
    return samples[0:n]

# Question 11


sample = norm_sample(0,100)
x_bar = sum(sample)/100
CI = [x_bar - 1.2815516/10, x_bar + 1.2815516/10]
if CI[0]<0 and CI[1]>0:
    print("Sample Mean within CI")
else:
    print("Sample Mean not within CI")
with open('data.txt', 'a') as file:
    file.write("Question 11\n")
    file.write("sample mean = " + str(x_bar) + "\n")
    file.write("The CI is: " + str(CI) + "\n")
    file.write("Does lie in CI:" + str(CI[0]<0 and CI[1]>0) + "\n\n")



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

table_25 = {"sample mean":[],
         "lower bound":[],
         "upper bound":[],
         "in confidence interval?":[]}

for _ in range(25):
    s = norm_sample(0,100)
    m = sum(s)/100
    lb = m - 1.2815516/10
    ub = m + 1.2815516/10
    in_ci = int(lb<0 and ub>0)
    table_25["sample mean"].append(m)
    table_25["lower bound"].append(lb)
    table_25["upper bound"].append(ub)
    table_25["in confidence interval?"].append(in_ci)


draw_table(table_25, "table_25", "|c||l|l|l|c|")

with open('data.txt', 'a') as file:
    file.write("n = 25\n")
    file.write("mu = 0\n")
    file.write("# mu not in CI: " + str(25-sum(table_25["in"])) + "\n\n")
               

# Question 13

def chi_sample(DoF, n):
    samples = []
    for _ in range(n):
        norm_samples = norm_sample(0, DoF)
        sqr = [x*x for x in norm_samples]
        samples.append(sum(sqr))
    return samples

chi_samples_1_100 = chi_sample(1, 100)
chi_samples_1_300 = chi_sample(1, 300)
chi_samples_1_500 = chi_sample(1, 500)

chi_samples_5_100 = chi_sample(5, 100)
chi_samples_5_300 = chi_sample(5, 300)
chi_samples_5_500 = chi_sample(5, 500)

chi_samples_40_100 = chi_sample(40, 100)
chi_samples_40_300 = chi_sample(40, 300)
chi_samples_40_500 = chi_sample(40, 500)



fig = plt.figure()
ax11 = fig.add_subplot(3, 3, 1)
ax11.hist(chi_samples_1_100, bins=10)
ax11.set_title("1 DoF, n = 100")

ax12 = fig.add_subplot(3, 3, 2)
ax12.hist(chi_samples_5_100, bins=10)
ax12.set_title("5 DoF, n = 100")

ax13 = fig.add_subplot(3, 3, 3)
ax13.hist(chi_samples_40_100, bins=10)
ax13.set_title("40 DoF, n = 100")

ax21 = fig.add_subplot(3, 3, 4)
ax21.hist(chi_samples_1_300, bins=10)
ax21.set_title("1 DoF, n = 300")

ax22 = fig.add_subplot(3, 3, 5)
ax22.hist(chi_samples_5_300, bins=10)
ax22.set_title("5 DoF, n = 300")

ax23 = fig.add_subplot(3, 3, 6)
ax23.hist(chi_samples_40_300, bins=10)
ax23.set_title("40 DoF, n = 300")

ax31 = fig.add_subplot(3, 3, 7)
ax31.hist(chi_samples_1_500, bins=10)
ax31.set_title("1 DoF, n = 500")

ax32 = fig.add_subplot(3, 3, 8)
ax32.hist(chi_samples_5_500, bins=10)
ax32.set_title("5 DoF, n = 500")

ax33 = fig.add_subplot(3, 3, 9)
ax33.hist(chi_samples_40_500, bins=10)
ax33.set_title("40 DoF, n = 500")

plt.tight_layout()
plt.show()

fig.savefig("Q13.pdf")