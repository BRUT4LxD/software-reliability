import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import math
import time

file_path = input('Podaj nazwe pliku z rozszerzeniem .csv znajdującego z w bieżącym folderze: ')
t = pd.read_csv(file_path, header=None).values.flatten()
alpha = float(input('Podaj dokładność: '))

start = time.time()

range_t = np.arange(len(t))
b = np.sum(range_t*t)

def j_m(x):
    return np.sum(1/(x - range_t)) - len(t)*np.sum(t)/(x*np.sum(t) - b)

i = len(t) 
N_jm  = fsolve(j_m, i, xtol=alpha)
while not (N_jm  > len(t) and N_jm  < 1000):
    N_jm  = fsolve(j_m, i, xtol=alpha)
    i+=1

phi_jm = len(t)/(N_jm*np.sum(t) - b)
T_241_jm = 1/(phi_jm*(N_jm - len(t) + 1))

end_jm = time.time() - start

start = time.time()

T = np.sum(t*t)
def s_w(x):
    return 2*len(t)/(x*T - np.sum(range_t*(t*t))) - 2*np.sum(1/((x - range_t)*T))

i = len(t) 
N_sw = fsolve(s_w, i, xtol=alpha)
while not (N_sw > len(t) and N_sw < 1000):
    N_sw = fsolve(s_w, i, xtol=alpha)
    i+=1

phi_sw = 2*np.sum(1/((N_sw - range_t)*T))
T_241_sw = np.sqrt(math.pi/(2*phi_sw*(N_sw - len(t) + 1)))   

end_sw = time.time() - start

np.set_printoptions(precision=8)
print("Model Jelińskiego-Morandy: \n")
print('N == ', N_jm)
print('Phi == ', phi_jm)
print('T_241 == ', T_241_jm)
print('Czas [ms] == ', end_jm*1000)
np.set_printoptions(precision=8)
print("\nModel Schicka-Wolvertona: \n")
print('N == ', N_sw)
print('Phi == ', phi_sw)
print('T_241 == ', T_241_sw)
print('Czas [ms] == ', end_sw*1000)

input('Press any key to conitunue...')