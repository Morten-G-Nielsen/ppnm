import numpy as np
from scipy.integrate import quad

f1 = lambda x: 1/np.sqrt(x)
f2 = lambda x: np.log(x)/np.sqrt(x)
f3 = lambda x: np.exp(-x**2)
f4 = lambda x: 1/(1+x**2)

q, err, info1 = quad(f1, 0, 1, full_output=True)
q, err, info2 = quad(f2, 0, 1, full_output=True)
q, err, info3 = quad(f3, -np.inf, np.inf, full_output=True)
q, err, info4 = quad(f4, 0, np.inf, full_output=True)

print("---Integration using QUADPCK from scipy---")
print("Evals for 1/sqrt(x) from 0 to 1:      ", info1['neval'])
print("Evals for ln(x)/sqrt(x) from 0 to 1:  ", info2['neval'])
print("Evals for exp(-x^2) from -inf to inf: ", info3['neval'])
print("Evals for 1/(1+x^2) from 0 to inf:    ", info4['neval'])
