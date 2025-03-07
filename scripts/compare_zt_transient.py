import numpy as np
import scipy.linalg as sp
import matplotlib.pyplot as plt

data_transient = np.loadtxt("../output/fields_nt300_n0.csv",delimiter=";",skiprows=1)
data_ztransfor = np.loadtxt("../output/fields_nt300_nz1200_n0.csv",delimiter=";",skiprows=1)

t_tr = data_transient[:,0]
t_zt = data_ztransfor[:,0]

A_tr = data_transient[:,1]
A_zt = data_ztransfor[:,1]

err = sp.norm(A_tr-A_zt)/sp.norm(A_tr)

plt.figure()
plt.plot(t_tr,A_tr,"g-")
plt.plot(t_zt,A_zt,"r--")
plt.grid(True,which="both")
plt.axhline()
plt.axvline()
plt.title("error = " + str(err))
plt.show()