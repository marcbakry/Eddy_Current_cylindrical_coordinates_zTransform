import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../output/quadrature_nodes.csv",delimiter=";")

plt.figure()
plt.plot(data[:,0],data[:,1],"k*")
plt.axhline()
plt.axvline()
plt.grid(True,which="both")
plt.show()