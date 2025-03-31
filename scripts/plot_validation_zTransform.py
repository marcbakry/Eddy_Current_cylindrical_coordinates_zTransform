import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../output/validation_zTransform.csv",delimiter=";")
t = data[:,0]
signal = data[:,1]
signal_backward = data[:,2]

plt.figure()
plt.plot(t,signal,"g-")
plt.plot(t,signal_backward,"r--")
plt.axhline()
plt.axvline()
plt.grid(True,which="both")
plt.show()