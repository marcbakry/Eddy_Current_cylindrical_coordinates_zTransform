#############################################
# PLOT THE OBSERVABLE AS A FUNCTION OF TIME #
#############################################
import sys
import numpy as np
import matplotlib.pyplot as plt

# get data
data_file = sys.argv[1]
data      = np.loadtxt(data_file,skiprows=1,delimiter=";")

time = data[:,0]
A    = data[:,1]
Br   = data[:,2]
Bz   = data[:,3]

xlim = (np.min(time),np.max(time))

# plot A
plt.figure()
plt.plot(time,A,"r-")
plt.xlim(xlim)
plt.axhline(color="k")
plt.axvline(color="k")
plt.grid(True,which="both")
plt.xlabel("Time (s)")
plt.ylabel("A")
plt.title("Potential A")

# plot Br
plt.figure()
plt.plot(time,Br,"r-")
plt.xlim(xlim)
plt.axhline(color="k")
plt.axvline(color="k")
plt.grid(True,which="both")
plt.xlabel("Time (s)")
plt.ylabel("Br")
plt.title("Radial component Br")

# plot Bz
plt.figure()
plt.plot(time,Bz,"r-")
plt.xlim(xlim)
plt.axhline(color="k")
plt.axvline(color="k")
plt.grid(True,which="both")
plt.xlabel("Time (s)")
plt.ylabel("Bz")
plt.title("Radial component Bz")

# plot
plt.show()