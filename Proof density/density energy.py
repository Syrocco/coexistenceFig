import numpy as np
import matplotlib.pyplot as plt


def sortAccordingToA(A, *args):
    arg = np.argsort(A)
    # [:] allows to change by pointer the arrays (and not overwrite them), so no need to return them
    A[:] = A[arg]
    for arr in args:
        arr[:] = arr[arg]
 
        
Lx = 23.690
T = 1.367

A = np.load("proof.npz")
LL = A["LL"]
LX = A["LX"]
L = [L/Lx for L, Lx in zip(LL, LX)]
phi = A["phi"]
density = A["dens"]
Es = A["ES"]/T
Eb = A["EB"]/T


sortAccordingToA(phi, density, Es, Eb)


# Fixing stupid mistake in the original npz
A = np.load("small.npz")
density[-1] = A["dens"][0]
Es[-1] = (A["ES"]/T)[0]
Eb[-1] = (A["EB"]/T)[0]

A = np.load("small2.npz")
density[0] = A["dens"][0]
Es[0] = (A["ES"]/T)[0]
Eb[0] = (A["EB"]/T)[0]

N = len(Es)
color = [plt.get_cmap("cool")(i/(N-1/N)) for i in range(N)]
# Create figure and axes
fig, axs = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

# First subplot
for i in range(N):
    axs[0].plot(L[i], density[i], color=color[i])
axs[0].set_ylabel('$\phi$')

# Second subplot
for i in range(N):
    axs[1].plot(L[i], Es[i], color=color[i])
axs[1].set_ylabel('$E_s/T$')

# Third subplot
for i in range(N):
    axs[2].plot(L[i], Eb[i], color=color[i])
axs[2].set_xlabel(r'$x/L_x$')
axs[2].set_ylabel('$E_b/T$')

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(hspace=0)

# Show plot
plt.show()
