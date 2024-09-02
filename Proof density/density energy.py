import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def sortAccordingToA(A, *args):
    arg = np.argsort(A)
    # [:] allows to change by pointer the arrays (and not overwrite them), so no need to return them
    A[:] = A[arg]
    for arr in args:
        arr[:] = arr[arg]


def downsample(x, n=3, m=3):
    
    x_padded = np.concatenate((x, x[:n - 1]))

    weights = np.ones(n)/n
    moving_avg = np.convolve(x_padded, weights, mode='valid')
    return moving_avg[:len(x)][::m]

fontsizeMain = 23
fontsize = 13
plt.rcParams.update({
    'font.size' : fontsizeMain,                   # Set font size to 11pt
    'axes.labelsize': fontsizeMain,               # -> axis labels
    'legend.fontsize': fontsizeMain,              # -> legends
})    
markersize = 10


n = 3
m = 2
        
color = ["#02ecee","#fa03ee", "#00fe01","#fa8173","black","#ffd41c"]
marker = ["^", "D","s", "o", "v", "<"]
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

A = np.load("proofq4.npz")
q4 = A["q4"]


sortAccordingToA(phi, density, Es, Eb, q4)


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
#color = [plt.get_cmap("cool")(i/(N-1/N)) for i in range(N)]
# Create figure and axes
fig, axs = plt.subplots(2, 1, figsize=(6, 5), sharex=True, layout = "constrained")


for i in range(N):
    axs[0].plot(L[i][::m], downsample(density[i], n = n, m = m), color=color[i], marker = marker[i], markersize = markersize, mfc = "none", linestyle = "--")
axs[0].set_ylabel('$\phi$')

for i in range(N):
    axs[1].plot(L[i][::m], downsample(q4[i], n = n, m = m), color=color[i], marker = marker[i], markersize = markersize, mfc = "none", linestyle = "--")
axs[1].set_ylabel('$q_4$')
axs[1].set_xlabel(r'$x/L_x$')
plt.tight_layout()
plt.subplots_adjust(hspace=0)
fig.savefig("EDMD_density_q4_profiles.pdf")

fig, axs = plt.subplots(2, 1, figsize=(6, 5), sharex=True, layout = "constrained")
for i in range(N):
    axs[0].plot(L[i][::m], downsample(Es[i], n = n, m = m), color=color[i], marker = marker[i], markersize = markersize, mfc = "none", linestyle = "--")
axs[0].set_ylabel('$E_s/T$')

# Third subplot
for i in range(N):
    axs[1].plot(L[i][::m], downsample(Eb[i], n = n, m = m), color=color[i], marker = marker[i], markersize = markersize, mfc = "none", linestyle = "--")
axs[1].set_xlabel(r'$x/L_x$')
axs[1].set_ylabel('$E_b/T$')


fig.savefig("EDMD_energy_profile.pdf")
plt.show()
