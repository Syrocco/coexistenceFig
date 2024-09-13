import numpy as np
import matplotlib.pyplot as plt

fontsizeMain = 33
fontsize = 18
plt.rcParams.update({
    'font.size' : fontsizeMain,                   # Set font size to 11pt
    'axes.labelsize': fontsizeMain,               # -> axis labels
    'legend.fontsize': fontsizeMain,              # -> legends
}) 
A = np.load("data.npz")
phi = A["phi"]
Et = A["et"]
Em  = A["em"]
arg = np.argsort(phi)
phi = phi[arg]
Et = Et[arg]
Em = Em[arg]

arg = phi > 0.8
phir = phi[arg]
Emr = Em[arg]
Etr = Et[arg]

Tb = 1.367
c1, c2 = "#fa03ee", "#00fe01"
S = 150
#######################


fig = plt.figure(layout="constrained", figsize=(11, 7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.scatter(phi, Et/Tb, label="Theory", facecolor=c1, marker="x", s=S)
ax1.scatter(phi, Em/Tb, label="Measured", facecolor="none", edgecolor=c2, marker="o", s=S)
ax1.set_xlabel(r"$\phi$")
ax1.set_ylabel(r"$T/T_b$")
ax1.legend(loc='lower left', frameon=False, fontsize = fontsizeMain*0.8)
ax1.axvspan(0.8476, 0.888, color='gray', alpha=0.2)
ax1.set_ylim(0, 1)
ax1.set_xlim(0, np.max(phi)*1.05)

inset_ax = ax1.inset_axes([0.55, 0.65, 0.4, 0.3])
inset_ax.scatter(phir, Etr/Tb,  facecolor=c1, marker="x", s=S)
inset_ax.scatter(phir, Emr/Tb, facecolor="none", edgecolor=c2, marker="o", s=S)
inset_ax.axvspan(0.8476, 0.888, color='gray', alpha=0.2)
inset_ax.set_xlabel(r"$\phi$")
inset_ax.set_ylabel(r"$T/T_b$")

ax2.scatter(phi, Em/Et, facecolor="none", edgecolor="#fa8173", marker="o", s=S)
ax2.set_xlabel(r"$\phi$")
ax2.set_ylabel(r"$T_\text{measured}/T_{\text{theory}}$")
ax2.axvspan(0.8476, 0.888, color='gray', alpha=0.2)
ax2.set_xlim(0, np.max(phi)*1.05)
ax2.set_ylim(1, None)

ax1.text(0.05, 0.88, 'a)', transform=ax1.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax2.text(0.05, 0.88, 'b)', transform=ax2.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')

fig.savefig("theo vs measuresement.pdf")