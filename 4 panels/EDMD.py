import numpy as np
import matplotlib.pyplot as plt
from helper import nanify, getY, getSolidAndLiquid,  getSolidFromLattice, getSolidAndLiquid2, Ncolors, postProcess
from matplotlib.colors import LogNorm


fig, ax = plt.subplots(1, 3, figsize=(15, 6), layout = "constrained")

fontsizeMain = 23
fontsize = 13
plt.rcParams.update({
    'font.size' : fontsizeMain,                   # Set font size to 11pt
    'axes.labelsize': fontsizeMain,               # -> axis labels
    'legend.fontsize': fontsizeMain,              # -> legends
})    
markersize = 10
c1 = "#00fe01"
c2 = "#fa03ee"
markersize = 100

epsilon = 0.2**2
if 1:
    
    A = np.load("varying res.npz")



    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])       


    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = False
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    color = Ncolors(7)
    so2 = []
    li3 = []
    so3 = []
    i = 0 # value of T
    for j in range(8, len(RES)):

        X3, Y3 = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], False)
        if fromLattice:
            so2.append(getSolidFromLattice(getY(L, XY[j, i], False)))
        li3.append(X3)
        so3.append(Y3)
    if plotCoexistenceOrWidth:
        if fromLattice:
            ax[0].scatter(so2, RES[8:], color=color[0], marker="s")
        ax[0].scatter(so3, RES[8:], color=c2, marker="o", s = markersize, facecolor = "none")
        ax[0].scatter(li3, RES[8:], color=c2, marker="s", s = markersize, facecolor = "none")
        #ax[0].legend()
        ax[0].set_xlabel(r"$\phi$")
        ax[0].set_ylabel(r"$\alpha$")
    else:
        ax[0].scatter(T, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
        #ax[0].legend()
        ax[0].set_xlabel(r"$T$")
        ax[0].set_ylabel(r"$\phi_s-\phi_l$")

if 1:
    A = np.load("arrays one res multiple temp.npz")
    E = np.load('arrays one res multiple temp E.npz')
    Eb = E["EB"]
    Es = E["ES"]
    E = E["E"]
    e = np.mean(E[0], axis = 1)
    es = np.nanmean(Es[0], axis = 1)
    eb = np.nanmean(Eb[0], axis = 1)
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])
    


        
    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = False
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    color = Ncolors(7)
    t = T/epsilon
    norm = LogNorm(vmin=np.min(t), vmax=np.max(t))
    cmap = plt.get_cmap('cool') 
    for j in range(1):
        so2 = []
        li3 = []
        so3= []
        for i in range(len(T)):
            
            X3, Y3 = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], False)
            if fromLattice:
                so2.append(getSolidFromLattice(getY(L, XY[j, i], False)))
            li3.append(X3)
            so3.append(Y3)
        if plotCoexistenceOrWidth:
            if fromLattice:
                ax[1].scatter(so2, T, color = color[j], marker = "s")
            sc_so3 = ax[1].scatter(so3, e/T, c= t, edgecolors=cmap(norm(t)), marker="o", s=markersize/1.5, cmap = cmap, norm=LogNorm())
            sc_li3 = ax[1].scatter(li3, e/T, c= t, edgecolors=cmap(norm(t)), marker="s", s=markersize/1.5, cmap = cmap, norm=LogNorm())
            sc_so3.set_facecolor('none')
            sc_li3.set_facecolor('none')
            
                        
            #ax[1].legend()
            ax[1].set_xlabel(r"$\phi$")
            ax[1].set_ylabel(r"$T/T_{b}$")
            ax[1].scatter([0.823], [1], marker = "*", color = "yellow", s = 200, edgecolors = "black")
            ax[1].scatter([0.87], [1], marker = "*", color = "yellow", s = 200, edgecolors = "black")
            ax[1].hlines(1, 0.80, 0.92, linestyle ="--", color = "gray", zorder = -20)
            ax[1].set_xlim(0.81, 0.9)
            ax[1].set_ylim(0, None)
            cbar = plt.colorbar(sc_so3, ax=ax[1], location='right', shrink=0.95, pad=-0.02)
            cbar.ax.set_title(r'$T_b/\varepsilon$', pad=10, fontsize = fontsizeMain*0.75) 
            
        else:
            ### Plot the the size of the region with respect to T instead of the coexistence vs T.
            ax[1].scatter(T/e, np.array(so3)-np.array(li3) , label = fr"$\alpha = ${RES[j]}")
            #ax[1].legend()
            ax[1].set_xlabel(r"$T$")
            ax[1].set_ylabel(r"$\phi_s-\phi_l$")


if 0:
    A = np.load("arraysCleaner.npz")
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])





    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = True
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    color = Ncolors(7)
    for j in range(len(RES) - 1):

        so2 = []
        li3 = []
        so3 = []
        for i in range(len(T)):
            X3, Y3 = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], False)
            if fromLattice:
                so2.append(getSolidFromLattice(getY(L, XY[j, i], False)))
            li3.append(X3)
            so3.append(Y3)
        if plotCoexistenceOrWidth:
            if fromLattice:
                ax[1, 0].scatter(so2, T, color=color[j+1], marker="s")
            ax[1, 0].scatter(so3, T, color=color[j+1], marker="v")
            ax[1, 0].scatter(li3, T, color=color[j+1], marker="v", label=fr"$\alpha = ${RES[j]}")
            ax[1, 0].legend()
            ax[1, 0].set_xlabel(r"$\phi$")
            ax[1, 0].set_ylabel("$T_{thermo}/\varepsilon$")
        else:
            ax[1, 0].scatter(T, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
            ax[1, 0].legend()
            ax[1, 0].set_xlabel(r"$T$")
            ax[1, 0].set_ylabel(r"$\phi_s-\phi_l$")
            
           
    A = np.load("Equilibrium.npz")
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])
    postProcess(XY, YY, liquid, solid)
   
   
   
   
    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = True
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    color = Ncolors(7)
    for j in range(len(RES)):
   
        so2 = []
        li3 = []
        so3 = []
        for i in range(len(T)):
            X3, Y3 = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], False)
            if fromLattice:
                so2.append(getSolidFromLattice(getY(L, XY[j, i], False)))
            li3.append(X3)
            so3.append(Y3)
        if plotCoexistenceOrWidth:
            if fromLattice:
                ax[1, 0].scatter(so2, T, color=color[j], marker="s")
            ax[1, 0].scatter(so3, T, color=color[j], marker="v")
            ax[1, 0].scatter(li3, T, color=color[j], marker="v", label=fr"$\alpha = ${RES[j]}")
            ax[1, 0].legend()
            ax[1, 0].set_xlabel(r"$\phi$")
            ax[1, 0].set_ylabel("$T_{thermo}/$")
        else:
            ax[1, 0].scatter(T, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
            ax[1, 0].legend()
            ax[1, 0].set_xlabel(r"$T$")
            ax[1, 0].set_ylabel(r"$\phi_s-\phi_l$")
           
   
    
    
if 1:
    A = np.load("final panel 3 new.npz")
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])

    E = nanify(A["E"])
    e = np.nanmean(E, axis = 2)


    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = False
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    
    
    color =  ["", "#02ecee","#fa03ee", "#00fe01","#fa8173","#ffd41c"]
    for j in range(len(RES)):
        li = []
        so = []
        so2 = []
        li3 = []
        so3 = []
        for i in range(len(T)):
            X3, Y3 = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], False)
            if fromLattice:
                so2.append(getSolidFromLattice(getY(L, XY[j, i], False)))
            li3.append(X3)
            so3.append(Y3)
        if plotCoexistenceOrWidth:
            if fromLattice:
                ax[2].scatter(so2, e[j]/T, color=color[j + 1], marker="s")
            ax[2].scatter(so3, e[j]/T, color=color[j + 1], marker="o", facecolor = "none", s = 100)
            ax[2].scatter(li3, e[j]/T, color=color[j + 1], marker="s", facecolor = "none", s = 100)
            ax[2].scatter([], [], c = color[j + 1], label=fr"{RES[j]}")
            
            ax[2].set_xlabel(r"$\phi$")
            ax[2].set_ylabel("$T/T_b$")
        else:
            ax[2].scatter(T/e[j], np.array(so3)-np.array(li3), label=fr"{RES[j]}")
            
            ax[2].set_xlabel(r"$T$")
            ax[2].set_ylabel(r"$\phi_s-\phi_l$")
            
        ax[2].legend(frameon = False, title = r"$\alpha=$", fontsize = fontsizeMain*0.8)
        ax[2].scatter([0.823], [1], marker = "*", color = "yellow", s = 200, edgecolors = "black")
        ax[2].scatter([0.87], [1], marker = "*", color = "yellow", s = 200, edgecolors = "black")
        ax[2].hlines(1, 0.80, 0.92, linestyle ="--", color = "gray", zorder = -20)
        ax[2].set_xlim(0.816, 0.90)
        
    




ax[0].text(0.05, 0.15, 'a)', transform=ax[0].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[1].text(0.05, 0.15, 'b)', transform=ax[1].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[2].text(0.05, 0.15, 'c)', transform=ax[2].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')

ax[0].text(0.4, 0.5, 'L+S', transform=ax[0].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[0].text(0.1, 0.5, 'L', transform=ax[0].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[0].text(0.85, 0.5, 'S', transform=ax[0].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')

ax[1].text(0.45, 0.45, 'L+S', transform=ax[1].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[1].text(0.1, 0.45, 'L', transform=ax[1].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
ax[1].text(0.89, 0.45, 'S', transform=ax[1].transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')

fig.savefig("EDMD_3panels.pdf")
plt.show()
