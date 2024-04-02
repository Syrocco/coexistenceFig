import numpy as np
import matplotlib.pyplot as plt
from helper import nanify, getY, getSolidAndLiquid,  getSolidFromLattice, getSolidAndLiquid2, Ncolors, postProcess



fig, ax = plt.subplots(2, 2, figsize=(14, 8))


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
    fromLattice = True
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
            ax[0, 0].scatter(so2, RES[8:], color=color[0], marker="s")
        ax[0, 0].scatter(so3, RES[8:], color=color[0], marker="v")
        ax[0, 0].scatter(li3, RES[8:], color=color[0], marker="v", label=fr"$T = ${T[i]}")
        ax[0, 0].legend()
        ax[0, 0].set_xlabel(r"$\phi$")
        ax[0, 0].set_ylabel(r"$\alpha$")
    else:
        ax[0, 0].scatter(T, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
        ax[0, 0].legend()
        ax[0, 0].set_xlabel(r"$T$")
        ax[0, 0].set_ylabel(r"$\phi_s-\phi_l$")

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
    fromLattice = True
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    color = Ncolors(7)
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
                ax[0, 1].scatter(so2, T/e, color = color[j], marker = "s")
            ax[0, 1].scatter(so3, T/e, color = color[j], marker = "v")
            ax[0, 1].scatter(li3, T/e, color = color[j], marker = "v" , label = fr"$\alpha = ${RES[j]}")
            ax[0, 1].legend()
            ax[0, 1].set_xlabel(r"$\phi$")
            ax[0, 1].set_ylabel("$T_{thermo}/E$")
            ax[0, 1].vlines(0.823, 1, 6)
            ax[0, 1].vlines(0.87, 1, 6)
        else:
            ### Plot the the size of the region with respect to T instead of the coexistence vs T.
            ax[0, 1].scatter(T/e, np.array(so3)-np.array(li3) , label = fr"$\alpha = ${RES[j]}")
            ax[0, 1].legend()
            ax[0, 1].set_xlabel(r"$T$")
            ax[0, 1].set_ylabel(r"$\phi_s-\phi_l$")


if 1:
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
            ax[1, 0].set_ylabel("$T_{thermo}/??$")
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
            ax[1, 0].set_ylabel("$T_{thermo}/??$")
        else:
            ax[1, 0].scatter(T, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
            ax[1, 0].legend()
            ax[1, 0].set_xlabel(r"$T$")
            ax[1, 0].set_ylabel(r"$\phi_s-\phi_l$")
           
   
    
    
if 1:
    A = np.load("arraysCleaner.npz")
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])
    E = np.load("arraysCleaner E.npz")

    E = E["E"]
    e = np.mean(E, axis = 2)


    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = True
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    
    
    color = Ncolors(7)
    for j in range(len(RES) - 1):
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
                ax[1, 1].scatter(so2, T/e[j], color=color[j + 1], marker="s")
            ax[1, 1].scatter(so3, T/e[j], color=color[j + 1], marker="v")
            ax[1, 1].scatter(li3, T/e[j], color=color[j + 1], marker="v", label=fr"$\alpha = ${RES[j]}")
            ax[1, 1].legend()
            ax[1, 1].set_xlabel(r"$\phi$")
            ax[1, 1].set_ylabel("$T_{thermo}/E$")
        else:
            ax[1, 1].scatter(T/e[j], np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
            ax[1, 1].legend()
            ax[1, 1].set_xlabel(r"$T$")
            ax[1, 1].set_ylabel(r"$\phi_s-\phi_l$")
    
    A = np.load("Equilibrium.npz")
    L = nanify(A["L"])
    RES = nanify(A["RES"])
    T = nanify(A["T"])
    XY = nanify(A["ratioxy"])
    YY = nanify(A["ratioyy"])
    liquid = nanify(A["liquid"])
    solid = nanify(A["solid"])
    e = np.copy(T)
    


    ### if True plot in square also the solid coexistence predicted from the lattice spacing value
    fromLattice = True
    ### If True plot the coexistence, if False plot the width of the coexistence zone instead
    plotCoexistenceOrWidth = True
    
    
    color = Ncolors(7)
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
                ax[1, 1].scatter(so2, T/e, color=color[j], marker="s")
            ax[1, 1].scatter(so3, T/e, color=color[j], marker="v")
            ax[1, 1].scatter(li3, T/e, color=color[j], marker="v", label=fr"$\alpha = ${RES[j]}")
            ax[1, 1].legend()
            ax[1, 1].set_xlabel(r"$\phi$")
            ax[1, 1].set_ylabel("$T_{thermo}/E$")
        else:
            ax[1, 1].scatter(T/e, np.array(so3)-np.array(li3), label=fr"$\alpha = ${RES[j]}")
            ax[1, 1].legend()
            ax[1, 1].set_xlabel(r"$T$")
            ax[1, 1].set_ylabel(r"$\phi_s-\phi_l$")


plt.tight_layout()

plt.show()
