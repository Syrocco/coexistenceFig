import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

fontsizeMain = 26
fontsize = 18
plt.rcParams.update({
    'font.size' : fontsizeMain,                   # Set font size to 11pt
    'axes.labelsize': fontsizeMain,               # -> axis labels
    'legend.fontsize': fontsizeMain,              # -> legends
}) 
def Ncolors(N):
    cmap = plt.get_cmap('tab20')  
    colors = [cmap(i) for i in np.linspace(0, 1, N)]
    return colors

end = 2
start = 2
c1 = "#00fe01"
c2 = "#fa03ee"
c3 = "#fa8173"
import numpy as np
from functions2 import Data, dataArray
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from functions2 import energy as EE

def Ncolors(N):
    cmap = plt.get_cmap('tab20')  
    colors = [cmap(i) for i in np.linspace(0, 1, N)]
    return colors


def sortAccordingToA(A, *args):
    arg = np.argsort(A)
    # [:] allows to change by pointer the arrays (and not overwrite them), so no need to return them
    A[:] = A[arg]
    for arr in args:
        arr[:] = arr[arg]

def elementsWithin(sortedArr, leftValue, rightValue):
    return np.searchsorted(sortedArr, leftValue, 'left'), np.searchsorted(sortedArr, rightValue, 'right')


def areaOccupied(x, rad, leftValue, rightValue, Ly):
    def circleSegmentArea(r, h):
        f = r**2 - (r - h)**2
        g = np.arccos(1 - h/r)
        f[f < 0] = 0
        g[np.isnan(g)] = np.pi
        return r**2*g - (r - h)*np.sqrt(f)
    totalArea = 0
    
    i, j = elementsWithin(x, leftValue + rad, rightValue - rad)
    totalArea += (j - i)*np.pi*rad**2
    
    i, j = elementsWithin(x, leftValue - rad, leftValue + rad)
    h = (x[i:j] - leftValue) + rad
    area = circleSegmentArea(rad, h)
    totalArea += np.sum(area)
    
    i, j = elementsWithin(x, rightValue - rad, rightValue + rad)
    h = (rightValue - x[i:j]) + rad
    area = circleSegmentArea(rad, h)
    totalArea += np.sum(area)
    return totalArea/(Ly*(rightValue - leftValue))





def makeCellLists(x, y, Lx, Ly, size = 1):
    lenX = int(Lx/size)
    lenY = int(Ly/size)
    xcoord = (x/Lx*lenX).astype(int)
    ycoord = (y/Ly*lenY).astype(int)
    xcoord[xcoord == lenX] = lenX - 1
    ycoord[ycoord == lenY] = lenY - 1
    X = -np.ones((lenX, lenY)).astype(int)
    for i in range(len(x)):
        X[xcoord[i], ycoord[i]] = i
    return X

def pbc(dx, dL):
    if dx > dL/2:
        return dx - dL
    elif dx < -dL/2:
        return dx + dL
    return dx
    
def getSpacing(x, y, Lx, Ly):
    xx = []
    yy = []
    X = makeCellLists(x, y, Lx, Ly)
    hx = []
    hy = []
    lenX = np.shape(X)[0]
    lenY = np.shape(X)[1]
    for i in range(lenX):
        for j in range(lenY):
            if X[i, j] != -1:
                xp = x[X[i, j]]
                yp = y[X[i, j]]
                for k in range(-5, 5, 1):
                    for l in range(-5, 5, 1):
                            num = X[(i+k)%lenX, (j + l)%lenY]
                            if num != -1:
                                xn = x[num]
                                yn = y[num]
                               
                                dy = np.abs(pbc(yn - yp, Ly))
                                dx = np.abs(pbc(xn - xp, Lx))
                                dr = np.sqrt(dx*dx + dy*dy)
                            
       
                                if  0< dr < 3:
                                    if 0 <= dy < 0.3:
                                        hx.append(dx)
                                        if np.abs(xn - xp)**2 + np.abs(yn - yp) < 10:
                                            xx.append([xp, xn])
                                            yy.append([yp, yn])
                                    if 0 <= dx < 0.3:
                                        hy.append(dy)
                                        if np.abs(xn - xp)**2 + np.abs(yn - yp) < 10:
                                            xx.append([xp, xn])
                                            yy.append([yp, yn])
    return np.array(hx), np.array(hy), np.array(xx), np.array(yy)
def density(x, rad, Lx, Ly, N):
    
    
    sortAccordingToA(x, rad)

    

   
    phi = np.zeros(N - 2)
    L = np.linspace(0, Lx, N)
    L = 0.5*(L[2:] + L[:-2])
    slab = Lx/N
    
    radS = np.min(rad)
    radB = np.max(rad)
    xS = x[rad == radS]
    xB = x[rad == radB]
    
    for i in range(N - 2):
        leftValue = slab*(i + 1)
        rightValue = slab*(i + 2)
        phi[i] = areaOccupied(xS, radS, leftValue, rightValue, Ly) + areaOccupied(xB, radB, leftValue, rightValue, Ly)
    return L, phi
        



def shift(data, frame, N = 20):
    dump = data.dump
    if frame < 0:
        dump.jump_to_frame(dump.nframes + frame)
    else:
        dump.jump_to_frame(frame)

    x = dump.get_atompropf("x")
    rad = dump.get_atompropf("radius")

    Lx = data.Lx
    Ly = data.Ly

    
    L, phi = density(np.copy(x), rad, Lx, Ly, N)
        
    
    dphi = np.max(phi)-np.min(phi)
    where = phi > (np.max(phi) - dphi/3)
    if np.sum(where[:N//10]) and np.sum(where[-N//10:]):
        w = np.where(where == False)[0]
        try:
            return (x + Lx - L[w[-1] + 1])%Lx
        except:
            return  (x - L[where][0])%Lx
    return (x - L[where][0])%Lx

def ratioSpacing(data, frame = -1, N = 20):
    dump = data.dump
    if frame < 0:
        dump.jump_to_frame(dump.nframes + frame)
    else:
        dump.jump_to_frame(frame)

    x = dump.get_atompropf("x")
    y = dump.get_atompropf("y")
    rad = dump.get_atompropf("radius")

    Lx = data.Lx
    Ly = data.Ly

    
    L, phi = density(np.copy(x), rad, Lx, Ly, N)
    
    
    typ = dump.get_atompropf("type").astype(bool)
    x = dump.get_atompropf("x")[typ]
    y = dump.get_atompropf("y")[typ]
    rad = dump.get_atompropf("radius")[typ]
    
    dphi = np.max(phi)-np.min(phi)
    where = phi > (np.max(phi) - dphi/3)
    
    if np.sum(where[:N//10]) and np.sum(where[-N//10:]):
        w = np.where(where == False)[0]
        try:
            where2 = np.logical_or(x < L[w[0] - 1], x > L[w[-1] + 1])
        except:
            print("full liquid?")
            print(data.loc)
            xmin = L[where][0]
            xmax = L[where][-1]
            where2 = np.logical_and(x > xmin, x < xmax)
    else:

        xmin = L[where][0]
        xmax = L[where][-1]
        where2 = np.logical_and(x > xmin, x < xmax)
    xnew = x[where2]
    ynew = y[where2]
    xx, yy, XX, YY = getSpacing(xnew, ynew, Lx, Ly)

    x = dump.get_atompropf("x")
    y = dump.get_atompropf("y")
    rad = dump.get_atompropf("radius")
    color = ["r"  if r > 0.7 else "b" for r in rad]
    zorder = [-1  if r < 0.7 else 0 for r in rad]
    ax = plt.gca()
    for i in range(len(x)):
        
        circle = plt.Circle((x[i], y[i]), rad[i], color=color[i], alpha = 0.7, zorder = zorder[i]) 
        ax.add_artist(circle)
    
    plt.plot(XX.T, YY.T, color = "darkslategray", alpha = 0.7)

    return xx, yy, XX, YY

def getY(x, y, plot = False, incredible = False):
    if incredible:
        A = 0.01
        B = 100
    else:
        A = 1
        B = 1
    mx = np.min(x)
    Mx = np.max(x)
    x = x[~np.isnan(y)][start:-end]
    y = y[~np.isnan(y)][start:-end]
    coeff = np.polyfit(x, y, 1)
    a, b = coeff
    res = (1 - b)/a
    return x, y, a*x + b, res
        

def getSolidFromLattice(res):
    return np.pi*(1 + 0.476**2)/(res**2)

def getSolidAndLiquid(x, y, liquidReduced, solidReduced, mean = True):

    x = x[~np.isnan(y)][start:-end]
    y = y[~np.isnan(y)][start:-end]
    
    if y[0] < 1 and y[-1] > 1:
        argAbove = np.argmax(y > 1)
        argBelow = argAbove - 1
        liqB = liquidReduced[argBelow]
        solB = solidReduced[argBelow]
        liqA = liquidReduced[argAbove]
        solA = solidReduced[argAbove]
        dyA = y[argAbove] -1
        dyB = 1-y[argBelow] 
        if mean:
            return (dyA*liqA + dyB*liqB)/(dyB + dyA), (dyA*solA + dyB*solB)/(dyB + dyA)
        else:
            if dyA> dyB:
                return liqB, solB 
            else:
                return liqA, solA
                
    else:
        return None, None

def getSolidAndLiquid2(x, y, liquidReduced, solidReduced, plot = False, incredible = False):
    x = x[~np.isnan(y)][:-1]
    liquidReduced = liquidReduced[~np.isnan(y)][:-1]
    solidReduced = solidReduced[~np.isnan(y)][:-1]
    y = y[~np.isnan(y)][:-1]
    aL, bL = np.polyfit(x, liquidReduced, 1)
    aS, bS = np.polyfit(x, solidReduced, 1)
    return x, aS*x + bS, solidReduced, aL*x + bL, liquidReduced
   
 


    
       

def nanify(arr):
    arr[arr == 0] = np.nan
    return arr

def getDensityShift(data, start = 0, n = -1, N = 20):

    if n == -1:
        n = data.dump.nframes
    if start == 0:
        start = n//2
    elif start < 0:
        start = n + start
    
    dump = data.dump
    dens = 0
    count = 0
    rad = dump.get_atompropf("radius")
    Lx = data.Lx
    Ly = data.Ly
    for i in range(start, n):
        count += 1
       
        
        
        dump.jump_to_frame(i)
        x = shift(data, i, N)

        L, phi = density(np.copy(x), np.copy(rad),  Lx, Ly, N)
       
        dens += phi
    return L, dens/count

A = np.load("../4 panels/final panel 3 new.npz")


L = nanify(A["L"])
RES = nanify(A["RES"])
T = nanify(A["T"])
XY = nanify(A["ratioxy"])
YY = nanify(A["ratioyy"])
liquid = nanify(A["liquid"])
solid = nanify(A["solid"])
E = A["E"]
e = np.mean(E[0], axis = 1)

total = False
    


i = 4
j = -1
import matplotlib.gridspec as gridspec
    
if 1:
    fig = plt.figure(figsize=(10, 8),)
    gs = gridspec.GridSpec(3, 2,hspace=0.3,wspace=0.45)  
    
    
    ax0 = fig.add_subplot(gs[0, :])  # This takes the whole top row
    data = Data("N_1944dtnoise_0.500res_0.900gamma_0.100T_0.500phi_0.840000rat_0.476vo_3.500ao_2.100delta_0.500Lx_117.961Ly_37.800q_0.000v_1.dump")
    my_dpi = 16
    plt.axis("off")
    ratioSpacing(data, 4, N = 30)
    plt.xlim(0, data.Lx)
    plt.ylim(0, data.Ly)
    ax0.annotate('', xy=(0.5, 0), xytext=(0.5, 11.5), arrowprops=dict(arrowstyle='<-', lw=1.5))
    ax0.annotate('', xy=(0, 0.5), xytext=(10, 0.5), arrowprops=dict(arrowstyle='<-', lw=1.5))
    ax0.text(-5, 6, 'y', fontsize=fontsizeMain, va='center')  # Label for y arrow
    ax0.text(5, -5, 'x', fontsize=fontsizeMain, ha='center') 
    
    img = Image.open("drawing.png")

    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    inset_ax = ax0.inset_axes([-0.2, 0.17, 0.65, 0.65]) 
    
    inset_ax.imshow(img, interpolation='none')
    inset_ax.axis('off')  
    
    ax3 = fig.add_subplot(gs[1, :]) 
    if 0:
        data = Data("N_1944dtnoise_0.500res_0.970gamma_0.100T_0.020phi_0.845500rat_0.476vo_3.500ao_2.100delta_0.500Lx_117.194Ly_37.800q_0.000v_1.dump")
        l, phi = getDensityShift(data, 100, 280, 50)
        ax3.spines[['right', 'top']].set_visible(False)
        PHI = np.roll(phi, 14)*1.02
        LL = l/data.Lx
        dl = np.diff(L)[0]
        np.savez("LL and phi.npz", LL = LL, PHI = PHI)
    A = np.load("LL and phi.npz")
    PHI = A["PHI"]
    LL = A["LL"]
    ax3.plot(LL, PHI, color = "darkslategray")
    ax3.set_xlabel(r"$x/L_x$", labelpad = -10)
    ax3.set_ylabel(r"$\phi(x)$")
    ax3.spines[['right', 'top']].set_visible(False)
    phil = 0.838
    phis = 0.888
    ax3.hlines(phil, 0, 1, linestyle = "--", alpha = 0.4, color = "darkslategray")
    ax3.hlines(phis, 0, 1, linestyle = "--", alpha = 0.4, color = "darkslategray")
    ax3.set_xlim(0, 1)

    
    # Set the updated ticks and labels
    ax3.set_yticks([phil, 0.86, phis, 0.9])
    ax3.set_yticklabels([r"$\phi_l$", "0.86",r"$\phi_s$" , "0.90"])
    
    # Bottom-left subplot
    ax1 = fig.add_subplot(gs[2, 0])
    x, y, z, res = getY(L, XY[j, i], True, True)
    ax1.scatter(x, y, facecolor="none", edgecolor=c1, s=150)
    ax1.plot(x, z, c=c1)
    ax1.vlines(res, -1000, 1000, color="darkslategray", linestyle="--")
    ax1.set_xlim(2.075, 2.114)
    ax1.set_ylim(0.99, 1.0097)
    ax1.set_xlabel("$a_y$")
    ax1.set_ylabel("$a_y/a_x$")
    ax1.scatter([2.0933], [1], marker="*", color="yellow", s=500, edgecolors="black", zorder=100)
    
    
    # Bottom-right subplot
    ax2 = fig.add_subplot(gs[2, 1])
    x, Y, solidReduced, Z, liquidReduced = getSolidAndLiquid2(L, XY[j, i], liquid[j, i], solid[j, i], True, True)
    ax2.plot(x, Y, color=c3)
    ax2.scatter(x, solidReduced, color=c3, facecolor="none", s=150, label=r"$\phi_s$")
    ax2.plot(x, Z, color=c2)
    ax2.scatter(x, liquidReduced, facecolor="none", edgecolor=c2, s=150, label=r"$\phi_l$")
    ax2.vlines(res, -1000, 1000, color="darkslategray", linestyle="--")
    ax2.set_xlim(2.0685, 2.114)
    ax2.set_ylim(0.8296, 0.8888)
    ax2.legend(frameon=False)
    ax2.set_xlabel("$a_y$")
    ax2.set_ylabel("$\phi$", labelpad = -3)
    ax2.scatter([2.0933], [0.8786], marker="*", color="yellow", s=500, edgecolors="black", zorder=100)
    ax2.scatter([2.0933], [0.8338], marker="*", color="yellow", s=500, edgecolors="black", zorder=100)
    
    ax0.text(-0.06, 0.75, 'a)', transform=ax0.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
    ax1.text(0.85, 0.3, 'c)', transform=ax1.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
    ax2.text(0.85, 0.3, 'd)', transform=ax2.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
    ax3.text(0.92, 0.5, 'b)', transform=ax3.transAxes, fontsize=fontsizeMain*1.3, ha='left', va='bottom')
    
    plt.tight_layout()
    plt.savefig("finding good lattice.pdf", bbox_inches='tight', pad_inches=0)
