import matplotlib.pyplot as plt
import numpy as np

def Ncolors(N):
    cmap = plt.get_cmap('tab20')  
    colors = [cmap(i) for i in np.linspace(0, 1, N)]
    return colors

end = 2
start = 2
def getY(x, y, plot = False):
    mx = np.min(x)
    Mx = np.max(x)
    x = x[~np.isnan(y)][start:-end]
    y = y[~np.isnan(y)][start:-end]
    coeff = np.polyfit(x, y, 1)
    a, b = coeff
    res = (1 - b)/a
    if plot:
        plt.scatter(x, y)
        plt.plot(x, a*x + b)
        plt.vlines(res, np.min(y), np.max(y), color = "gray")
        plt.xlim(mx, Mx)
    return res

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

def getSolidAndLiquid2(x, y, liquidReduced, solidReduced, plot = False):
    mx = np.min(x)
    Mx = np.max(x)
    x = x[~np.isnan(y)][:-1]
    liquidReduced = liquidReduced[~np.isnan(y)][:-1]
    solidReduced = solidReduced[~np.isnan(y)][:-1]
    y = y[~np.isnan(y)][:-1]
    aL, bL = np.polyfit(x, liquidReduced, 1)
    aS, bS = np.polyfit(x, solidReduced, 1)
    X = getY(x, y)
    if plot:
        plt.plot(x,  aL*x + bL, color = "C0")
        plt.scatter(x, liquidReduced, color = "C0")
        plt.plot(x,  aS*x + bS, color = "C1")
        plt.scatter(x, solidReduced, color = "C1")
        plt.vlines(X, np.min(liquidReduced), np.max(solidReduced), color = "gray")
        plt.xlim(mx, Mx)
    return aL*X + bL, aS*X + bS
    
def nanify(arr):
    arr[arr == 0] = np.nan
    return arr 

def postProcess(XY, YY, liquid, solid):
    
    XY[0, :, :10] = np.nan
    YY[0, :, :10] = np.nan
    solid[0, :, :10] = np.nan
    liquid[0, :, :10] = np.nan
    
    XY[0, :, -3:] = np.nan
    YY[0, :, -3:] = np.nan
    solid[0, :, -3:] = np.nan
    liquid[0, :, -3:] = np.nan
