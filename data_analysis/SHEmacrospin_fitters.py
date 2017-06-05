# -*- coding: utf-8 -*-
"""
Created on Thu Oct 08 18:42:19 2015

@author: admin
"""
import spinmob as sm
import SHEmacrospin_calculators as SHEc
import numpy as np
import scipy as sp
import scipy.optimize as spop
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Magnet parameters
M_s  = 0.12
dims = [4.0, 0.5, 0.05]
vol  = dims[0]*dims[1]*dims[2]
demag = SHEc.demagFactors(dims)
Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]   

# Various constants
pi = np.pi
hbar = 1.054572e-34
e_ch = 1.602e-19
gamma = 1.76e11
kappa = hbar*gamma/(2*vol*e_ch*M_s**2)

   
def fit_spectrum(f, V, sigma = None, absolute_sigma = True, startingpt = None, plot = True, subplot = False, gcf = False, ylabel = r'$\delta R$' + r'$(\Omega)$'):
    """
    This fits a spectrum to a single asymmetric lorentzian + symmetric lorentzian
    """
    def fitfunc(f, a1, a2, f0, G, off):
        V_th = a1 * 1 / ((f**2-f0**2)**2 + (2*pi*G * f)**2) + a2 * (f**2-f0**2) / ((f**2-f0**2)**2 + (2*pi*G * f)**2) + off
        return V_th
        
    # Figure out a reasonable starting point for the fitting
    if startingpt == None:
        maxV = np.argmax(V); minV = np.argmin(V);
        start_off = (V[0]+V[-1])/2     
        start_f0 = (f[maxV] + f[minV])/2
        start_G = np.abs(f[maxV] - f[minV])/2
        start_a1  = max(V)*np.sign(max(V)+min(V)-2*start_off) * (2*pi*start_G * start_f0)**2
        start_a2  = max(V)*np.sign(maxV-minV) * (2*pi*start_G)**2
    
        startingpt = [start_a1, start_a2, start_f0, start_G, start_off]
        print startingpt
    
    try:       
        popt, pcov = spop.curve_fit(fitfunc, f, V, p0 = startingpt, maxfev = 10000)
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        popt = np.zeros(5)
        pcov = np.zeros([5,5])
        perr = np.zeros(5)
        
    # Take the absolute value of the damping constant
    popt[3] = np.abs(popt[3])   
   
    # Plot the result
    if plot:
        if not gcf: 
            plt.figure()
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
            plt.subplot(gs[0])        
        plt.errorbar(f,V, yerr = sigma, fmt='o')
        if np.abs(popt[2]) > 1e8  and np.abs(popt[2]) < 20e9:
            fit_x   = (np.linspace(0, f[-1], 200))
            fit_y   = fitfunc(fit_x, popt[0], popt[1], popt[2], popt[3], popt[4])
            fit_sy  = fitfunc(fit_x, popt[0],       0, popt[2], popt[3], popt[4])
            fit_asy = fitfunc(fit_x,       0, popt[1], popt[2], popt[3], popt[4])
            plt.plot(fit_x, fit_y)
            if subplot:
                plt.plot(fit_x, fit_sy, fit_x, fit_asy)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel(ylabel)
        plt.figtext(0,0.02,r'$f_0 =$' + str(popt[2]/1e9) + "GHz")
        plt.tight_layout()
        
        # Residuals
        if not gcf and np.abs(popt[2]) > 1e8 and np.abs(popt[2]) < 20e9:      
            plt.subplot(gs[1])
            plt.plot(f, fitfunc(f, popt[0], popt[1], popt[2], popt[3], popt[4]) - V)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel(ylabel)
            plt.tight_layout() 
    
    return popt, perr
   
 
def fit_spectrum_kink(f, V, startingpt = None, plot = True, subplot = False, gcf = False):
    """
    This fits a kink spectrum to two asymmetric lorentzians + symmetric lorentzians
    """
    def fitfunc(f, a0, b0, f0, G0, a1, b1, f1, G1, off):
        V_th = a0 * 1 / ((f**2-f0**2)**2 + (2*pi*G0 * f)**2) + b0 * (f**2-f0**2) / ((f**2-f0**2)**2 + (2*pi*G0 * f)**2) + \
        a1 * 1 / ((f**2-f1**2)**2 + (2*pi*G1 * f)**2) + b1 * (f**2-f1**2) / ((f**2-f1**2)**2 + (2*pi*G1 * f)**2) + \
        off
        return V_th
        
    # Figure out a reasonable starting point for the fitting
    if startingpt == None:
        maxV = np.argmax(V); minV = np.argmin(V);
        print([maxV, minV])
        start_off = (V[0]+V[-1])/2     
        start_f0 = f[maxV]
        start_G0 = np.abs(f[maxV] - f[minV])/2
        start_a0  = max(V)*np.sign(max(V)+min(V)-2*start_off) * (2*pi*start_G0 * start_f0)**2
        start_b0  = max(V)*np.sign(maxV-minV) * (2*pi*start_G0)**2
        start_a1  = max(V)*np.sign(max(V)+min(V)-2*start_off) * (2*pi*start_G0 * start_f0)**2
        start_b1  = max(V)*np.sign(maxV-minV) * (2*pi*start_G0)**2       
        start_f1 = f[minV]
        start_G1 = start_G0
    
        startingpt = [start_a0, start_b0, start_f0, start_G0, start_a1, start_b1, start_f1, start_G1, start_off]
        print startingpt
        
    popt, pcov = spop.curve_fit(fitfunc, f, V, p0 = startingpt, maxfev = 10000)
    perr = np.sqrt(np.diag(pcov))
 
    # Take the absolute value of the frequency and damping constant
    popt[2] = np.abs(popt[2])  
    popt[3] = np.abs(popt[3])
    popt[6] = np.abs(popt[6])  
    popt[7] = np.abs(popt[7])   
    if popt[2] > popt[6]:
        temp = popt[0]; popt[0] = popt[4]; popt[4] = temp   
        temp = popt[1]; popt[1] = popt[5]; popt[5] = temp
        temp = popt[2]; popt[2] = popt[6]; popt[6] = temp
        temp = popt[3]; popt[3] = popt[7]; popt[7] = temp
        
   
     # Plot the result
    if plot:
        fit_x   = (np.linspace(f[1], f[-1], 200))
        fit_y    = fitfunc(fit_x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8])
        fit_sy0  = fitfunc(fit_x, popt[0],       0, popt[2], popt[3],       0,       0, popt[6], popt[7], popt[8])
        fit_asy0 = fitfunc(fit_x,       0, popt[1], popt[2], popt[3],       0,       0, popt[6], popt[7], popt[8])
        fit_sy1  = fitfunc(fit_x,       0,       0, popt[2], popt[3], popt[4],       0, popt[6], popt[7], popt[8])
        fit_asy1 = fitfunc(fit_x,       0,       0, popt[2], popt[3],       0, popt[5], popt[6], popt[7], popt[8])        
        if not gcf: 
            plt.figure()
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
            plt.subplot(gs[0])
        plt.plot(f, V, 'o', fit_x, fit_y)
        if subplot:
            plt.plot(fit_x, fit_sy0, fit_x, fit_asy0, fit_x, fit_sy1, fit_x, fit_asy1)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel(r'$\delta R$' + r'$(\Omega)$')
        plt.tight_layout()
        
        # Residuals
        if not gcf:
            plt.subplot(gs[1])
            plt.plot(f, fitfunc(f, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]) - V)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel(r'$\delta R$' + r'$(\Omega)$')
            plt.tight_layout()
    
    return popt, perr

  
def fit_spectrum_multiple(startingpt = None, xval = None, yerr = True, plot = True, subplot = False, plotres = True):
    """
    This fits a bunch of spectra located in the same folder.
    """    
    
    # Import the files (user input)
    a = sm.data.load_multiple()   
    
    # Fit 'em
    popt = []; perr = [];
    for k in range(0,len(a)):        
        if not yerr: 
            err = None
        else:
            err = a[k][5] 
        poptt, perrt = fit_spectrum(a[k][0], a[k][3], sigma = err, absolute_sigma = True, startingpt = startingpt, plot = plot, subplot = subplot)
        popt.append(poptt); perr.append(perrt)
    popt = np.array(popt); perr = np.array(perr)

    # Plot the results
    if plotres:
        if xval == None:
            xval = range(0,len(a))
        plt.figure(figsize = (20,5), dpi = 80); ax = []
        
        ax.append(plt.subplot(1,2,1)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,2], yerr = perr[:,2], color = 'b')
        ax[0].set_ylabel('Resonance Frequency (Hz)', color='b')
        for tl in ax[0].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[0].twinx())   
        plt.errorbar(x = xval, y = popt[:,3], yerr = perr[:,3], color = 'r')
        ax[1].set_ylabel('Resonance Linewidth (Hz)', color='r')
        for tl in ax[1].get_yticklabels(): tl.set_color('r');
        
        ax.append(plt.subplot(1,2,2)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,0], yerr = perr[:,0], color = 'b')
        ax[2].set_ylabel('Symmetric Peak Amplitude', color='b')
        for tl in ax[2].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[2].twinx())  
        plt.errorbar(x = xval, y = popt[:,1], yerr = perr[:,1], color = 'r')
        ax[3].set_ylabel('Antisymmetric Peak Amplitude', color='r')
        for tl in ax[3].get_yticklabels(): tl.set_color('r')

    return popt, perr


def fit_spectrum_kink_multiple(startingpt = None, xval = None, plot = True, subplot = False, plotres = True):
    """
    This fits a bunch of spectra located in the same folder.
    """    
    
    # Import the files (user input)
    a = sm.data.load_multiple()   
    
    # Fit 'em
    popt = []; perr = [];
    for k in range(0,len(a)):
        poptt, perrt = fit_spectrum_kink(a[k][0], a[k][3], startingpt = startingpt, plot = plot, subplot = subplot)
        popt.append(poptt); perr.append(perrt)
    popt = np.array(popt); perr = np.array(perr)

    # Plot the results
    if plotres:
        if xval == None:
            xval = range(0,len(a))
        plt.figure(figsize = (20,10), dpi = 80); ax = []
        
        ax.append(plt.subplot(2,2,1)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,2], yerr = perr[:,2], color = 'b')
        ax[0].set_ylabel('Resonance Frequency (Hz)', color='b')
        for tl in ax[0].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[0].twinx())   
        plt.errorbar(x = xval, y = popt[:,3], yerr = perr[:,3], color = 'r')
        ax[1].set_ylabel('Resonance Linewidth (Hz)', color='r')
        for tl in ax[1].get_yticklabels(): tl.set_color('r');
        
        ax.append(plt.subplot(2,2,2)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,0], yerr = perr[:,0], color = 'b')
        ax[2].set_ylabel('Symmetric Peak Amplitude', color='b')
        for tl in ax[2].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[2].twinx())  
        plt.errorbar(x = xval, y = popt[:,1], yerr = perr[:,1], color = 'r')
        ax[3].set_ylabel('Antisymmetric Peak Amplitude', color='r')
        for tl in ax[3].get_yticklabels(): tl.set_color('r')
            
        ax.append(plt.subplot(2,2,3)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,6], yerr = perr[:,6], color = 'b')
        ax[4].set_ylabel('Resonance Frequency (Hz)', color='b')
        for tl in ax[4].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[4].twinx())   
        plt.errorbar(x = xval, y = popt[:,7], yerr = perr[:,7], color = 'r')
        ax[5].set_ylabel('Resonance Linewidth (Hz)', color='r')
        for tl in ax[5].get_yticklabels(): tl.set_color('r');
        
        ax.append(plt.subplot(2,2,4)); plt.hold(True)
        plt.errorbar(x = xval, y = popt[:,4], yerr = perr[:,4], color = 'b')
        ax[6].set_ylabel('Symmetric Peak Amplitude', color='b')
        for tl in ax[6].get_yticklabels(): tl.set_color('b');
            
        ax.append(ax[6].twinx())  
        plt.errorbar(x = xval, y = popt[:,5], yerr = perr[:,5], color = 'r')
        ax[7].set_ylabel('Antisymmetric Peak Amplitude', color='r')
        for tl in ax[7].get_yticklabels(): tl.set_color('r')
            
        plt.tight_layout()

    return popt, perr


def fit_f0_vs_demag(H, th_H, f0, thoff = None, demag = None, startingpt = [1.0, 1.0, 0.0], eqstart = None):
    """
    This fits the measured resonance frequencies to the magnitude of the applied
    external field and the demagnetizing factors.
    """    
    # Convert to radians if it's obvious that th_H is in degrees
    if np.any(th_H)>2*pi: th_H = np.array(th_H)*pi/180    
    
    # Definition of the fitting function
    def fitfunc(x, N1, N2, th0):
        f0_th = []  
        print x
        for k in x:   
            print H, N1, N2, th0
            
            # Get starting value for equilibrium angle calculation
            if eqstart == 'th_H': startpt = (M_s*np.cos(k), M_s*np.sin(k))
            elif isinstance(eqstart,float): startpt = (M_s*np.cos(eqstart), M_s*np.sin(eqstart))
            else: startpt = None
            th = SHEc.eqAngle2(M_s, H_mag = k, th_H = th_H, demag = [N1, N2, 4*pi-N1-N2])    
            
            # Calculate resonance frequency
            NzxP = 1 - 2*N1 - N2
            NzyP = 1 - N1 - 2*N2  
            NyxP = N2 - N1
            f0_th.append( gamma * np.real( \
                                  np.sqrt(  ( ( k*np.cos(th_H) + NzxP*M_s*np.cos(th) ) * np.cos(th) + ( k*np.sin(th_H) + NzyP*M_s*np.sin(th) ) * np.sin(th)  ) * \
                                         (    ( k*np.cos(th_H) + NyxP*M_s*np.cos(th) ) * np.cos(th) + ( k*np.sin(th_H) - NyxP*M_s*np.sin(th) ) * np.sin(th)  ) + 0j) / (2*pi) ) )   
                
             
        return f0_th
   
    # Initially fit the data with the supplied offset angle
    if thoff == None: thfit = 0
    else: thfit = thoff 
    popt, pcov = spop.curve_fit(lambda x, N1, N2: fitfunc(x,N1,N2,thfit), th_H, f0, p0 = startingpt[0:2:1], maxfev = 10000)
    popt = np.append(popt,thoff)   
    
    # Now also fit over the offset angle using the previous fit as starting values
    if thoff == None:
        if len(startingpt)<4: startingpt.append(0)
        popt, pcov = spop.curve_fit(fitfunc, th_H, f0, p0 = [popt[0], popt[1], startingpt[3]], maxfev = 10000)

    # Calculate errors        
    perr = np.sqrt(np.diag(pcov))

    # Plot the result
    fit_x = (np.linspace(min(H), max(H), 200))
    fit_y = fitfunc(fit_x, popt[0], popt[1], popt[2])
    plt.figure()
    plt.plot(np.array(H), f0, 'o', fit_x, fit_y)
    plt.xlabel('Magnetic Field (mT)')
    plt.ylabel('Resonance Frequency (Hz)')
    #plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Nx = ' + str(np.round(popt[1],2)) + ', Ny = ' + str(np.round(popt[2],2)) + ', Nz = ' + \
    #          str(np.round(4*pi-popt[2]-popt[1],2)))        
    plt.title('Nx = ' + str(np.round(popt[0],2)) + ', Ny = ' + str(np.round(popt[1],2)) + ', Nz = ' + \
              str(np.round(4*pi-popt[1]-popt[0],2)) + r', $\theta_{off}$ = ' + str(np.round(180/pi*popt[2],2)))
    plt.tight_layout()
    
    return popt, perr
      
    
def fit_f0_vs_th_Hdemag(th_H, f0, thoff = None, demag = None, startingpt = [0.0, 1.0, 1.0, 0.0], eqstart = None):
    """
    This fits the measured resonance frequencies to the magnitude of the applied
    external field and the demagnetizing factors.
    """    
    # Convert to radians if it's obvious that th_H is in degrees
    if max(th_H)>2*pi: th_H = np.array(th_H)*pi/180    
    
    # Definition of the fitting function
    def fitfunc(x, H, N1, N2, th0):
        f0_th = []     
        for k in x:              
            # Get starting value for equilibrium angle calculation
            if eqstart == 'th_H': startpt = (M_s*np.cos(k), M_s*np.sin(k))
            elif isinstance(eqstart,float): startpt = (M_s*np.cos(eqstart), M_s*np.sin(eqstart))
            else: startpt = None
            th = SHEc.eqAngle2(M_s, H_mag = H, th_H = k, demag = [N1, N2, 4*pi-N1-N2])    
            
            # Calculate resonance frequency
            NzxP = 1 - 2*N1 - N2
            NzyP = 1 - N1 - 2*N2  
            NyxP = N2 - N1
            kp   = k - th0
            f0_th.append( gamma * np.real( \
                                  np.sqrt(  ( ( H*np.cos(kp) + NzxP*M_s*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) + NzyP*M_s*np.sin(th) ) * np.sin(th)  ) * \
                                         (    ( H*np.cos(kp) + NyxP*M_s*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) - NyxP*M_s*np.sin(th) ) * np.sin(th)  ) + 0j) / (2*pi) ) )   
             
        print [H, N1, N2]
        return f0_th
   
    # Initially fit the data with the supplied offset angle
    if thoff == None: thfit = 0
    else: thfit = thoff 
    popt, pcov = spop.curve_fit(lambda x, H, N1, N2: fitfunc(x,H,N1,N2,thfit), th_H, f0, p0 = startingpt[0:3:1], maxfev = 10000)
    popt = np.append(popt,thoff)   
    
    # Now also fit over the offset angle using the previous fit as starting values
    if thoff == None:
        print "thoff..."
        if len(startingpt)<4: startingpt.append(0)
        popt, pcov = spop.curve_fit(fitfunc, th_H, f0, p0 = [popt[0], popt[1], popt[2], startingpt[3]], maxfev = 10000)

    # Calculate errors        
    perr = np.sqrt(np.diag(pcov))

    # Plot the result
    fit_x = (np.linspace(min(th_H), max(th_H), 200))
    fit_y = fitfunc(fit_x, popt[0], popt[1], popt[2], popt[3])
    plt.figure()
    plt.plot(np.array(th_H)*180/pi, f0, 'o', fit_x*180/pi, fit_y)
    plt.xlabel('Magnetic Field Angle')
    plt.ylabel('Resonance Frequency (Hz)')
    #plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Nx = ' + str(np.round(popt[1],2)) + ', Ny = ' + str(np.round(popt[2],2)) + ', Nz = ' + \
    #          str(np.round(4*pi-popt[2]-popt[1],2)))        
    plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Nx = ' + str(np.round(popt[1],2)) + ', Ny = ' + str(np.round(popt[2],2)) + ', Nz = ' + \
              str(np.round(4*pi-popt[2]-popt[1],2)) + r', $\theta_{off}$ = ' + str(np.round(180/pi*popt[3],2)))
    plt.tight_layout()
    
    return popt, perr
    
    
def fit_f0_vs_th_Msdemag(th_H, f0, sf0 = None, thoff = None, H = None, demag = None, startingpt = [0.0, 1.0, 1.0, 0.0], eqstart = None):
    """
    This fits the measured resonance frequencies to the magnitude of the applied
    external field and the demagnetizing factors.
    """    
    # Convert to radians if it's obvious that th_H is in degrees
    if max(th_H)>2*pi: th_H = np.array(th_H)*pi/180    
    
    # Definition of the fitting function
    def fitfunc(x, M_s, N1, N2, th0):
        f0_th = []     
        for k in x:        
            print M_s, N1, N2, th0
            
            # Get starting value for equilibrium angle calculation
            if eqstart == 'th_H': startpt = (M_s*np.cos(k), M_s*np.sin(k))
            elif isinstance(eqstart,float): startpt = (M_s*np.cos(eqstart), M_s*np.sin(eqstart))
            else: startpt = None
            th = SHEc.eqAngle2(M_s, H_mag = H, th_H = k, demag = [N1, N2, 4*pi-N1-N2])    
            
            # Calculate resonance frequency
            NzxP = 1 - 2*N1 - N2
            NzyP = 1 - N1 - 2*N2  
            NyxP = N2 - N1
            kp   = k - th0
            f0_th.append( gamma * np.real( \
                                  np.sqrt(  ( ( H*np.cos(kp) + NzxP*M_s*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) + NzyP*M_s*np.sin(th) ) * np.sin(th)  ) * \
                                         (    ( H*np.cos(kp) + NyxP*M_s*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) - NyxP*M_s*np.sin(th) ) * np.sin(th)  ) + 0j) / (2*pi) ) )   
            
             
        return f0_th
   
    # Initially fit the data with the supplied offset angle
    if thoff == None: thfit = 0
    else: thfit = thoff 
    popt, pcov = spop.curve_fit(lambda x, H, N1, N2: fitfunc(x,H,N1,N2,thfit), th_H, f0, sigma = sf0, p0 = startingpt[0:3:1], maxfev = 10000)
    popt = np.append(popt,thoff)   
    
    # Now also fit over the offset angle using the previous fit as starting values
    if thoff == None:
        if len(startingpt)<4: startingpt.append(0)
        popt, pcov = spop.curve_fit(fitfunc, th_H, f0, p0 = [popt[0], popt[1], popt[2], startingpt[3]], maxfev = 10000)

    # Calculate errors        
    perr = np.sqrt(np.diag(pcov))

    # Plot the result
    fit_x = (np.linspace(min(th_H), max(th_H), 200))
    fit_y = fitfunc(fit_x, popt[0], popt[1], popt[2], popt[3])
    plt.figure()
    plt.plot(np.array(th_H)*180/pi, f0, 'o', fit_x*180/pi, fit_y)
    plt.xlabel('Magnetic Field Angle')
    plt.ylabel('Resonance Frequency (Hz)')
    #plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Nx = ' + str(np.round(popt[1],2)) + ', Ny = ' + str(np.round(popt[2],2)) + ', Nz = ' + \
    #          str(np.round(4*pi-popt[2]-popt[1],2)))        
    plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Nx = ' + str(np.round(popt[1],2)) + ', Ny = ' + str(np.round(popt[2],2)) + ', Nz = ' + \
              str(np.round(4*pi-popt[2]-popt[1],2)) + r', $\theta_{off}$ = ' + str(np.round(180/pi*popt[3],2)))
    plt.tight_layout()
    
    return popt, perr

    
def fit_f0_vs_th_H_HMs(th_H, f0, thoff = None, demag = None, startingpt = [0.1, 0.08, 0], eqstart = None):
    """
    This fits the measured resonance frequencies to the magnitude of the applied
    external field and the demagnetizing factors.
    """    
    # Convert to radians if it's obvious that th_H is in degrees
    if max(th_H)>2*pi: th_H = np.array(th_H)*pi/180    
    
    # Definition of the fitting function
    def fitfunc(x, H, Ms, th0):
        global demag
        f0_th = []     
        for k in x:
            # Get starting value for equilibrium angle calculation
            if eqstart == 'th_H': startpt = (M_s*np.cos(k), M_s*np.sin(k))
            elif isinstance(eqstart,float): startpt = (M_s*np.cos(eqstart), M_s*np.sin(eqstart))
            else: startpt = None
            th = SHEc.equilibriumAngle(M_s = Ms, H_mag = H, th_H = k, demag = demag, startingpt = startpt)       
            
            # Calculate resonance frequency
            NzxP = 1 - 2*demag[1] - demag[2]
            NzyP = 1 - demag[1] - 2*demag[2]
            NyxP = demag[2] - demag[1]
            kp   = k - th0
            f0_th.append( gamma * np.real( \
                                  np.sqrt(  ( ( H*np.cos(kp) + NzxP*Ms*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) + NzyP*Ms*np.sin(th) ) * np.sin(th)  ) * \
                                         (    ( H*np.cos(kp) + NyxP*Ms*np.cos(th) ) * np.cos(th) + ( H*np.sin(kp) - NyxP*Ms*np.sin(th) ) * np.sin(th)  ) + 0j) / (2*pi) ) )   
                                    
        return f0_th

    # Initially fit the data with the supplied offset angle
    if thoff == None: thfit = 0
    else: thfit = thoff 
    popt, pcov = spop.curve_fit(lambda x, H, Ms: fitfunc(x,H,Ms,thfit), th_H, f0, p0 = startingpt[0:2:1], maxfev = 10000)
    popt = np.append(popt,thoff)   
    
    # Now also fit over the offset angle using the previous fit as starting values
    if thoff == None:
        if len(startingpt)<4: startingpt.append(0)
        popt, pcov = spop.curve_fit(fitfunc, th_H, f0, p0 = [popt[0], popt[1], startingpt[2]], maxfev = 10000)

    # Calculate errors        
    perr = np.sqrt(np.diag(pcov))

    # Plot the result
    fit_x = (np.linspace(min(th_H), max(th_H), 200))
    fit_y = fitfunc(fit_x, popt[0], popt[1], popt[2])
    plt.figure()
    plt.plot(np.array(th_H)*180/pi, f0, 'o', fit_x*180/pi, fit_y)
    plt.plot(fit_x*180/pi, fitfunc(fit_x, 0.03, 0.03, 0))
    plt.xlabel('Magnetic Field Angle')
    plt.ylabel('Resonance Frequency (Hz)')       
    plt.title('H = ' + str(np.round(popt[0]*1000,2)) + 'mT, Ms = ' + str(np.round(popt[1],2)) + r', $\theta_{off}$ = ' + str(np.round(180/pi*popt[2],2)))
    plt.tight_layout() 
    
    return popt, perr
    
    
def plot_spectrum_multiple(dbs = None, dbcol = 3, dberr = 5):
    # Import the files (user input)
    if dbs == None: dbs = sm.data.load_multiple();
    
    offset = 0
    plt.figure(); plt.hold(True)
    for k in range(len(dbs)): 
        offset += -np.min(dbs[k][dbcol])
        plt.errorbar(dbs[k][0], dbs[k][dbcol] + offset, dbs[k][dberr], fmt='-o')
        offset += np.max(dbs[k][dbcol])
    plt.xlabel("Frequency (Hz)")
    plt.ylabel(r'$\delta R \; (\Omega$,' + " relative)")
    
    
def plot_fit(params, gcf = True, subplot = False):
    """
    Plots a fit of the anti-sym + sym lorentzians according to the provided parameters.
    """
    