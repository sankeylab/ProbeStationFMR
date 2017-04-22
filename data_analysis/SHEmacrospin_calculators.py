# -*- coding: utf-8 -*-
"""
Created on Thu Oct 08 15:36:26 2015

@author: admin
"""

import numpy as np
import scipy as sp
import scipy.optimize as spop
import scipy.special as spsp 
import matplotlib.pyplot as plt

# Various constants
pi = np.pi
hbar = 1.054572e-34
e_ch = 1.602e-19
gamma = 1.76e11
mu_0 = 1.2566370e-6

#vol  = 4*pi/3 * dims(1) * dims(2) * dims(3)
#kappa = vol**-1 * hbar*gamma / (2*e_ch)             # Spin Hall prefactor


def demagFactors(dims):
    """
    Calculates the demagnetization factors for the elongated ellipsoid
    """
    a = float(dims[0])
    b = float(dims[1])
    c = float(dims[2])
    
    e = np.sqrt(1-b**2/a**2)
    f = np.sqrt(1-e**2)
    K = spsp.ellipk(e)
    E = spsp.ellipe(e)
    
    case = np.argmax([a/b, b/(2*c)])
    demag = []
    
    if case == 0:
        demag.append( ((b*c/a**2) * (np.log(4*a/(b+c)) - 1)) )
        demag.append( (c/(b+c) - 1/2*(b*c/a**2)*np.log(4*a/(b+c)) + b*c*(3*b+c)/(4*a**2*(b+c))) )
        demag.append( (b/(b+c) - 1/2*(b*c/a**2)*np.log(4*a/(b+c)) + b*c*(b+3*c)/(4*a**2*(b+c))) )
        
    elif case == 1:
        demag.append( (c/a * f * (K - E)/e**2) )
        demag.append( (c/a * (E - K*f**2)/(e**2*f) ) )
        demag.append( (1 - c*E / (a*f) ) )
        
    return np.array(demag)
 

def equilibriumAngleEnergy(M_s = 8e5, alpha = 0.0075, dims = [1.0, 1.0, 1.0], H_mag=0, th_H=0, eta = 1, I_dc = 0, demag = None, fig = False):
    """
    Calculates the equilibrium angle by using the more reliable energy formalism
    """
    if demag == None: demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]
    
    bnds = ((0, pi/2),)
    
    Energy = lambda x: -M_s * H_mag * np.cos(x-th_H) + mu_0/2 * demag[1] * M_s**2 * np.sin(x)**2
    th_M = sp.optimize.minimize(Energy, th_H, bounds = bnds).x
    
    if fig:
        plt.figure()
        x = np.linspace(-pi, pi, 100)
        y = Energy(x)
        plt.plot(x,y)
    
    return th_M
 

def equilibriumAngle(M_s = 0.08, alpha = 0.0075, dims = [1.0, 1.0, 1.0], H_mag=0, th_H=0, eta = 1, I_dc = 0, demag = None, startingpt = None):
    """
    Calculates the equilibrium angle for a macrospin in an in-plane magnetic field.
    No DC current is applied in this case. Either the relative dimensions of the 
    nanowire or its demagnetization factors are required.
    """
    if demag == None: demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]
    
    vol  = 4*pi/3 * dims[0] * dims[1] * dims[2]
    kappa = 1/(vol) * hbar*gamma / (2*e_ch*M_s**2)
    lambd = alpha*gamma/M_s
    
    def equations(p):
        Mx, My = p 
        f1 = H_mag*np.cos(th_H) * My / (H_mag*np.sin(th_H) + kappa*eta*I_dc/lambd - Nyx * My) - Mx
        f2 = Mx**2 + My**2 - M_s**2
        return f1, f2
    
    if th_H == 0 and I_dc == 0:
        th_m = 0
    else:
        if startingpt == None: 
            th_est = H_mag*np.sin(th_H)/(H_mag*np.cos(th_H)+Nyx*M_s)            
            startingpt = (M_s*np.cos(th_est), M_s*np.sin(th_est))    
        M_x, M_y   = spop.fsolve(equations, startingpt, factor=1e-6);
        th_m = np.arctan2(M_y,M_x)

    th_m = np.mod(th_m, 2*pi);
    return float(th_m)
        

def eqAngle2(M_s = 0.08, alpha = 0.0075, dims = [1.0, 1.0, 1.0], H_mag=0, th_H=0, eta = 1, I_dc = 0, demag = None, nsteps = 36):
    """
    Calculates the equilibrium angle for a macrospin in an in-plane magnetic field.
    Supports cases including a dc current. Either the relative dimensions of the 
    nanowire or its demagnetization factors are required.
    """

    if demag == None: demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]    
   
    vol  = 4*pi/3 * dims[0] * dims[1] * dims[2]
    kappa = 1/(vol) * hbar*gamma / (2*e_ch*M_s**2)
    lambd = alpha*gamma/M_s
    
    def equations(p):
        Mx, My = p 
        f1 = H_mag*np.cos(th_H) * My - (H_mag*np.sin(th_H) + kappa*eta*I_dc/lambd - Nyx * My) * Mx
        f2 = Mx**2 + My**2 - M_s**2
        return f1, f2
        
    M_x = []; M_y = []; angs = np.linspace(0, 2*pi, nsteps+1)
    for k in angs:
         sol   = spop.root(equations, x0 = np.array([M_s*np.cos(k), M_s*np.sin(k)]), method = 'lm');       
         M_x.append(sol.x[0]); M_y.append(sol.x[1]);
       
    un, unin  = np.unique(np.round(M_x,5), return_index = True)
    sol = np.vstack((np.array(M_x)[unin], np.array(M_y)[unin]))
    
    # Filter the solution according to their norm (must be close to M_s)
    normsol = np.sqrt(sol[0,:]**2 + sol[1,:]**2)
    sol = sol[:, np.round(normsol,5) == np.round(M_s, 5)]
    
    # Compute the angle
    solang = np.arctan2(sol[1,:], sol[0,:])
    
    # Choose the solution that corresponds to the magnetic field quadrant, if we so please
    quad_H   = np.floor(th_H/(pi/2))
    solang = solang[np.floor( solang / (pi/2) ) == quad_H]

    if len(solang) == 0:
        solang = np.arctan2(sol[1,:], sol[0,:])
        solmin = np.argmin(np.abs(solang-th_H))
        solang = solang[solmin]
    else:
        solang = solang[0]
    
    return solang
    

def resonanceParameters(M_s = 8e5, alpha = 0.0075, dims = [1.0, 1.0, 1.0], H_mag = 0, th_H = 0, eta = 1, I_dc = 0, print_fields=False):
    """
    Calculates the resonance frequency and the linewidth for a given external field and DC current
    """
    # Calculate demagnetization factors and initialize variables   
    demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]
    
    Hox = H_mag * np.cos(th_H)
    Hoy = H_mag * np.sin(th_H)
    
    # Calculate the equilibrium angle   
    th0 = equilibriumAngleEnergy(M_s, alpha, dims, H_mag, th_H, eta, I_dc)
    #th = eqAngle2(M_s, alpha, dims, H_mag, th_H, eta, I_dc)    
    
    # Calculate the various magnetic field terms
    Hczx  = ( Hox + mu_0*Nzx*M_s*np.cos(th0) ) * np.cos(th0)
    Hszy  = ( Hoy + mu_0*Nzy*M_s*np.sin(th0) ) * np.sin(th0)
    Hcyx  = ( Hox + mu_0*Nyx*M_s*np.cos(th0) ) * np.cos(th0)
    Hsyx  = ( Hoy - mu_0*Nyx*M_s*np.sin(th0) ) * np.sin(th0)
    if print_fields: print Hczx, Hszy, Hcyx, Hsyx
    
    # Calculate the change in damping constant due to the dc current
    d_alpha = - 1 / (dims[0]*dims[1]*dims[2]) * hbar/(2*e_ch*M_s) * \
    (Hczx + Hszy + Hcyx + Hsyx) / ( (Hczx + Hszy)**2 + (Hcyx + Hsyx)**2 ) * eta * I_dc * np.sin(th0)
    
    # Calculate the frequency and the linewidth
    f0 = gamma * np.sqrt(  ( Hczx + Hszy ) * ( Hcyx + Hsyx ) ) / (2*pi)
    G =  gamma * (alpha + d_alpha) * np.sqrt( (Hczx + Hszy)**2 + (Hcyx + Hsyx)**2 ) / (2*pi)   
                                 
    return f0, G, th0, d_alpha



def resonanceVsHmag(M_s=8e5, alpha = 0.0075, dims = [1.0, 1.0, 1.0], th_H=0, eta = 1, I_dc=0, numsteps = 100, finalmag = 0.1, plot = True, plotNV = True):
    """
    Calculates and plots the resonance frequency and linewidth as a function of the magnitude of an external
    magnetic field. Accepts an arbitrary number of magnetic field angles.
    """
    # Calculate demagnetization factors and initialize variables   
    demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]

    H_mag = np.linspace(0, finalmag, numsteps)    
    
    try: iter(th_H)
    except TypeError: th_H = [th_H];
        
    # Calculate the equilibrium angle and the resonance frequency 
    th_m = []; f0 = []; G = []; d_alpha = []; cnt = 0
    for k in th_H:
        th_m.append([]); f0.append([]); G.append([]); d_alpha.append([]);
        
        for l in H_mag:
            f, g, th0, da = resonanceParameters(M_s = M_s, alpha = alpha, dims = dims, H_mag = l, th_H = k, eta = eta, I_dc = I_dc)
            th_m[cnt].append(th0)
            f0[cnt].append(f)
            G[cnt].append(g)
            d_alpha[cnt].append(da)
                                   
        cnt += 1    
        
    # Plot.    
    if plot == True:
        plt.figure(figsize = (15,5), dpi = 80); plt.hold(True)
        
        plt.subplot(141)
        plt.xlabel('Magnetic Field Magnitude (mT)'); plt.ylabel('Resonance Frequency (Hz)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        for m in range(0,cnt): 
            plt.plot(H_mag*1e3, f0[m]);
        if plotNV:
            x = np.linspace(0,finalmag,200)
            y1 = 2.87e9 + 28e9*x
            y2 = 2.87e9 - 28e9*x
            plt.plot(x*1e3, y1, '--k')
            plt.plot(x*1e3, y2, '--k')           
        plt.tight_layout()
        
        plt.subplot(142)
        plt.xlabel('Magnetic Field Magnitude (mT)'); plt.ylabel('Resonance Linewidth (Hz)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        for m in range(0,cnt): 
            plt.plot(H_mag*1e3, G[m]);
        plt.tight_layout()       
        
        plt.subplot(143)
        plt.xlabel('Magnetic Field Magnitude (mT)'); plt.ylabel('Magnetization Angle')
        for m in range(0,cnt): 
            #x = np.linspace(0, finalangle*180/pi, 200)
            #y = 180/pi * H_mag[m]*np.sin(x*pi/180)/(H_mag[m]*np.cos(x*pi/180) + Nyx*M_s)
            #plt.plot(x,y,'--k')            
            plt.plot(H_mag*1e3, np.array(th_m[m])*180/pi, label = str(th_H[m]*180/pi) + r'$^\circ$');
            plt.ylim([0, 90])
        plt.legend(prop={'size':12}, loc='best')
        plt.tight_layout()
                
        plt.subplot(144)
        plt.xlabel('Magnetic Field Magnitude (mT)'); plt.ylabel(r'$\sin(2\theta_0^M)$')
        for m in range(0,cnt): 
            #x = np.linspace(0, finalangle*180/pi, 200)
            #y = 180/pi * H_mag[m]*np.sin(x*pi/180)/(H_mag[m]*np.cos(x*pi/180) + Nyx*M_s)
            #plt.plot(x,y,'--k')            
            plt.plot(H_mag*1e3, np.sin(2.0*np.array(th_m[m])));
            plt.ylim([-1, 1])
        plt.legend(prop={'size':12}, loc='best')
        plt.tight_layout()       
        
    return th_H, th_m, f0, G, d_alpha
      
     
        
def resonanceVsThH(M_s=8e5, alpha = 0.0075, dims = [1.0, 1.0, 1.0], H_mag=0, eta = 1, I_dc=0, numsteps = 91, finalangle = pi/2-0.00001, plot = True):    
    """
    Calculates and plots the resonance frequency and linewidth as a function of the angle of an external
    magnetic field. Accepts an arbitrary number of magnetic field amplitudes.
    """
    # Calculate demagnetization factors and initialize variables   
    demag = demagFactors(dims);
    Nzx = demag[2]-demag[0]; Nyx = demag[1]-demag[0]; Nzy = demag[2]-demag[1]

    th_H = np.linspace(0, finalangle, numsteps)    
    
    try: iter(H_mag)
    except TypeError: H_mag = [H_mag];
          
    # Calculate the equilibrium angle and the resonance frequency 
    th_m = []; f0 = []; G = []; d_alpha = []; cnt = 0
    for k in H_mag:
        th_m.append([]); f0.append([]); G.append([]); d_alpha.append([]);
        
        for l in th_H:
            f, g, th, da = resonanceParameters(M_s = M_s, alpha = alpha, dims = dims, H_mag = k, th_H = l, eta = eta, I_dc = I_dc)
            th_m[cnt].append(th0)
            f0[cnt].append(f)
            G[cnt].append(g)
            d_alpha[cnt].append(da)
                                   
        cnt += 1
    
    # Plot.    
    if plot == True:
        plt.figure(figsize = (15,5), dpi = 80); plt.hold(True)
        
        plt.subplot(141)
        plt.xlabel('Magnetic Field Angle'); plt.ylabel('Resonance Frequency (Hz)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        for m in range(0,cnt): 
            plt.plot(th_H*180/pi, f0[m]);
        plt.tight_layout()
        
        plt.subplot(142)
        plt.xlabel('Magnetic Field Angle'); plt.ylabel('Resonance Linewidth (Hz)')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        for m in range(0,cnt): 
            plt.plot(th_H*180/pi, G[m]);
        plt.tight_layout()       
        
        plt.subplot(143)
        plt.xlabel('Magnetic Field Angle'); plt.ylabel('Magnetization Angle')
        for m in range(0,cnt): 
            #x = np.linspace(0, finalangle*180/pi, 200)
            #y = 180/pi * H_mag[m]*np.sin(x*pi/180)/(H_mag[m]*np.cos(x*pi/180) + Nyx*M_s)
            #plt.plot(x,y,'--k')            
            plt.plot(th_H*180/pi, np.array(th_m[m])*180/pi, label = str(np.round(H_mag[m]*1000,2)) + 'mT');
            plt.ylim([0, np.ceil(max(th_m[m])/(pi/2))*90])
        plt.legend(prop={'size':12}, loc='best')
        plt.tight_layout()
        
        plt.subplot(144)
        plt.xlabel('Magnetic Field Angle'); plt.ylabel(r'$\sin(2\theta_0^M)$')
        for m in range(0,cnt): 
            #x = np.linspace(0, finalangle*180/pi, 200)
            #y = 180/pi * H_mag[m]*np.sin(x*pi/180)/(H_mag[m]*np.cos(x*pi/180) + Nyx*M_s)
            #plt.plot(x,y,'--k')            
            plt.plot(th_H*180/pi, np.sin(2.0*np.array(th_m[m])));
            plt.ylim([-1, 1])
        plt.tight_layout()               
        
    return th_H, th_m, f0, G, d_alpha
    
    
def device_resistance(dims = [1.0, 1.0, 1.0], tpt = 0.01):
    
    rpt = 21.9*10**-6*10**4 * dims[0] / (dims[1] * tpt)
    rpy = 65.2*10**-6*10**4 * dims[0] / (dims[1] * dims[2])
    
    rtot = (1/rpt + 1/rpy)**-1
    
    return rtot
    