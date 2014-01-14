#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import array
import numpy
from scipy.constants import Boltzmann, eV, Avogadro
from numpy.linalg import norm
from pylab import *

from scipy.integrate import simps
from numpy import sinh, tanh, log


THz2meV=1/0.241799
cminv2meV=1/8.0655
k_B=Boltzmann/eV

def get_EOS(d, comment=""):
    # Fitting functions
    def BMEOS(v,v0,b0,b0p):
        return (b0/b0p)*(pow(v0/v,b0p) - 1)
    fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]
    errfunc = lambda p, x, y: fitfunc(p, x) - y
    
    pv=array([d[:,0]**3,d[:,1]]).T
    
    # Estimate the initial guess assuming b0p=1
    # Limiting volumes
    v1=min(pv[:,0])
    v2=max(pv[:,0])
    # The pressure is falling with the growing volume
    p2=min(pv[:,1])
    p1=max(pv[:,1])
    b0=(p1*v1-p2*v2)/(v2-v1)
    v0=v1*(p1+b0)/b0
    # Initial guess
    p0=[v0,b0,1]
    #Fitting
    #print p0
    fit, succ = optimize.leastsq(errfunc, p0[:], args=(pv[:,0],pv[:,1]))
    
    # Ranges - the ordering in pv is not guarateed at all!
    # In fact it may be purely random.
    x=numpy.array([min(pv[:,0]),max(pv[:,0])])
    y=numpy.array([min(pv[:,1]),max(pv[:,1])])
    
    # Plot the P(V) curves and points for the crystal
    # Plot the points
    plot(pv[:,0]**1/3,pv[:,1],'.',label='Calc. '+comment)
    
    # Mark the center P=0 V=V0
    axvline(fit[0]**1/3,ls='--')
    axhline(0,ls='--')
    
    # Plot the fitted B-M EOS through the points
    xa=numpy.linspace(x[0],x[-1],20)
    plot(xa**1/3,fitfunc(fit,xa),'-', label="B-M fit:\n$V_0$=%f ($A_0$=%f),\n$B_0$=%f kBar,\n$B'_0$=%f  " 
         % (fit[0], fit[0]**(1.0/3), fit[1], fit[2]) )
    legend()
    ylabel('Pressure (kBar)')
    xlabel('Lattice constant ($\mathrm{\AA}$)')
    draw();
    return fit, pv


def plot_phonons(freq=None, dos=None, 
                    qpath=array([]), qpname=[], 
                    exper=None, ax=None, label=None, **kwargs):

#idx, ex=True, bdir='/home/jochym/Desktop/Fizyka/', qp=array([G1,X,G2,L]), lbl=None, ax=None, castep=False):
    s=0
    t=[]
    for x in [0]+map(norm,qpath[1:]-qpath[:-1]):
        s+=x
        t.append(s)
    
    if ax :
        ax1=ax[0]
    else :
        ax1=subplot2grid([1,4],[0,3],colspan=1)
    
    clr=ax1.plot( dos[1], dos[0])[0].get_color()
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.xaxis.tick_top()
    
    if ax :
        ax2=ax[1]
    else :
        ax2=subplot2grid([1,4],[0,0],colspan=3,sharey=ax1)
    
    clr=None
    
    for i in range(freq.shape[0]-1):
        if clr :
            ax2.plot(freq[0,1:], freq[i+1,1:],'.',color=clr)
        else :
            clr=ax2.plot(freq[0,1:], freq[i+1,1:],'.',label=label)[0].get_color()
        
        
    if exper :
        clr=ax2.plot(exper[0],exper[1],'o',label='Experiment')[0].get_color();

    
    xlim(0,max(freq[0,1:]))
    ylim(0,ylim()[1])
    if qpname :
        #print t, qpname, min(freq[0]), max(freq[0])
        xticks(t,qpname)
        vlines(xticks()[0][1:-1],ylim()[0],ylim()[1],linestyles='--')
        
    subplots_adjust(wspace=0)
    return [ax1,ax2]
    
def get_thermodynamic_functions(dos,Tmin=None,Tmax=1500,Tstep=10):
        # If Tmin=None get the step size
        if Tmin==None :
            Tmin=Tstep
        # get the ndf from the normalization of dos
        ndf=round(simps(dos[1],dos[0]))
        
        # frequencies/energies
        # QE outputs dos in cm^-1, let's convert to sane units 
        # Let's convert to energy (eV) to include the hbar
        nu=1e-3*cminv2meV*dos[0]
        dos=dos[1]
        
        # We need to cut the nu to positive range
        no=array([[o,g] for o,g in zip(nu,dos) if o>0 ])
        nu=no[:,0]
        dos=no[:,1]
    
        # correct the normalization for units, negative cut and numerical errors
        dos=dos/simps(dos,nu)
    
        # Zero-point energy - $\hbar\omega/2$ for each degree of freedom integrated over $\omega$
        F0=ndf*simps(0.5*dos*nu,nu)
    
        # Put in special case for the T=0. Just ZPV energy.
        # T, Free energy, Cv, Cp, S
        tfun=[[0,F0,0]]
        
        for T in arange(Tmin,Tmax,Tstep):
            # Maradudin book formula - the only one important for thermal expansion
            # The result is in eV, temperature in K
            e=0.5*nu/(k_B*T)
            Fph=ndf*T*k_B*simps(dos*log(2*sinh(e)),nu)
            Cv=eV*Avogadro*ndf*k_B*simps(e*e*dos/sinh(e)**2,nu)
            #S=k_B*(simps(2*dos*e/(exp(2*e)-1),nu)-simps(dos*(1-exp(-2*e)),nu))
            # Alternative formula from QE
            #Sqe=k_B*simps(dos*(e/tanh(e)-log(2*sinh(e))),nu)
            # Formulas lifted from the QHA code in QE - not correct with current units
            #q=0.5*a3*nu/T
            #E_int=simps(a1*dos*nu/tanh(q),nu)
            #S0=simps(dos*(q/tanh(q)-log(2*sinh(q))),nu)
            tfun.append([T, Fph, Cv])
            #print " %5.2f  %12f  %12f  %12f" % (T, Cv, S, Sqe)
        return array(tfun).T
