#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import array
import numpy
from numpy.linalg import norm
from pylab import *

from scipy.integrate import simps
from numpy import sinh, tanh, log
from scipy import optimize

from scipy.constants import Boltzmann, Avogadro, electron_volt
from ase.units import Rydberg,  eV, GPa
from pyspglib import spglib
from scipy.integrate import simps
from numpy import sinh, tanh, log

k_B=Boltzmann/electron_volt
THz2meV=1/0.241799
cminv2meV=1/8.0655

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
                    exper=None, ax=None, label=None, positive_only=True, **kwargs):

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
            ax2.plot(freq[0], freq[i+1],'-',color=clr)
        else :
            clr=ax2.plot(freq[0], freq[i+1],'-',label=label)[0].get_color()
        
        
    if exper :
        clr=ax2.plot(exper[0],exper[1],'o',label='Experiment')[0].get_color();

    
    xlim(0,max(freq[0]))
    
    if positive_only :
        ylim(0,ylim()[1])
    else :
        ylim(ylim()[0],ylim()[1])

    if qpname :
        #print t, qpname, min(freq[0]), max(freq[0])
        xticks(t,qpname)
        vlines(xticks()[0][1:-1],ylim()[0],ylim()[1],linestyles='--')
        
    subplots_adjust(wspace=0)
    return [ax1,ax2]
    

def plot_bands(result=None, show_gap=True,
                    qpath=array([]), qpname=[], 
                    ax=None, label=None, **kwargs):
    
    assert(result)
    kpoints=result['bands_kpt']
    energies=result['bands']
    edos=result['edos']
    pk=kpoints[0]
    kpth=[0]
    for k in kpoints:
        kpth.append(norm(k-pk)+kpth[-1])
        pk=k
    kpth=array(kpth[1:])

    ax1,ax2=plot_phonons(insert(energies,0,kpth,axis=0),edos,
                qpath=qpath,qpname=qpname,ax=ax,label=label,positive_only=False, **kwargs)
    
    if show_gap and 'HOL' in result.keys() and 'LUL' in result.keys():
        hol=result['HOL']
        lul=result['LUL']
        sca(ax2)
        axhspan(hol,lul,alpha=0.15, color='k')
        sca(ax1)
        axhspan(hol,lul,alpha=0.15, color='k')

    return [ax1,ax2]


def get_thermodynamic_functions(dos,Tmin=None,Tmax=1500,Tstep=10):
    '''
    Calculate free energy contribution from phonons 
    and heat capacity as a function of temperature using PDOS.
    '''
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
    


def analyze_QHA_run(Calc, Tmin=1, Tmax=1500, Tsteps=100):
    '''
    Run the QHA analysis procedure on the set of volume calculation Calc.
    This procedure produces an array with thermodynamic functions for
    all calculated volumes and given temperature range (the zero temperature is always included).
    
    Input
    =====
    
        Calc   - List of calculation data
        Tmin   - Minimum temperature of the scan (K)
        Tmax   - Maximum temperature of the scan (K)
        Tsteps - Number of temperature steps
        bdir   - base directory (subdirectories host and calculation are located below)
        
    Output
    ======
        
        qha  - array of N x Tsteps x #params containing resulting data:
               V, T, Etot, Fph, P*V, Cv, S, ... 
               for each calculation and each temperature step.
    '''
    phdos={}
    qha=[]
    
    for n,c in enumerate(Calc):
        PCMUL=len(c.get_atomic_numbers())/len(spglib.find_primitive(c)[2])
        # Read the calculated dos
        clc=c.calc.results
        phdos[n]=clc['phdos']
        #Etot=c.get_potential_energy()
        Etot=clc['etotal']*Rydberg

        # Pressure returned in kbar we need to change it to eV/A^3 to get p*V in eV
        # eV/A^3=160.2176487GPa
        # kbar=0.1 GPa
        # eV/A^3 = 1602.176487 kbar
        P=clc['pressure']/1602.176487
        #P=c.get_isotropic_pressure(c.get_stress())

        # Scan the temperatures at this volume 
        # Unit cell (primitive) volume in Bohr^3 => (A^3)
        #V=(calcQHA[idx][2]['A'])**3
        
        #V=clc['volume']*(Bohr**3)
        V=c.get_volume()/PCMUL
    
        # get the ndf from the normalization of dos
        dos=phdos[n]
        ndf=round(simps(dos[1],dos[0]))
        
        # frequencies/energies
        # QE outputs dos in cm^-1, let's convert to sane units 
        # Let's convert to energy to include the hbar
        nu=1e-3*cminv2meV*phdos[n][0]
        dos=phdos[n][1]
        
        # We need to cut the nu to positive range
        no=array([[o,g] for o,g in zip(nu,dos) if o>0 ])
        nu=no[:,0]
        dos=no[:,1]
    
        # correct the normalization for units, negative cut and numerical errors
        dos=dos/simps(dos,nu)
    
        # plot in meV
        #plot(1e-3*nu,dos,label=n);
        #legend()
        #show()
        
        # plot the integrand at 300K
        #T=300
        #plot(nu,dos*log(2*sinh(0.5*nu/(k_B*T))),label=idx)
        
        # Zero-point energy - $\hbar\omega/2$ for each degree of freedom integrated over $\omega$
        F0=ndf*simps(0.5*dos*nu,nu)
    
        # Put in special case for the T=0. Just ZPV energy.
        qhav=[[V,0,Etot,F0,0,0,0,0]]
        qha.append(qhav)
        print '# %s: V=%.2f A^3, Etot=%.3f eV, P=%6.1f GPa' % (n, V, Etot, P/GPa)
        for T in linspace(Tmin,Tmax,Tsteps):
            # Maradudin book formula - the only one important for thermal expansion
            # The result is in eV, temperature in K
            e=0.5*nu/(k_B*T)
            Fph=ndf*T*k_B*simps(dos*log(2*sinh(e)),nu)
            Cv=ndf*k_B*simps(e*e*dos/sinh(e)**2,nu)
            S=k_B*(simps(2*dos*e/(exp(2*e)-1),nu)-simps(dos*(1-exp(-2*e)),nu))
            # Alternative formula from QE
            Sqe=k_B*simps(dos*(e/tanh(e)-log(2*sinh(e))),nu)
            # Formulas lifted from the QHA code in QE - not correct with current units
            #q=0.5*a3*nu/T
            #E_int=simps(a1*dos*nu/tanh(q),nu)
            #S0=simps(dos*(q/tanh(q)-log(2*sinh(q))),nu)
            qhav.append([V, T, Etot, Fph, P*V, Cv, S, Sqe])
            #print " %5.2f  %12f  %12f  %12f" % (T, Cv, S, Sqe)
        
    return array(qha)
    #legend();


def bm_eos(v,p):
    e0,v0,b0,b0p=p
    x=(v0/v)**(2.0/3)
    return e0+(9*v0*b0/16)*(b0p*(x - 1)**3+(6-4*x)*((x-1)**2))

def fit_and_plot_QHA(qha,P=0.0,nop=20, eos_fun=bm_eos):
    '''
    Run a QHA fitting procedure on data generated by analyze_QHA_run.
    Runs the fits and plots nop temperature lines with fitting results.
    Returns fitted parameters along the isobar. 
    
    Input
    =====
    
        qha - array produced by analyze_QHA_run function
        nop - number of temperature lines ploted. The first and last temperature is alwayes included.
    
    Output
    ======
        
        array containing: T, V(T), E(T), B(T), B'(T) along the isobar
    
    '''
    
    eos_err = lambda p, x, y: eos_fun(x,p) - y
    V0=[]
    plots=range(0,qha.shape[1],qha.shape[1]//nop)+[qha.shape[1]-1]
    # Plots in lattice constants
    for t in range(0,qha.shape[1]):
        #e=qha[:,t,4]+qha[:,t,3]+qha[:,t,2]-qha[:,t,5]
        # Calculate E_tot+F_ph
        v=qha[:,t,0]
        G=qha[:,t,2]+qha[:,t,3]+P*v
        if t in plots :
            plot(v,G,'+',label="T=%.0fK"%(qha[0,t,1],));
        p0=[-1800,(max(v)+min(v))/2,200,1]
        fit,succ=optimize.leastsq(eos_err, p0[:], args=(v,G))
        if not succ : 
            print succ
            continue
        xmi=min(v)*0.995
        xma=max(v)*1.005
        x=linspace(xmi, xma, 100)
        if t in plots:
            plot(x,eos_fun(x,fit),'-')
            plot(fit[1],fit[0],'r.')
        # Add the fitting results T, V(T), E(T), B(T), B'(T) to the collection
        V0.append([qha[0,t,1],fit[1],fit[0],fit[2],fit[3]])
        #print succ, fit
    xlabel('Primitive Cell volume ($\AA^3$)')
    ylabel('Gibbs free energy (eV)')
    legend(loc='upper right')
    V0=array(V0).T
    
    # Fitted equilibrium volumes
    plot(V0[1],V0[2],'r-');
    
    return copy(V0)

    
#def analyze_QHA_run(Calc, Tmin=1, Tmax=1500, Tsteps=100):
#    '''
#    Run the QHA analysis procedure on the set of volume calculation Calc.
#    This procedure produces an array with thermodynamic functions for
#    all calculated volumes and given temperature range (the zero temperature is always included).
#    
#    Input
#    =====
#    
#        Calc   - List of calculation data
#        Tmin   - Minimum temperature of the scan (K)
#        Tmax   - Maximum temperature of the scan (K)
#        Tsteps - Number of temperature steps
#        bdir   - base directory (subdirectories host and calculation are located below)
#        
#    Output
#    ======
#        
#        qha  - array of N x Tsteps x #params containing resulting data:
#               V, T, Etot, Fph, P*V, Cv, S, ... 
#               for each calculation and each temperature step.
#    '''
#    phdos={}
#    qha=[]
#    
#    for idx in [k for k in sorted(Calc.keys())]:
#        # Read the calculated dos
#        host=Calc[idx][0]
#        try :
#            phdos[idx]=loadtxt(bdir + host + idx +'/'+Calc[idx][2]['prefix']+'.dos').T
#            clc=qeio.read_quantumespresso_textoutput(bdir + host + idx + '/pw.out')
#            # QE returns Ry, now in eV
#            Etot=clc['etotal']*Ry
#            # Pressure returned in kbar we need to change it to eV/A^3 to get p*V in eV
#            # eV/A^3=160.2176487GPa
#            # kbar=0.1 GPa
#            # eV/A^3 = 1602.176487 kbar
#            P=clc['pressure']/1602.176487
#        except IOError :
#            print "No DOS for: ",idx
#            continue

#        # Scan the temperatures at this volume 
#        # Unit cell (primitive) volume in Bohr^3 => (A^3)
#        #V=(calcQHA[idx][2]['A'])**3
#        V=clc['volume']*(Bohr**3)
#    
#        # get the ndf from the normalization of dos
#        dos=phdos[idx]
#        ndf=round(simps(dos[1],dos[0]))
#        
#        # frequencies/energies
#        # QE outputs dos in cm^-1, let's convert to sane units 
#        # Let's convert to energy to include the hbar
#        nu=cminv2meV*phdos[idx][0]
#        dos=phdos[idx][1]
#    
# 
#        # Now it is in meV switch to eV
#        nu=1e-3*nu
#    
#        # We need to cut the nu to positive range
#        no=array([[o,g] for o,g in zip(nu,dos) if o>0 ])
#        nu=no[:,0]
#        dos=no[:,1]
#    
#        # correct the normalization for units, negative cut and numerical errors
#        dos=dos/simps(dos,nu)
#    
#        # plot in meV
#        #plot(1e-3*nu,dos,label=idx);
#        #legend()
#        #show()
#        
#        # plot the integrand at 300K
#        #T=300
#        #plot(nu,dos*log(2*sinh(0.5*nu/(k_B*T))),label=idx)
#        
#        # Zero-point energy - $\hbar\omega/2$ for each degree of freedom integrated over $\omega$
#        F0=ndf*simps(0.5*dos*nu,nu)
#    
#        # Put in special case for the T=0. Just ZPV energy.
#        qhav=[[V,0,Etot,F0,0,0,0,0]]
#    
#        qha.append(qhav)
#        print '# %s: V=%f A^3, Etot=%f eV, P=%f eV/A^3' % (idx, V, Etot, P)
#        for T in linspace(Tmin,Tmax,Tsteps):
#            # Maradudin book formula - the only one important for thermal expansion
#            # The result is in eV, temperature in K
#            e=0.5*nu/(k_B*T)
#            Fph=ndf*T*k_B*simps(dos*log(2*sinh(e)),nu)
#            Cv=ndf*k_B*simps(e*e*dos/sinh(e)**2,nu)
#            S=k_B*(simps(2*dos*e/(exp(2*e)-1),nu)-simps(dos*(1-exp(-2*e)),nu))
#            # Alternative formula from QE
#            Sqe=k_B*simps(dos*(e/tanh(e)-log(2*sinh(e))),nu)
#            # Formulas lifted from the QHA code in QE - not correct with current units
#            #q=0.5*a3*nu/T
#            #E_int=simps(a1*dos*nu/tanh(q),nu)
#            #S0=simps(dos*(q/tanh(q)-log(2*sinh(q))),nu)
#            qhav.append([V, T, Etot, Fph, P*V, Cv, S, Sqe])
#            #print " %5.2f  %12f  %12f  %12f" % (T, Cv, S, Sqe)
#        
#    return array(qha)
#    #legend();
