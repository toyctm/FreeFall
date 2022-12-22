#!/bin/env python

import thermo
from thermo import chemical as chem

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors
from matplotlib.ticker import StrMethodFormatter
from freefallmethod import *

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \cdot 10^{{{}}}$'.format(a, b)

def mfp(air):
    return np.sqrt(np.pi/8.)*air.mu/0.4987445*1/np.sqrt(air.P*air.rho) # Jennings 1988

def knudsen(air,d):
    return 2.*mfp(air)/d   # Seinfeld Eq. 8.1

def ccun(air,d):
    kn=knudsen(air,d)
    return 1. + kn*(1.257+0.4*np.exp(-1.1/kn))

def terminalspeed(air,d,g,rhop,deltafunc, freeslipcor=True):
   
   Vd = ( (rhop - air.rho) * g * d**2. ) / (18.*air.mu)  # Vd is the Stokes free fall velocity
   R = air.rho * Vd * d / (2.* air.mu)
   if(freeslipcor):
       Cc=ccun(air,d)
       R = R * Cc
       Vd = Vd* Cc
       
   delta = deltafunc(R)
   Vdcor = Vd * (1+delta)
   Re    = R  * (1+delta)
   return Vd,Vdcor,Re


cg = FreeFallMethod('Clift-Gauvin','cg',f_cliftgauvin)
stokes = FreeFallMethod('Stokes','stokes',f_stokes)

temp=298.15
pres=101325.

air=chem.Mixture('air',T=temp,P=pres)
print('density',air.rho)
print('pressure',air.P)
print('mu',air.mu)

dtup= [ 10**logd for logd in np.arange(-6.,-2.99,0.01) ]
dar= np.asarray(dtup)
printtup = (1.e-6, 1.e-5, 1.e-4)

rhop=2650.   # kg / m3 Quartz
g=9.81

Vdtup_cg=()
Vdtup_stokes=()
Vdtup_cgcc=()
Retup_cg=()
Retup_stokes=()
Retup_cgcc=()

for d in dtup:
   Vd, Vdcor_cgcc, Re_cgcc  = terminalspeed(air,d,g,rhop,cg.delta, freeslipcor=True)
   Vd, Vdcor_cg, Re_cg  = terminalspeed(air,d,g,rhop,cg.delta, freeslipcor=False)
   Vd, Vdcor_stokes, Re_stokes  = terminalspeed(air,d,g,rhop,stokes.delta, freeslipcor=False)

   Ccun = ccun(air,d)
   
   Vdtup_cg=Vdtup_cg + (Vdcor_cg,)
   Vdtup_stokes=Vdtup_stokes + (Vdcor_stokes,)
   Vdtup_cgcc=Vdtup_cgcc + (Vdcor_cgcc,)
   Retup_cg=Retup_cg + (Re_cg,)
   Retup_stokes=Retup_stokes + (Re_stokes,)
   Retup_cgcc=Retup_cgcc + (Re_cgcc,)

Vdstokes_ar=np.asarray(Vdtup_stokes)
Vdcg_ar=np.asarray(Vdtup_cg)
Vdcgcc_ar=np.asarray(Vdtup_cgcc)
   
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(dar*1.e6,Vdtup_cg,label='Clift-Gauvin')
plt.plot(dar*1.e6,Vdtup_cgcc,label='Clift-Gauvin+slip-correction')
plt.plot(dar*1.e6,Vdtup_stokes,label='Stokes')
plt.yscale('log')
plt.xscale('log')
ax.set_xlim([1., 1.e3])
plt.ylabel (r'$v_{\infty{}}\ (m/s)$')
plt.xlabel (r'$D\ (\mu{}m)$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('fig1_vinf.png',bbox_inches='tight',dpi=300)
plt.close()

fig,ax=plt.subplots(figsize=(4,3))
plt.plot(dar*1.e6,Retup_cg, label='Clift-Gauvin')
plt.plot(dar*1.e6,Retup_cgcc,label='Clift-Gauvin+slip-correction')
plt.plot(dar*1.e6,Retup_stokes, label='Stokes')
plt.yscale('log')
plt.xscale('log')
ax.set_xlim([1., 1.e3])
plt.ylabel ('Re')
plt.xlabel (r'$D\ (\mu{}m)$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.savefig('fig1_re.png',bbox_inches='tight',dpi=300)
plt.legend()
plt.close()



fig,ax=plt.subplots(figsize=(4,3))
plt.plot(dar*1.e6,(Vdcg_ar - Vdstokes_ar)/Vdstokes_ar,label='Clift-Gauvin')
plt.plot(dar*1.e6,(Vdcgcc_ar - Vdstokes_ar)/Vdstokes_ar,label='Clift-Gauvin+slip-correction')
plt.xscale('log')
ax.set_xlim([1., 1.e3])
ax.set_ylim([-1.03, 0.2])
plt.ylabel (r'$\left(v_{\infty{}}-v_{\infty{}}^{stokes}\right)/v_{\infty{}}^{stokes}$')
plt.xlabel (r'$D\ (\mu{}m)$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend(loc='lower left')
plt.savefig('fig1_corrections.png',bbox_inches='tight',dpi=300)
plt.close()


# Perform pressure-diameter speed diagram
# Set the list of virtual Reynolds numbers
rtup = [ 10**logr for logr in np.arange(-6.,5.,0.1) ]
rarray = np.asarray(rtup)
cg.evaluate(rarray)
cg.logistic_fit()


temp=298.15

dar    = np.asarray([ 10**logd for logd in np.arange(-6.,-3,0.1) ])
presar = np.linspace(20000,100000.,41)

npp=presar.shape[0]
nd=dar.shape[0]
Rear=np.zeros((nd,npp))
Vinf=np.zeros((nd,npp))
Vfit=np.zeros((nd,npp))
cc_ar=np.zeros((nd,npp))

for ip in range(npp):
   pres=presar[ip]
   air=chem.Mixture('air',T=temp,P=pres)
   print(pres)
   for id in range(nd):
      d = dar[id]
      Vd, vdcor, re  = terminalspeed(air,d,g,rhop,cg.delta, freeslipcor=True)
      Vd, vdcor_fit, re  = terminalspeed(air,d,g,rhop,cg.logistic_delta, freeslipcor=True)
      Cc = ccun(air,d)
      Rear[id,ip] = re
      Vinf[id,ip] = vdcor
      Vfit[id,ip] = vdcor_fit
      cc_ar[id,ip] = Cc
      
fig,ax=plt.subplots(figsize=(4,3))      
plt.contourf(presar/100.,dar*1.e6,Vinf, norm=colors.LogNorm(), levels=[1.e-5,1.e-4, 1.e-3, 1.e-2, 1.e-1, 1., 10., 14.])
plt.xlabel (r'$P\ (hPa)$')
plt.ylabel (r'$D\ (\mu{}m)$')
plt.colorbar(format=ticker.FuncFormatter(fmt),label=r'$v_\infty{}\mathrm{(\,m\,s^{-1}})$')
CS = plt.contour(presar/100.,dar*1.e6,Rear, levels=[0.1, 1., 10., 100., 1000.], colors='black')
CS2 = plt.contour(presar/100.,dar*1.e6,cc_ar, levels=[1.01, 1.1, 1.5, 2., 10.], colors='gray')
ax.clabel(CS, CS.levels, inline=True, fontsize=6)
ax.clabel(CS2, CS2.levels, inline=True, fontsize=6)
plt.yscale('log')
ax.set_yticks([1.,10.,100.,1000.])
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.savefig('fig3_vinf.png',bbox_inches='tight',dpi=300)
plt.close()

fig,ax=plt.subplots(figsize=(4,3))      
plt.contourf(presar/100.,dar*1.e6,(Vfit-Vinf)/Vinf*100, levels=[-2., -1., -0.5, -0.1, 0.1, 0.5, 1., 2.], cmap='bwr')
plt.xlabel (r'$P\ (hPa)$')
plt.ylabel (r'$D\ (\mu{}m)$')
plt.colorbar(label='Fit error (%)')
plt.yscale('log')
ax.set_yticks([1.,10.,100.,1000.])
plt.grid(color='black', linestyle=':', linewidth=0.5)
plt.savefig('fig3_fiterror.png',bbox_inches='tight',dpi=300)
plt.close()

      
