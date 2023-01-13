#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from scipy.optimize import curve_fit
from freefallmethod import *
from scipy.optimize import fsolve


# Set the list of virtual Reynolds numbers
rtup = [ 10**logr for logr in np.arange(-6.,5.,0.1) ]
rarray = np.asarray(rtup)


# Initialize the Clift-Gauvin method
cg=FreeFallMethod('Clift-Gauvin', 'cg', f_cliftgauvin)
fb=FreeFallMethod('Flemmer-Banks', 'fb', f_flemmerbanks)
cheng=FreeFallMethod('Cheng', 'cheng', f_cheng)
methodtup=(cg,cheng)

for method in methodtup:
  print('Evaluate for method',method.name)
  method.evaluate(rarray)
  
fig,ax=plt.subplots(figsize=(4,3))
for method in methodtup:
   plt.plot(method.rarray, [12./method.rarray[i]*(1+method.f(method.rarray[i])) for i in range(method.n) ],label=method.name, linewidth=0.5)
plt.xscale('log')
plt.yscale('log')
ax.set_xlim([1.e-5,1.e+5])
plt.legend()
plt.savefig('fig2_Cd.png',bbox_inches='tight',dpi=300)

fig,ax=plt.subplots(figsize=(4,3))
for method in methodtup:
   plt.plot(method.rarray, method.deltavec,label=method.name, linewidth=2.)
plt.xscale('log')
plt.ylabel (r'$\delta{}\left(R\right)$')
plt.xlabel (r'$R$')
ax.set_xlim([1.e-5,1.e+5])
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('fig2_delta.png',bbox_inches='tight',dpi=300)

for method in methodtup:
  print('Fitting for method',method.name)
  method.logistic_fit()

def f(x):
  return cg.logistic_delta(x)+0.01
Re_0 = fsolve(f, 0.03)
print('==============================================')
print('Virtual Reynolds number for delta=0.01', Re_0)
print('double-check', cg.logistic_delta(Re_0))
print('==============================================')





# Plot solutions with fit done on the virtual Reynolds number
for method in methodtup:
   fig,ax=plt.subplots(figsize=(4,3))
   plt.plot(method.rarray, method.deltavec,label=method.name, linewidth=2.)
   plt.plot(method.rarray, method.fitdelta, color='red',label='Generalized logistic', linewidth=1.0)
   plt.xscale('log')
   plt.ylabel (r'$\delta{}\left(R\right)$')
   plt.xlabel (r'$R$')
   ax.set_xlim([1.e-5,1.e+5])
   plt.grid(color='b', linestyle=':', linewidth=1)
   plt.legend()
   plt.savefig('fig2_delta_'+method.shortname+'.png',bbox_inches='tight',dpi=300)
   plt.close()

# Plot solutions with fit done on the Reynolds number
   fig,ax=plt.subplots(figsize=(4,3))
   plt.plot(method.reyarray, method.deltavec,label=method.name, linewidth=1.)
   plt.plot(method.reyarray, method.fitdelta, color='black',label='Generalized logistic', linewidth=0.5)
   plt.plot(method.reyarray, method.reyfitdelta, color='red',label='Fitted on Re', linewidth=0.5)
   plt.xscale('log')
   plt.ylabel (r'$\delta{}\left(Re\right)$')
   plt.xlabel (r'$Re$')
   ax.set_xlim([1.e-5,1.e+3])
   plt.grid(color='b', linestyle=':', linewidth=1)
   plt.legend()
   plt.savefig('fig2_delta_Re_'+method.shortname+'.png',bbox_inches='tight',dpi=300)
   plt.close()


# Plot vinfty error with Reynolds
   fig,ax=plt.subplots(figsize=(4,3))
   plt.plot(method.reyarray,method.fiterror,label='Fit error')
   plt.xscale('log')
   plt.ylabel (r'$\left(v_{\infty{}}^{fit}-v_{\infty{}}\right)/v_{\infty{}}\ (\%)$')
   plt.xlabel (r'$Re$')
   ax.set_xlim([1.e-5,1.e+3])
   ax.set_ylim([-3.5,3.5])
   plt.grid(color='b', linestyle=':', linewidth=1)
   plt.legend()
   plt.savefig('fig2_vinfty_error_'+method.shortname+'.png',bbox_inches='tight',dpi=300)
   plt.close()

# Plot vinfty model discrepancy with Reynolds
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(cg.reyarray,(cg.deltavec-cheng.deltavec)/(1.+cheng.deltavec)*100,label='CG-Cheng discrepancy')
# plt.plot(cg.reyarray,(cg.fitdelta-cheng.deltavec)/(1.+cheng.deltavec)*100,label='CG(fit)-Cheng discrepancy')
plt.xscale('log')
   #plt.yscale('log')
plt.ylabel (r'$\left(v_{\infty{}}^{cg}-v_{\infty{}}^{Cheng}\right)/v_{\infty{}^{Cheng}}\ (\%)$')
plt.xlabel (r'$Re$')
ax.set_xlim([1.e-5,1.e+3])
ax.set_ylim([-3.5,3.5])
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('fig2_vinfty_discrepancy.png',bbox_inches='tight',dpi=300)
plt.close()

# Compare our method to the method of Cheng 2009
cheng.evaluate_chengwstar()

# Plot wstar in our method and in the Cheng formula
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(cheng.reyarray,cheng.wstar,label='Numerical solution')
plt.plot(cheng.reyarray,cheng.wstarcheng,label='Cheng formula', linewidth=0.5)
plt.plot(cheng.reyarray,cheng.wstarfit,label='Logistic formula', linewidth=0.5)
plt.xscale('log')
plt.ylabel (r'$w^{*}$')
plt.xlabel (r'$Re$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('fig2_wstar_cheng.png',bbox_inches='tight',dpi=300)
plt.close()

# Plot wstar error in our method and in the Cheng formula
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(cheng.reyarray,(cheng.wstarcheng-cheng.wstar)/cheng.wstar*100,label='Cheng formula error')
plt.plot(cheng.reyarray,(cheng.wstarfit-cheng.wstar)/cheng.wstar*100,label='Logistic fit error')
plt.xscale('log')
plt.ylabel (r'$w^{*}\ \mathrm{error\ (\%)}$')
plt.xlabel (r'$Re$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('fig2_wstar_error.png',bbox_inches='tight',dpi=300)
plt.close()

# Compare to fig2 of Cheng
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(cheng.dstar,cheng.cdcheng,label='Cheng formula')
plt.xscale('log')
plt.yscale('log')
plt.ylabel (r'$C_D$')
plt.xlabel (r'$d_{*}$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('chengfig2.png',bbox_inches='tight',dpi=300)
plt.close()

# Just reproduce fig2 of cheng
dvec = np.asarray([ 10**logr for logr in np.arange(-.5,4.,0.1) ])
fig,ax=plt.subplots(figsize=(4,3))
plt.plot(dvec,432/dvec**3*(1+0.022*dvec**3)**0.54+0.47*(1-np.exp(-0.15*dvec**0.45)),label='Cheng eq. 2')
plt.xscale('log')
plt.yscale('log')
plt.ylabel (r'$C_D$')
plt.xlabel (r'$d_{*}$')
plt.grid(color='b', linestyle=':', linewidth=1)
plt.legend()
plt.savefig('chengfig2-plain.png',bbox_inches='tight',dpi=300)
plt.close()
