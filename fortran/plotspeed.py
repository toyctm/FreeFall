#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

expintervals=np.asarray([-7.,-6.,-5.,-4.,-3.])
expcenters=10.**(np.asarray([-6.5, -5.5, -4.5, -3.5]))

niters={}
niters['fixpoint']=[ 0., 0., 1.40372717, 7.11220407 ]
niters['bisection']=[ 0., 0., 4.41905737, 7.47149038 ]

exetime={}
exetime['fixpoint']=[ 12.400000000000000, 12.699999999999999, 69.000000000000000, 272.69999999999999 ]
exetime['bisection']=[ 12.480000000000000, 12.730000000000000, 132.16000000000000, 226.11000000000001 ]
exetime['aersett']=[ 12.529999999999999, 12.890000000000001, 38.549999999999997, 35.579999999999998 ]

labeldict={}
labeldict['fixpoint']='Fixed-point method'
labeldict['bisection']='Bisection method'
labeldict['aersett']='AerSett'

colordict={}
colordict['fixpoint']='red'
colordict['bisection']='orange'
colordict['aersett']='green'

avgdict={}
for method in ('aersett', 'bisection', 'fixpoint'):
    avgdict[method]=np.sum(np.asarray(exetime[method]))/4.

fig,ax=plt.subplots(figsize=(4,3))
for method in ('aersett', 'bisection', 'fixpoint'):
  plt.stairs(edges=10**expintervals, values=exetime[method], label=labeldict[method], color=colordict[method])
  plt.plot([10**-7., 10**-3], [avgdict[method], avgdict[method]], color=colordict[method], linestyle=':', linewidth=1.)
  if(method=='bisection'):
      vpos=avgdict[method]+4.
  if(method=='aersett'):
      vpos=avgdict[method]-12.
  if(method=='fixpoint'):
      vpos=avgdict[method]-12.

  plt.text(1.3e-4, vpos, "avg={:.1f} ns".format(avgdict[method]), color=colordict[method], fontsize='x-small')
plt.legend(loc='upper left')
plt.xscale('log')
ax.set_xlim([1.e-7,1.e-3])
ax.set_ylim([0.,300.])
plt.ylabel ('Execution time (ns)')
plt.xlabel (r'$D\ (\mathrm{m})$')
plt.savefig('exetime.png',bbox_inches='tight',dpi=300)
