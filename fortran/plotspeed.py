#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

expintervals=np.asarray([-7.,-6.,-5.,-4.,-3.])
expcenters=10.**(np.asarray([-6.5, -5.5, -4.5, -3.5]))

niters={}
niters['fixpoint']=[ 0., 0., 1.23, 7.11 ]
niters['bisection']=[ 0., 0., 3.38, 7.47 ]

exetime={}
exetime['fixpoint']=[ 47.35, 47.41, 84.62, 291. ]
exetime['bisection']=[ 48.32, 48.27, 135.63, 260. ]
exetime['aersett']=[ 45.92, 46.03, 42.76, 34.15 ]

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
      vpos=avgdict[method]+4.
  if(method=='fixpoint'):
      vpos=avgdict[method]-18.

  plt.text(10**-4., vpos, "avg={:.1f} ns".format(avgdict[method]), color=colordict[method], fontsize='x-small')
plt.legend(loc='upper left')
plt.xscale('log')
ax.set_xlim([1.e-7,1.e-3])
ax.set_ylim([0.,500.])
plt.ylabel ('Execution time (ns)')
plt.xlabel (r'$D\ (\mu{}\mathrm{m})$')
plt.savefig('exetime.png',bbox_inches='tight',dpi=300)
