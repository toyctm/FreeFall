#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

expintervals=np.asarray([-7.,-6.,-5.,-4.,-3.])
expcenters=10.**(np.asarray([-6.5, -5.5, -4.5, -3.5]))

niters={}
niters['fixpoint']=[ 0., 0., 1.47902966, 7.40053558 ]
niters['bisection']=[ 0., 0., 4.45684910, 7.83219337 ]

exetime={}
exetime['fixpoint']=[ 12.070000000000000, 12.410000000000000, 69.890000000000001, 281.07999999999998 ]
exetime['bisection']=[ 12.100000000000000, 12.369999999999999, 129.59000000000000, 232.86000000000001 ]
exetime['aersett']=[ 11.990000000000000, 12.199999999999999, 37.649999999999999, 35.159999999999997 ]

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
