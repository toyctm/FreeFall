#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

expintervals=np.asarray([0.,1.,2.,3.,4.])
expcenters=10.**(np.asarray([.5, 1.5, 2.5, 3.5]))

niters={}
niters['fixpoint']=[ 3.40601826, 5.91496563, 9.50085926, 14.9136591 ]
niters['bisection']=[ 6.24640656, 7.11500883, 8.24732399, 9.46562958 ]

exetime={}
exetime['fixpoint']=[ 96.719999999999999, 177.25000000000000, 288.95999999999998, 455.00999999999999 ]
exetime['bisection']=[ 151.72000000000000, 176.13000000000000, 198.63000000000000, 216.66999999999999 ]
exetime['aersett']=[ 17.800000000000001, 18.109999999999999, 18.270000000000000, 18.260000000000002 ]

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
  plt.plot([10**0., 10**4], [avgdict[method], avgdict[method]], color=colordict[method], linestyle=':', linewidth=1.)
  plt.text(10**3., avgdict[method]+4., "avg={:.1f} ns".format(avgdict[method]), color=colordict[method], fontsize='x-small')
plt.legend(loc='upper left')
plt.xscale('log')
ax.set_xlim([1.,1.e+4])
ax.set_ylim([0.,500.])
plt.ylabel ('Execution time (ns)')
plt.xlabel (r'$R$')
plt.savefig('exetime.png',bbox_inches='tight',dpi=300)
