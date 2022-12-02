#!/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

expintervals=np.asarray([0.,1.,2.,3.,4.])
expcenters=10.**(np.asarray([.5, 1.5, 2.5, 3.5]))

niters={}
niters['fixpoint']=[ 3.51564956, 6.18300438, 9.87728405, 15.4518681 ]
niters['bisection']=[ 6.83266497, 7.42684412, 8.51116753, 9.70776367 ]

exetime={}
exetime['fixpoint']=[ 103.2399999999999, 188.25999999999999, 309.86000000000001, 473.31000000000000 ]
exetime['bisection']=[ 167.60000000000002, 184.69000000000000, 210.28000000000000, 226.96000000000001 ]
exetime['aersett']=[ 18.149999999999999, 18.039999999999999, 18.410000000000000, 18.780000000000001]

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
  plt.plot(expcenters, exetime[method], label=labeldict[method], color=colordict[method])
  plt.plot([10**0.5, 10**3.5], [avgdict[method], avgdict[method]], color=colordict[method], linestyle=':', linewidth=1.)
  plt.text(10**3., avgdict[method]+4., "avg={:.1f} ns".format(avgdict[method]), color=colordict[method], fontsize='x-small')
plt.legend()
plt.xscale('log')
ax.set_xlim([1.,1.e+4])
ax.set_ylim([0.,500.])
plt.savefig('exetime.png',bbox_inches='tight',dpi=300)
