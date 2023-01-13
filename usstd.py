#!/bin/env python

def usstd(zvals):

  zbot = [0.      , 11000. , 20000. , 32000. , 47000. , 51000.  , 71000. , 100000.]
  lr   = [-0.0065 , 0.0    , 0.001  , 0.0028 , 0.0    , -0.0028 , -0.002 , 0.]

  tprev= 288.15
  pprev= 101325.
  ptup = (pprev,)
  ttup = (tprev,)
  zprev=zvals[0]
  g=9.81
  R=287.1
  islice=0
  for z in zvals:
      if (z >= zbot[islice+1]):
          islice+=1
      lapser=lr[islice]
      
      dz=z-zprev
      tcur=tprev+lapser*dz
      pcur=pprev-g/(R*0.5*(tcur+tprev))*pprev*dz
      if(z<zvals[-1]):
        ttup=ttup+(tcur,)
        ptup=ptup+(pcur,)

        zprev=z
        tprev=tcur
        pprev=pcur
  return ptup,ttup
