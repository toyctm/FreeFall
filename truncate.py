#!/bin/env python

import math
import numpy as np

def truncate_plus(x,ndigits):
    # Truncate to ndigits significant numbers
    nmag = int(math.floor(np.log10(x)))
    x1 = x/10**nmag    # bw 1.00000 and 9.999999
    return round(x1,ndigits-1)*10**nmag

def truncate(x,ndigits):
    if(x>0):
        return truncate_plus(x,ndigits)
    elif(x<0):
        return -truncate_plus(-x,ndigits)
    else:
        return 0.


