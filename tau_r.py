import numpy as np
from math import pi,sqrt
import scipy.constants as sc 


def returnR(star):
    tau = star.tau
    tau_inf = tau[-1]
    l = len(tau) - 1

    for i in range(l, -1, -1):
        diff = tau_inf - tau[i]

        if diff > (2/3):
            break

    return i
