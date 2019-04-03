from time import time

import numpy as np
import matplotlib.pyplot as plt

from starClass import bisection
import pickle

massLimit = 100
radiusLimit = 25
X = 0.78
Y = 0.2
Z = 0.02
gamma = 5 / 3

rho_c_min = 0.3
rho_c_max = 500
N = 100

T_c = 8.23e6
L_sun = 3.828e26
R_sun = 6.9951e8

# Range of T_c for Main Sequence
rng_T = np.linspace(3e6,30e6,50)

M_bh = 1e-6

def solveStar(T):
    
    print("********** Creating star with Tc = {} **************".format(T))
    star_optimized = bisection(rho_c_min, rho_c_max, T, massLimit, radiusLimit, N, M_bh)
    if star_optimized is not None:
        pickle.dump( star_optimized, open( "stars/{}.p".format(T), "wb" ) )
    
    print("********** Finished star with Tc = {} **************".format(T))

def solveAllStars(T_lst):
    from multiprocessing import Pool
    pool = Pool(processes=4)
    
    start = time()
    
    pool.map(solveStar, T_lst)
    
    print("Time elapsed is {} s".format(time() - start))

def main_sequence(path="MS_stars", fig=None, axes=None):
    import os
    star_lst = os.listdir(path)
    name = os.path.basename(path)
    
    BV = np.array([])
    L = np.array([])
    stars = np.array([])

    for s in star_lst:
        tmp_star = pickle.load( open("{}/{}".format(path, s), "rb") )
        BV = np.append(BV, tmp_star.T[-1])
        L = np.append(L, tmp_star.L[-1])
        stars = np.append(stars, tmp_star)

    if axes is None:
        fig, axes = plt.subplots()
    
    axes.plot(np.log10(BV), np.log10(L/L_sun), '.', label=name)
    
    axes.set_xlim(4.5,3.0)
    axes.set_xlabel(r'T (K)')
    axes.set_ylabel(r'$\frac{L}{L_\odot}$')
    axes.grid()
    if axes is None:
        fig.savefig('{}'.format(os.path.basename(path)))
    return axes

def main_sequence_TZ(wMS=False, overplot=False):
    template = "stars(1e-{})"
    seqs = [2,3,4,5,9]
    
    fig, axes = plt.subplots()
    
    if wMS:
        main_sequence(fig=fig, axes=axes)
    
    for s in seqs:
        if overplot:
            main_sequence(template.format(s), fig, axes)
        else:
            main_sequence(template.format(s))
    
    axes.legend()
    if overplot:
        fig.savefig('Overplot.png')
    plt.show()

def load_star(path, Tc):
    return pickle.load( open("{}/{}.p".format(path, Tc), "rb"))

def plotAll(star):
    plot(star)
    plotPressure(star)
    plotOpacity(star)
    plotdLdr(star)
    plotdlogPdlogT(star)

    plt.show()
    print("Tc is {} K".format(star.T[0]))

def plot(star):
    r_surf = star.r[-1]
    rho0 = star.rho[0]
    T0 = star.T[0]
    L_surf = star.L[-1]
    M_enc = star.M[-1]

    fig, axes = plt.subplots()
    axes.plot(star.r/r_surf, star.rho/rho0,label=r'$\frac{\rho}{\rho_c}$ vs $\frac{R}{R_*}$')
    axes.plot(star.r/r_surf, star.T/T0,label=r'$\frac{T}{T_c}$ vs $\frac{R}{R_*}$')
    axes.plot(star.r/r_surf, star.L/L_surf,label=r'$\frac{L}{L_*}$ vs $\frac{R}{R_*}$')
    axes.plot(star.r/r_surf, star.M/M_enc,label=r'$\frac{M}{M_*}$ vs $\frac{R}{R_*}$')
    axes.set_xlabel(r'$\frac{R}{R_*}$')
    axes.axvspan(star.r[star.convectiveRegion]/r_surf, 1, facecolor='grey', alpha=0.5)
    axes.legend()

def plotPressure(star):
    r_surf = star.r[-1]
    P_max = star.P[0]

    fig, axes = plt.subplots()
    axes.plot(star.r/r_surf, star.P_gamma/P_max, label=r'$P_{\gamma}$')
    axes.plot(star.r/r_surf, star.P_deg/P_max, label=r'$P_{degenerate}$')
    axes.plot(star.r/r_surf, star.P_ideal/P_max, label=r'$P_{ideal}$')
    axes.plot(star.r/r_surf, star.P/P_max, label=r'$P_{total}$')
    axes.set_xlabel(r'$\frac{R}{R_*}$')
    axes.axvspan(star.r[star.convectiveRegion]/r_surf, 1, facecolor='grey', alpha=0.5)
    axes.legend()

def plotOpacity(star):
    r_surf = star.r[-1]
    
    fig, axes = plt.subplots()
    axes.plot(star.r/r_surf, np.log10(star.kappa), label=r'$\kappa$')
    axes.plot(star.r/r_surf, np.log10(star.k_es), label=r'$\kappa_{es}$')
    axes.plot(star.r/r_surf, np.log10(star.k_ff), label=r'$\kappa_{ff}$')
    axes.plot(star.r/r_surf, np.log10(star.k_H), label=r'$\kappa_{H^-}$')
    axes.set_xlabel(r'$\frac{R}{R_*}$')
    axes.set_ylim(-2,10)
    axes.axvspan(star.r[star.convectiveRegion]/r_surf, 1, facecolor='grey', alpha=0.5)
    axes.legend()

def plotdLdr(star):
    r_surf = star.r[-1]
    
    fig, axes = plt.subplots()
    axes.plot(star.r/r_surf, star.dL, 'k.', label=r'$dL/dr$')
    axes.plot(star.r/r_surf, star.dL_pp, 'r-', label=r'$dL_{pp}/dr$')
    axes.plot(star.r/r_surf, star.dL_CNO, label=r'$dL_{cno}/dr$')
    axes.set_xlabel(r'$\frac{R}{R_*}$')
    axes.axvspan(star.r[star.convectiveRegion]/r_surf, 1, facecolor='grey', alpha=0.5)
    axes.legend()

def plotdlogPdlogT(star):
    r_surf = star.r[-1]
    
    fig, axes = plt.subplots()
    axes.plot(star.r/r_surf, (star.dP / star.dT) * (star.dP / star.dT), 'k.', label=r'$dlogP/dlogT$')
    axes.set_xlabel(r'$\frac{R}{R_*}$')
    axes.axvspan(star.r[star.convectiveRegion]/r_surf, 1, facecolor='grey', alpha=0.5)
    axes.legend()

#solveAllStars(rng_T)
#S = main_sequence()
