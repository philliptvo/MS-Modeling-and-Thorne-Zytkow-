import numpy as np
#from mpmath import *
import matplotlib.pyplot as plt
from tau_r import returnR


# Constants 
gamma = 5/3
X = 0.78
Y = 0.20
Z = 0.02
X_CNO = 0.03*X
mu = 1 / (2*X + 0.75*Y + 0.5*Z)

G = 6.67e-11
hbar = 1.05457e-34
m_e = 9.109e-31
m_p = 1.672e-27
k_B = 1.38e-23
sigma = 5.67e-8
c = 3e8
a = (4 * sigma) / c
R_sun = 7e8
M_sun = 1.989e30
L_sun = 3.828e26


# ALl opacity functions (sum and individual components)
def opacity(rho, T):
    K_es = 0.02 * (1 + X)                                           # m^2/kg
    K_ff = 1e24 * (Z + 0.0001) * (rho / 1e3)**(0.7) * T**(-3.5)     # m^2/kg
    K_h = 2.5e-32 * (Z / 0.02) * (rho / 1e3)**(0.5) * T**9          # m^2/kg

    return 1 / (1/K_h + 1/max(K_es,K_ff))

def kappa_es():
    return 0.02 * (1 + X)

def kappa_ff(rho, T):
    return 1e24 * (Z + 0.0001) * ((rho / 1e3) ** 0.7) * T ** (-3.5)

def kappa_H(rho, T):
    return 2.5e-32 * (Z / 0.02) * ((rho / 1e3) ** 0.5) * T ** 9

# Al lenergy generation functions (components and sum)
def e_pp(rho, T):
    return 1.07e-7 * (rho / 1e5) * X ** 2 * (T / 1e6) ** 4

def e_CNO(rho, T):
    return 8.24e-26 * (rho / 1e5) * X * (X_CNO) * (T / 1e6) ** (19.9)

def energy_gen(rho, T):
    return e_pp(rho, T) + e_CNO(rho, T)


# Pressure function
def pressure(rho, T):
    P_deg = ((3 * np.pi ** 2) ** (2 / 3) * hbar ** 2 / (5 * m_e))*(rho / m_p) ** (5 / 3)
    P_ideal = (rho * k_B * T) / (mu * m_p) 
    P_photon = (1 / 3) * a * T ** 4
    
    return[P_deg, P_ideal, P_photon]

# Differential Equations 
def drhodr(r, rho, T, L, M, P):
    a = -((G * M * rho) / r ** 2 + (dPdT(T, P, rho) * dTdr(r, T, rho, L, M, P)))
    b =  dPdrho(rho, P, T)
    
    return a/b

def dPdr(r, P, rho, M):
    return -(G * M * rho) / r ** 2

def dPdrho(rho, P, T):
    a = (3 * np.pi ** 2) ** (2 / 3) * hbar ** 2 
    b =  (3 * m_e * m_p)
    c = (rho / m_p) ** (2 / 3)
    d = ((k_B * T) / (mu * m_p))
    
    return a*c/b + d

def dPdT(T, P, rho):
    term1 = (rho * k_B) / (mu * m_p)
    term2 = (4 / 3) * a * T ** 3
    
    return term1 + term2

def dTdr(r, T, rho, L, M, P):
    rad = (3 * opacity(rho, T) * rho * L) / (16 * np.pi * a * c * T ** 3 * r ** 2)
    con = (1 - (1 / gamma)) * (G * T * M * rho) / (P * r ** 2)
    
    return -1*min(rad, con)
    
def dMdr(r, M, rho):
    return 4 * np.pi * r ** 2 * rho

def dLdr(r, L, rho, energyRate):
    return 4 * np.pi * r ** 2 * rho * energyRate

def dTaudr(r, tau, rho, T):
    return opacity(rho, T) * rho

# Function to track 
def deltaTau(r, tau, rho, T, L, M, P):
    return (opacity(rho, T) * rho ** 2) / abs(drhodr(r, rho, T, L, M, P))

# RK4 METHOD

def findY(fn, x, y, dx, *args):
    dy1 = dx * fn(x, y, *args)
    dy2 = dx * fn(x + 0.5 * dx, y + 0.5 * dy1, *args)
    dy3 = dx * fn(x + 0.5 * dx, y + 0.5 * dy2, *args)
    dy4 = dx * fn(x + dx, y + dy3, *args)
    
    y = y + (1 / 6) * (dy1 + 2 * dy2 + 2 * dy3 + dy4)
    
    return y

def bisection(lowLimit, highLimit, T_c, massLimit, radiusLimit, N, M_BH = None):
#    print("-------------------------")
#    print("Checking rho_c {} g/cm^3".format(lowLimit))
    
    trialStar = Star(lowLimit * 1000, T_c, M_BH)
    f_lower = trialStar.solve(1e-1, massLimit, radiusLimit, False)
    
#    print("-------------------------")
#    print("Checking rho_c of {} g/cm^3".format(highLimit))
    
    trialStar = Star(highLimit * 1000, T_c, M_BH)
    f_higher = trialStar.solve(1e-1, massLimit, radiusLimit, False)    
    
    #print(f_lower, f_higher)
    
    if f_lower * f_higher >= 0:
#        print("Bisection method fails.")
        return None
    
    a_n = lowLimit
    b_n = highLimit
    
    for n in range(1,N+1):
        m_n = (a_n + b_n) / 2
        
#        print("-------------------------")
#        print("Checking rho_c of {} g/cm^3".format(m_n))
        
        trialStar = Star(m_n * 1000, T_c, M_BH)
        f_midpoint = trialStar.solve(1e-1, massLimit, radiusLimit, False)
        
        #print(f_midpoint)
        
        if abs(f_midpoint) < 0.01:
            break
        
        if f_lower * f_midpoint < 0:
            a_n = a_n
            b_n = m_n
            
            f_higher = f_midpoint
        elif f_higher * f_midpoint < 0:
            a_n = m_n
            b_n = b_n
            
            f_lower = f_midpoint
        elif f_midpoint == 0:
#            print("Found exact solution.")
            break
        else:
#            print("Bisection method fails2.")
            return None
    
#    print(n)
#    print("-------------------------")
    trialStar = Star(m_n * 1000, T_c)
    trialStar.solve(1e-3, massLimit, radiusLimit, True)
    
    return trialStar

class Star:
    def __init__(self, rho_c, T_c, M_BH = None):
                
        self.r = np.empty(1)
        self.rho = np.empty(1)
        self.P = np.empty(1)
        self.T = np.empty(1)
        self.M = np.empty(1)
        self.L = np.empty(1)
        self.tau = np.empty(1)
        self.deltaTau = np.empty(1)
        self.kappa = np.empty(1)
        self.k_H = np.empty(1)
        self.k_es = np.empty(1)
        self.k_ff = np.empty(1)
        self.dP = np.empty(1)
        self.dT = np.empty(1)
        self.dL = np.empty(1)
        self.dL_pp = np.empty(1)
        self.dL_CNO = np.empty(1)
        self.P_deg = np.empty(1)
        self.P_ideal = np.empty(1)
        self.P_gamma = np.empty(1)
        
        self.rho[0] = rho_c
        self.T[0] = T_c
              
        if M_BH == None:   # This means it's not a BH
            self.r[0] = 0.00001
            self.M[0] = (4 * np.pi * self.r[0] ** 3 * self.rho[0]) / 3
            self.L[0] = (4 * np.pi * self.r[0] ** 3 * self.rho[0] * 
              energy_gen(self.rho[0], self.T[0])) / 3
        else:
            self.M[0] = M_BH
            self.L[0] = (1.3e31) * (self.M[0]/M_sun)
            self.r[0] = np.sqrt(self.L[0] * mu * m_p/(4 * np.pi * c * k_B * self.T[0] * self.rho[0]))
        
        P_deg, P_ideal, P_gamma = pressure(self.rho[0], self.T[0])
        
        self.P[0] = P_deg + P_ideal + P_gamma
        self.M[0] = (4 * np.pi * self.r[0] ** 3 * self.rho[0]) / 3
        self.L[0] = (4 * np.pi * self.r[0] ** 3 * self.rho[0] * 
              energy_gen(self.rho[0], self.T[0])) / 3
       
        self.deltaTau[0] = 1
        self.tau[0] = 1
                     
        self.kappa[0] = opacity(self.rho[0], self.T[0])
        self.k_H[0] = kappa_H(self.rho[0], self.T[0])
        self.k_es[0] = kappa_es()
        self.k_ff[0] = kappa_ff(self.rho[0], self.T[0])
        
        self.dP[0] = dPdr(self.r[0], self.P[0], self.rho[0], self.M[0])
        self.dT[0] = dTdr(self.r[0], self.T[0], self.rho[0], self.L[0], 
               self.M[0], self.P[0])
        self.dL[0] = dLdr(self.r[0], self.L[0], self.rho[0], 
               energy_gen(self.rho[0], self.T[0]))
        self.dL_pp[0] = dLdr(self.r[0], self.L[0], self.rho[0], 
                  e_pp(self.rho[0], self.T[0]))
        self.dL_CNO[0] = dLdr(self.r[0], self.L[0], self.rho[0], 
                  e_CNO(self.rho[0], self.T[0]))
        
        
        self.P_deg[0] = P_deg
        self.P_ideal[0] = P_ideal
        self.P_gamma[0] = P_gamma
        
        self.convectiveRegion = -1

        self.h = 1000000
        
    def solve(self, tauLimit, massLimit, radiusLimit, solveFull):
        import time
        
        i = 0
        
        while self.deltaTau[i] > tauLimit and self.M[i] < massLimit * M_sun and self.r[i] < radiusLimit * R_sun:
            P_next = sum(pressure(self.rho[i], self.T[i]))
            self.P = np.append(self.P, P_next)
            
            rho_next = findY(drhodr, self.r[i], self.rho[i], self.h, self.T[i], self.L[i], self.M[i], 
                             self.P[i])
            self.rho = np.append(self.rho, rho_next)
            
            T_next = findY(dTdr, self.r[i], self.T[i], self.h, self.rho[i], self.L[i], self.M[i], 
                           self.P[i])
            self.T = np.append(self.T, T_next)
            
            M_next = findY(dMdr, self.r[i], self.M[i], self.h, self.rho[i])
            self.M = np.append(self.M, M_next)
            
            L_next = findY(dLdr, self.r[i], self.L[i], self.h, self.rho[i], energy_gen(self.rho[i], 
                           self.T[i]))
            self.L = np.append(self.L, L_next)
            
            
            Tau_next = findY(dTaudr, self.r[i], self.tau[i], self.h, self.rho[i], self.T[i])
            self.tau = np.append(self.tau, Tau_next)
            
            deltaTau_next = deltaTau(self.r[i], self.deltaTau[i], self.rho[i], self.T[i], self.L[i], 
                                     self.M[i], self.P[i])
            
            self.deltaTau = np.append(self.deltaTau, deltaTau_next)
                       
            
            if solveFull:
                self.kappa = np.append(self.kappa, opacity(self.rho[i], self.T[i]))
                self.k_H = np.append(self.k_H, kappa_H(self.rho[i], self.T[i]))
                self.k_es = np.append(self.k_es, kappa_es())
                self.k_ff = np.append(self.k_ff, kappa_ff(self.rho[i], self.T[i]))
                
                self.dP = np.append(self.dP, dPdr(self.r[i], self.P[i], self.rho[i], self.M[i]))
                self.dT = np.append(self.dT, dTdr(self.r[i], self.T[i], self.rho[i], self.L[i], 
                                                  self.M[i], self.P[i]))
                
                self.dL = np.append(self.dL, dLdr(self.r[i], self.L[i], self.rho[i], 
                                                 energy_gen(self.rho[i], self.T[i])))
                
                self.dL_pp = np.append(self.dL_pp, dLdr(self.r[i], self.L[i], self.rho[i], 
                                                 e_pp(self.rho[i], self.T[i])))
                
                self.dL_CNO = np.append(self.dL_CNO, dLdr(self.r[i], self.L[i], self.rho[i], 
                                                 e_CNO(self.rho[i], self.T[i])))
                
                P_deg, P_ideal, P_gamma = pressure(self.rho[i], self.T[i])
                self.P_deg = np.append(self.P_deg, P_deg)
                self.P_ideal = np.append(self.P_ideal, P_ideal)
                self.P_gamma = np.append(self.P_gamma, P_gamma)
            
                if self.convectiveRegion == -1:
                    a = dTdr(self.r[i], self.T[i], self.rho[i], self.L[i],
                             self.M[i], self.P[i])
                    
                    convective = - (1 - (1 / gamma))
                    convective = convective * (G * self.T[i] * self.M[i] * self.rho[i])
                    convective = convective / (self.P[i] * self.r[i] ** 2)
                
                    if a == convective:
                        self.convectiveRegion = i
            
            self.r = np.append(self.r, self.r[i] + self.h)

        
            if self.rho[i] - self.rho[i + 1] > self.rho[i + 1]:
                self.h = 0.01 * self.h

            i = i + 1
            
        self.kappa = self.kappa * 10
        self.k_H = self.k_H * 10
        self.k_es = self.k_es * 10
        self.k_ff = self.k_ff * 10
        
        if self.M[i] >= massLimit * M_sun or self.r[i] >= radiusLimit * R_sun:
            pass
#            print("Convergence failed! This process took {}s.".format(end - start))
        else:
            pass
#            print("Convergence successful! This process took {}s.".format(end - start))
        
        r_star = returnR(self)
        
        # truncating the arrays for stars at the radius
        
        self.r = self.r[:r_star]
        self.rho = self.rho[:r_star]
        self.P = self.P[:r_star]
        self.T = self.T[:r_star]
        self.M = self.M[:r_star]
        self.L = self.L[:r_star]
        self.tau = self.tau[:r_star]
        self.deltaTau = self.deltaTau[:r_star]
        self.kappa = self.kappa[:r_star]
        self.k_H = self.k_H[:r_star]
        self.k_es = self.k_es[:r_star]
        self.k_ff = self.k_ff[:r_star]
        self.dP = self.dP[:r_star]
        self.dT = self.dT[:r_star]
        self.dL = self.dL[:r_star]
        self.dL_pp = self.dL_pp[:r_star]
        self.dL_CNO = self.dL_CNO[:r_star] 
        self.P_deg = self.P_deg[:r_star]
        self.P_ideal = self.P_ideal[:r_star]
        self.P_gamma = self.P_gamma[:r_star]
        
        return self.rho_c_fn(self.L[-1], self.T[-1], self.r[-1])
        
    def rho_c_fn(self, L_surf, T_surf, R_star):
        x = L_surf - (4 * np.pi * sigma * R_star ** 2 * T_surf ** 4)
        x = x / np.sqrt(4 * np.pi * sigma * R_star ** 2 * T_surf ** 4 * L_surf)
    
        return x





