import numpy as np
from scipy.integrate import quad

#spacelike
p = 0.8

def rp_l(z):
    return (1 - np.exp(-z))**p

def b0(n):
    return 11 - 2/3*n

def K0(n):
    return b0(n) / (4 * np.pi)

def rho0_l(z, n):
    return K0(n) * (-np.pi) * rp_l(z)

def e0_l(z, n):
    def integrand_l(x, z, n):
        return rho0_l(x, n) / ((x + z) * (x + 1))
    
    return (1 - z) / np.pi * quad(integrand_l, 0, np.inf, args=(z, n))[0]

def epsilon0_l(t, L, n):
    return e0_l(t / L**2, n) - e0_l(0, n)

def b1(n):
    return 102 - 38/3*n

def K1(n):
    return b1(n) / (4 * np.pi * b0(n))

def F(x):
    if 0 < x <= 1:
        return -np.pi - np.arctan(np.pi / np.log(x))
    else:
        return -np.arctan(np.pi / np.log(x))

def rho1_l(z, n):
    return K1(n) * rp_l(z) * F(z)

def e1_l(z, n):
    def integrand_l(x, z, n):
        return rho1_l(x, n) / ((x + z) * (x + 1))
    
    return [(1 - z) / np.pi * quad(integrand_l, 0, np.inf, args=(z, n))[0],
            (1 - z) / np.pi * quad(integrand_l, 0, np.inf, args=(z, n))[1]]
    
def epsilon1_l(t, L, n):
    return e1_l(t / L**2, n)[0] - e1_l(0, n)[0]

#print(epsilon0_l(2,0.3,3))
#print(epsilon1_l(2,0.3,3)) 

def epsilon_l(t, L, n):
    return epsilon0_l(t, L, n) + epsilon1_l(t, L, n)

#print(epsilon_l(2, 0.3, 3))

##############################################################à

#timelike

def G_l(z, n):
    def integrand_l(x, z, n):
        return rho0_l(x, n) / ((x - z) * (x + 1))
    integral_left, _ = quad(integrand_l, 0, z - 1e-6, args=(z, n))
    integral_right, _ = quad(integrand_l, z + 1e-6, np.inf, args=(z, n))
    
    #print(integral_left)
    #print(integral_right)

    # Somma i risultati degli integrali delle due parti
    return (1 + z) / np.pi * (integral_left + integral_right)

    #integrale a valore principale intero con inf (no python)
    #return (1 + z) / np.pi * quad(integrand_l, 0, np.inf, args=(z, n), weight='cauchy', wvar=z)[0]
'''
def G_l(z,n):
    def integrand(x, z, n):
        return rho0_l(x, n) / ((x - z) * (x + 1))
    
    # Definisci gli estremi degli intervalli finiti che circondano la singolarità
    a = 0
    b = 1
    
    # Suddividi l'intervallo di integrazione intorno alla singolarità
    z1 = (a + z) / 2
    z2 = (b + z) / 2
    
    # Calcola l'integrale su ciascun sotto-intervallo
    integral_left, _ = quad(integrand, a, z1, args=(z, n))
    integral_center, _ = quad(integrand, z1, z2, args=(z, n))
    integral_right, _ = quad(integrand, z2, b, args=(z, n))
    integral_inf, _ = quad(integrand, b, np.inf, args=(z, n))
    
    print(integral_left)
    print(integral_center)
    print(integral_right)
    print(integral_inf)

    # Somma i risultati degli integrali dei sotto-intervalli
    return (1 + z) / np.pi * (integral_left + integral_center + integral_right+integral_inf)
'''
def real_l(z, n):
    return G_l(z, n) - e0_l(0, n)

def rhoa_l(z, n):
    return -rho0_l(z, n) / (real_l(z, n)**2 + rho0_l(z, n)**2)

def A_l(s, L, n):
    def integrand_l(sigma, L, n):
        return 1 / sigma * rhoa_l(sigma / L**2, n)
    
    return[1 / np.pi * quad(integrand_l, s, np.inf, args=(L, n))[0],1 / np.pi * quad(integrand_l, s, np.inf, args=(L, n))[1]]
    #return 1 / np.pi * quad(integrand_l, s, 100000*s, args=(L, n))[0]

#print(G_l(2,3))
#print(rhoa_l(100,3))
#print(A_l(0.001, 0.3, 3))
