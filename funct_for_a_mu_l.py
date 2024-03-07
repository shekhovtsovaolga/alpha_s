from decimal import *
from numpy import sqrt, sin, cos, pi, radians, arctan, arctan2, log, sinh, exp,arccosh, inf, arccos
from scipy import integrate
from scipy.integrate import quad, romberg, fixed_quad
from alpha_s_lorenzo import*

###################################### Quark parameter values #######################################
# Quark masses
def quark_mass(m_u=0.26, m_d=0.26, m_s=0.45, m_c=1.35, m_b=4.4, m_t=174.):
    return {'m_u': m_u, 'm_d': m_d, 'm_s': m_s, 'm_c': m_c,  'm_b': m_b, 'm_t': m_t}

# Quark charges
def quark_charge():
    return {'q_u': 2./3.,  'q_c': 2./3.,  'q_t': 2./3,
            'q_d': -1./3., 'q_s': -1./3,  'q_b': -1./3.}

# Quark parameters
def quark_parameter():
    return {'u': {'m_f': quark_mass()['m_u'], 'q_f': quark_charge()['q_u']},
            'd': {'m_f': quark_mass()['m_d'], 'q_f': quark_charge()['q_d']},
            's': {'m_f': quark_mass()['m_s'], 'q_f': quark_charge()['q_s']},
            'c': {'m_f': quark_mass()['m_c'], 'q_f': quark_charge()['q_c']},
            'b': {'m_f': quark_mass()['m_b'], 'q_f': quark_charge()['q_b']},
            't': {'m_f': quark_mass()['m_t'], 'q_f': quark_charge()['q_t']}}

###################################### REPRESENTATION FOR FUNCTIONS K(s) ##########################
# Analytical representation of function K(s)
def funct_k_analyt(s, m_mu=0.1057, m_pi=0.139):
    if s <= 4.*m_pi**2:
        return 0
    else:
        beta = sqrt(1. - 4.*m_mu**2/s)
        x = (1. - beta)/(beta + 1.)
        return ( x**2/2.*(2. - x**2) + (1. + x**2)*(1. + x)**2/x**2*(log(1. + x) - x + x**2/2.)
                 + (1. + x)/(1. - x)*x**2*log(x) )

# Decimal type for analytical representation of function K(s)
def decimal_funct_k_analyt(s, m_mu=0.1057, m_pi=0.139):
    if s <= 4.*m_pi**2:
        return 0
    else:
        beta = (Decimal(1.) - Decimal(4.*m_mu**2/s)).sqrt()
        x = (Decimal(1.) - beta)/(beta + Decimal(1.))
        return ( x**2/Decimal(2.)*(Decimal(2.) - x**2) +
                (Decimal(1.) + x**2)*(Decimal(1.) + x)**2/x**2*((Decimal(1.) + x).ln() - x + x**2/Decimal(2.))
                 + (Decimal(1.) + x)/(Decimal(1.) - x)*x**2*x.ln() )
    
# Functions K(s) from p.1 , eq. (1).
# The two parametrizations for K(s) are tested:
# Function K(s) written as integral_0^1
def funct_k(s, m_mu=0.105):
    expr = lambda x: x**2*(1. - x)/(x**2 + (1. - x)*s/m_mu**2)
    # Integration of expr from 0 to 1 
    return quad(expr, 0, 1)

# Function hat_k is (0.69; 1) for the whole region, it is used in a_mu calculation
def funct_hat_k_analyt(s, s_limit_1=40000, m_mu=0.1057, m_pi=0.139):
    # The threshold s_limit_1=40000 is gotten emperically,
    # see studies in a_mu_alpha_0
    if s <= s_limit_1:
        # Integration of expr from 0 to 1
        return funct_k_analyt(s, m_mu, m_pi)*3.*s/m_mu**2
    else:
        return 1.

def decimal_funct_hat_k_analyt(s, s_limit_1=40000, m_mu=0.1057, m_pi=0.139):
    # The threshold s_limit_1=40000 is gotten emperically,
    # see studies in a_mu_alpha_0
    if s <= s_limit_1:
        # Integration of expr from 0 to 1
        return decimal_funct_k_analyt(s, m_mu, m_pi)*Decimal(3.*s/m_mu**2)
    else:
        return 1.

def funct_hat_k(s, m_mu=0.1057):
    return funct_k(s)[0]*3.*s/m_mu**2 

###################################### FUNCTIONS FROM PAGE 2 ##########################################
#
def funct_v(m_f, s):
    if s > 4.*m_f**2:
        return sqrt(1. - 4.*m_f**2/s)
    else:
        return 0

def funct_chi(m_f, s):
    if s > 4.*m_f**2:
        return arccosh(sqrt(s)/(2.*m_f))
    else:
        return 0
    
def funct_t(m_f, s):
    if s > 4.*m_f**2:
        return (3. - funct_v(m_f, s)**2)/2.*funct_v(m_f, s)
    else:
        return 0

def funct_g(m_f, s):
    if s > 4.*m_f**2:
        return 4./3.*pi*( pi/(2.*funct_v(m_f, s)) - 
                     (3. + funct_v(m_f, s))/4.*(pi/2. - 3./(4.*pi)) ) 
    else:
        return 0
                     
def funct_x(m_f, s, a=0.2):
    if s > 4.*m_f**2:
        return pi*a/sinh(funct_chi(m_f, s))
    else:
        return 0
    
def funct_s(m_f, s, a=0.2):
    if s > 4.*m_f**2:
        return funct_x(m_f, s, a)/(1. - exp(-funct_x(m_f, s, a)))
    else:
        return 0

###################################### DEPENDENCE OF ACTIVE QUARKS ON  s ##########################
#
def n_f(s, fixed_3=False, nf_low_3=False):
    m_u = quark_mass()['m_u']
    m_s = quark_mass()['m_s']
    m_c = quark_mass()['m_c']
    m_b = quark_mass()['m_b']
    m_t = quark_mass()['m_t']

    if fixed_3:
        return 3.

    if s < 4.*m_u**2:
        if nf_low_3:
            return 3.
        else:
            return 0.
    if 4.*m_u**2 <= s < 4.*m_s**2:
        if nf_low_3:
            return 3.
        else:
            return 2.

    if 4.*m_s**2 <= s < 4.*m_c**2:
        return 3.
    if 4.*m_c**2 <= s < 4.*m_b**2:
        return 4.
    if 4.*m_b**2 <= s < 4.*m_t**2:
        return 5.
    if s>= 4.*m_t**2:
        return 6.
    
######################################  PARAMETRIZATION FOR alpha_s ###########################################
#
def beta_0(s, fixed_3=False, nf_low_3=False):
    return 11. - 2./3.*n_f(s, fixed_3, nf_low_3)

# arXiv:1604.08082, Eq. (4.54)
def funct_alpha_0(s, lambda1=0.3, fixed_3=False, nf_low_3=False):
    f = lambda1**2/(lambda1**2 - s) + 1./log(s/lambda1**2)
    return 4.*pi/beta_0(s, fixed_3, nf_low_3)*f

# arXiv:0611229, Eq. (2.14)
def funct_alpha_1(s, lambda1=0.3, fixed_3=False, nf_low_3=False):
    L1 = log(s / lambda1 ** 2)
    return 4./(beta_0(s, fixed_3, nf_low_3))*arccos(L1/sqrt(L1**2 + pi **2))

def funct_alpha(s, alpha_order=0,  lambda1=0.3, fixed_3=False, nf_low_3=False):
    if alpha_order==0:
        return funct_alpha_0(s, lambda1, fixed_3, nf_low_3)
    elif alpha_order==1:
        return funct_alpha_1(s, lambda1, fixed_3, nf_low_3)
    else:
        raise ValueError('Wrong value alpha_order')

######################################  FUNCTION R_f(s): PAGE 1  ###########################################
#
# Function R_f(s), last Equation from p.1
def funct_r_f(m_f, q_f, s, alpha_order=0, lambda1=0.3,fixed_3=False, nf_low_3=False, a=0.2):
    return 3.*q_f**2*funct_t(m_f, s)*(funct_s(m_f, s, a) +
                             A_l(s, lambda1,n_f(s,fixed_3,nf_low_3)) * funct_g(m_f, s) -
                             funct_x(m_f, s, a)/2. )
    
def funct_r(s, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2):
    m_u = quark_mass()['m_u']
    m_d = quark_mass()['m_d']
    m_s = quark_mass()['m_s']
    m_c = quark_mass()['m_c']
    m_b = quark_mass()['m_b']
    m_t = quark_mass()['m_t']

    q_u = quark_charge()['q_u']
    q_d = quark_charge()['q_d']
    q_s = quark_charge()['q_s']
    q_c = quark_charge()['q_c']
    q_b = quark_charge()['q_b']
    q_t = quark_charge()['q_t']

    if s < 4.*m_u**2:
        return 0
    
    if 4.*m_u**2 <= s < 4.*m_s**2:
        return (funct_r_f(m_u, q_u, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_d, q_d, s, alpha_order, lambda1, fixed_3, nf_low_3, a))
    
    if 4.*m_s**2 <= s < 4.*m_c**2:
        return (funct_r_f(m_u, q_u, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_d, q_d, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_s, q_s, s, alpha_order, lambda1, fixed_3, nf_low_3, a))
    
    if 4.*m_c**2 <= s < 4.*m_b**2:
        return (funct_r_f(m_u, q_u, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_d, q_d, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_s, q_s, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_c, q_c, s, alpha_order, lambda1, fixed_3, nf_low_3, a))
        
    if 4.*m_b**2 <= s < 4.*m_t**2:
        return (funct_r_f(m_u, q_u, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_d, q_d, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_s, q_s, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_c, q_c, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_b, q_b, s, alpha_order, lambda1, fixed_3, nf_low_3, a))
    
    if s >= 4.*m_t**2:
        return (funct_r_f(m_u, q_u, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_d, q_d, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_s, q_s, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_c, q_c, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_b, q_b, s, alpha_order, lambda1, fixed_3, nf_low_3, a) +
                funct_r_f(m_t, q_t, s, alpha_order, lambda1, fixed_3, nf_low_3, a))

######################################  a_mu_had(s)  ###########################################
# TODO: to check formulae for amu_s
# method: 'quad' or 'romberg'
# f: type of a quark u, d, s, c, t, b
# s1, s2: low and high limit of integration
# s_limit_1: value of s where K_hat(s > s_limit_1) = 1
# alpha_order: order of alpha_s value calculation 0/1
# lambda1: value of alpha  ->  log(s/Lambda1**2)
# fixed_3: True for active quark number is 3
# nf_low_3: for low energy s<= 4*m_s**2 the number of active quarks is 3
# a:    default value a=0.2, is used in functions funct_s(m_f, s, a) and funct_x(m_f, s, a)
# integral: True for integral representation of K_hat

def d_amu_s(s,  s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2,
            integral=True):
    m_mu = 0.1057
    alpha = 1./137.
    m_pi = 0.139
    if integral:
        return ((alpha*m_mu/3./pi)**2*
                funct_r(s, alpha_order, lambda1, fixed_3, nf_low_3, a) * funct_hat_k(s, m_mu)
                /s**2)
    else:
        return ((alpha*m_mu/3./pi)**2*
                funct_r(s, alpha_order, lambda1, fixed_3, nf_low_3, a)*funct_hat_k_analyt(s, s_limit_1, m_mu, m_pi)
                /s**2)

def amu_s1_s2(method, s1, s2, s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2,
              integral=True):
    expr = lambda s: d_amu_s(s,  s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral)
    if method == 'quad':
        int1 = quad(expr, s1, s2, limit=500, epsabs=1.5e-12, epsrel=1.5e-12)
    elif method == 'romberg':
        int1 = romberg(expr, s1, s2, tol=1.48e-12, rtol=1.48e-12, divmax=50)
    else:
        raise ValueError('Wrong integral method')

    return int1

def d_amu_s_quark_f(f, s, s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2,
                    integral=True):
    m_f = quark_parameter()[f]['m_f']
    q_f = quark_parameter()[f]['q_f']
    m_mu = 0.1057
    alpha = 1. / 137.
    m_pi = 0.139
    if integral:
        return ((alpha*m_mu/3./pi)**2*
                funct_r_f(m_f, q_f, s, alpha_order, lambda1, fixed_3, nf_low_3, a) * funct_hat_k(s, m_mu)
                /s**2)
    else:
        return ((alpha*m_mu/3./pi)**2*
                funct_r_f(m_f, q_f, s, alpha_order, lambda1, fixed_3, nf_low_3, a)*
                funct_hat_k_analyt(s, s_limit_1, m_mu, m_pi)
                /s**2)


def amu_s1_s2_quark_f(method, f, s1, s2, s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False,
                      a=0.2, integral=True):
    expr = lambda s: d_amu_s_quark_f(f, s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral)
    if method == 'quad':
        int1 = quad(expr, s1, s2, limit=500, epsabs=1.5e-12, epsrel=1.5e-12)
    elif method == 'romberg':
        int1 = romberg(expr, s1, s2, tol=1.48e-12, rtol=1.48e-12, divmax=50)
    else:
        raise ValueError('Wrong integral method')

    return int1


# Another way to calculate the total a_mu contribution summed over quarks:
def d_amu_s_quark_total(s, s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2,
                        integral=True):
    return (d_amu_s_quark_f('u', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) +
            d_amu_s_quark_f('d', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) +
            d_amu_s_quark_f('s', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) +
            d_amu_s_quark_f('c', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) +
            d_amu_s_quark_f('b', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) +
            d_amu_s_quark_f('t', s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral) )

def amu_s1_s2_quark_total(method, s1, s2, s_limit_1=40000, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False,
                          a=0.2, integral=True):
    expr = lambda s: d_amu_s_quark_total(s, s_limit_1, alpha_order, lambda1, fixed_3, nf_low_3, a, integral)
    if method == 'quad':
        int1 = quad(expr, s1, s2, limit=500, epsabs=1.5e-12, epsrel=1.5e-12)
    elif method == 'romberg':
        int1 = romberg(expr, s1, s2, tol=1.48e-12, rtol=1.48e-12, divmax=50)
    else:
        raise ValueError('Wrong integral method')

    return int1

            
# TODO: # TO TEST HERE DIVISION ON REGION NUMBER AND OTHER METHOD OF INTEGRATION  AND FIND CONTRIBUTION FROM  s and c quarks
# Estimation of high energy integral over t = 1/s  -> dt = ds/s**2, K(s) is fixed to 1 -> to be discussed as
#  K_hat_integrand is 0.98 rather than 1. for high s
def integrand_high_energy(t, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2):
    m_mu = 0.1057
    alpha = 1. / 137.
    m_pi = 0.139
    return funct_r(1./t, alpha_order, lambda1, fixed_3, nf_low_3, a)*(alpha*m_mu/3./pi)**2

def amu_high_energy_estimation(method, t_max, alpha_order=0, lambda1=0.3, fixed_3=False, nf_low_3=False, a=0.2, t_min=0):
    expr = lambda t: integrand_high_energy(t, alpha_order, lambda1, fixed_3, nf_low_3, a)
    if method == 'quad':
        int1 = quad(expr, t_min, t_max, limit=500, epsabs=1.5e-12, epsrel=1.5e-12)
    elif method == 'romberg':
        int1 = romberg(expr, t_min, t_max, tol=1.48e-12, rtol=1.48e-12, divmax=50)
    else:
        raise ValueError('Wrong integral method')

    return int1











