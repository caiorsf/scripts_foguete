import numpy as np
from scipy.optimize import fsolve, minimize,LinearConstraint
import math
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

#OXIDIZER

#VARIABLES
r_kso0=1e-3      
r_bxo0=0.3e-3    
r_ko0=2e-3       

#CONSTANTS
##Condições de projeto
###Parâmetros do liquido
delta_p = 4.4e5         #Queda de pressão 
rho_o = 1140            #Densidade
mu_o = 7.64e-6          #Viscosidade dinâmica
uk_o = 1.94e-6          #Viscosidade cinemática
mass_ro_target=0.1729   #Vazão mássica requerida
###Parâmetros geométricos
Xi_yo=0.3               #Razão entre o raio de entrada do canal tangencial e seu diâmetro
l_bxo=6*r_bxo0          #Comprimento do canal tangencial
n_o=6                   #Número de canais tangenciais

bnds=[(1e-3,10e-3),(0.5e-3,2.5e-3),(1e-3,10e-3)]

def func(u):
    r_kso = u[0]
    r_bxo = u[1]
    r_ko = u[2]
    
    # phi_i_o
    c=2**0.5*r_kso*(r_ko-r_bxo)
    f=c/(2*n_o*r_bxo**2)

    phi_i_o=float(fsolve(lambda x: f*x**1.5+x-1, 1.0))

    # lambda_o
    mu_io=phi_i_o*(phi_i_o/(2-phi_i_o))**0.5
    print('mu_io='+str(mu_io))
    f_kso=np.pi*r_kso**2
    mass_io=mu_io*f_kso*(2*rho_o*delta_p)**0.5
    W_o=mass_io / ( n_o * np.pi * r_bxo**2 * rho_o )
    Re_o=W_o*2*r_bxo*uk_o**-1
    lambda_o=0.3164*Re_o**(-0.25)

    # A_io
    A_io=r_kso*(r_ko-r_bxo) / (n_o*r_bxo**2)

    # phi_omega_o
    d=2*( ( (r_ko-r_bxo)*(r_ko-r_kso-r_bxo)*lambda_o / 2 ) +n_o*r_bxo**2 )
    e=c/d

    phi_omega_o=float(fsolve(lambda x: e*x**1.5+x-1,1.0))
    print('phi_omega_o='+str(phi_omega_o))

    # A_omega_o
    A_omega_o=((1-phi_omega_o)*2**0.5) / (phi_omega_o**1.5)

    # mass_ro
    K_0=A_io/A_omega_o
    C_o=r_ko-r_bxo/ r_kso
    f=np.pi*r_kso**2*(2*rho_o*delta_p)**0.5
    g=(phi_omega_o**-2 + (A_io**2*K_0**2) / (1-phi_omega_o) +Xi_yo*n_o*(A_io/C_o)**2)**0.5
    mass_ro=f/g 
    print("mass_ro = {}, u = {}".format(mass_ro, u))
    delta_mass = mass_ro - mass_ro_target
    return delta_mass**2

u0 = np.array([r_kso0, r_bxo0, r_ko0])

con=LinearConstraint([-1,0,1],[0],[np.inf])

result = minimize(func, u0, bounds=bnds, method="SLSQP", tol=1e-8,constraints=con)

result_new=[i for i in result.x]

print('r_kso={:.1E},r_bxo={:.1E},r_ko={:.1E},l_bxo={:.1E}'.format(result_new[0],result_new[1],result_new[2],6*result_new[1]))
