import numpy as np
from scipy.optimize import fsolve, minimize
import math
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper

#FUEL

#VARIABLES
r_kc0=4.5e-3
r_ksc=r_kc0
r_bxc0=0.3e-3

#CONSTANTS
##Condições de projeto
###Parâmetros do liquido
delta_p=7.1e5           #Queda de pressão
rho_c=800               #Densidade 
u_kc=2e-6               #Viscosidade cinemática
mass_rc_target=0.0722   #Vazão mássica requerida
###Parâmetros geométricos
Xi_yc=1.1               #Razão entre o raio de entrada do canal tangencial e seu diâmetro
l_bxc=6*r_bxc0          #Comprimento do canal tangencial
n_c=6                   #Número de canais tangenciais

K_c=1

def func(u):
    r_kc=u[0]
    r_bxc=u[1]
    R_bxc=r_kc-r_bxc
    A_c=R_bxc*r_ksc / (n_c*r_bxc**2)
    c=2**0.5*r_ksc*(r_kc-r_bxc)
    f=c/(2*n_c*r_bxc**2)
    phi_c=float(fsolve(lambda x: f*x**1.5+x-1, 1.0))
    print('phi_c:'+str(phi_c))
    mu_c=phi_c*(phi_c/(2-phi_c))**0.5
    print('mu_c:'+str(mu_c))
    #r_asc=((1-phi_c)*r_ksc)**0.5
    f_ksc=np.pi*r_ksc**2
    mass_ic=mu_c*f_ksc*(2*rho_c*delta_p)**0.5
    W_c=mass_ic/(n_c*np.pi*r_bxc**2*rho_c)
    Re_c=W_c*2*r_bxc*u_kc**-1
    lambda_c=0.3164*Re_c**(-0.25)
    Xi_c=Xi_yc+lambda_c*(l_bxc/(2*r_bxc))
    C_c=R_bxc/r_ksc
    f=np.pi*r_kc**2*(2*rho_c*delta_p)**0.5
    g=(phi_c**-2 + (A_c**2*K_c**2) / (1-phi_c) +Xi_c*n_c*(A_c/C_c)**2)**0.5
    mass_rc=f/g 
    print("mass_rc = {}, u = {}".format(mass_rc, u))
    delta_mass = mass_rc - mass_rc_target
    return delta_mass**2


u0 = np.array([r_kc0,r_bxc0])

bnds=[(1e-3,10e-3),(0.5e-3,2.5e-3)]

result = minimize(func, u0, bounds=bnds, method="SLSQP", tol=1e-8)

result_new=[i for i in result.x]

print('r_kc={:.1E},r_ksc={:.1E},r_bxc={:.1E},l_bxc={:.1E}'.format(result_new[0],result_new[0],result_new[1],result_new[1]*6))


