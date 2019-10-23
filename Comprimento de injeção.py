#Oxidante
mass_ro = 0.1729343632791287
mu_io=0.5141532076547852
phi_omega_o=0.7020209970839623
delta_po=4.4e5
rho_o =1140
#Combustível
mass_rc = 0.0721982629534907
phi_c=0.1339960640658725
mu_c=0.035907220197069774
delta_pc=7.1e5
rho_c=800

km=mass_ro/mass_rc
tau=0.11e-3

a=( km*mu_c / ( (km+1)*phi_c ) ) * (delta_pc/rho_c)**0.5
b=( mu_io / ( ( km+1 )*phi_omega_o ) ) *  (delta_po/rho_o)**0.5

l_inj=(2**0.5)*tau*(a+b)

print('l_inj(m)={:.2E}'.format(l_inj))