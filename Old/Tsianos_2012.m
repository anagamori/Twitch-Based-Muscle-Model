
close all
clear all
clc

V = 0;
f = 1;
Lce = 1;

e_1 = -76.6;
e_2 = -792;
e_3 = 124;
e_4 = 0.72;

a = 1/3;
R = 1.5;

a_f = 0.56;
n_f_0 = 2.1;
n_f_1 = 5;
n_f = n_f_0+n_f_1*(1/Lce-1);

omega = 1.12;
beta = 2.3;
rho = 1.62;

E_initial_f_tet_Lce_0_V_0 =( e_1*V^2+e_2*V+e_3)/(e_4-V);

E_a_f_tet  = a*E_initial_f_tet_Lce_0_V_0;
E_xb_f_tet_Lce_0_V_0  = E_initial_f_tet_Lce_0_V_0 - E_a_f_tet;

Af = 1-exp(-(f/(a_f*n_f))^n_f);
FL = exp(-((Lce^beta-1)/omega)^rho);

E_xb_f_Lce_Vce = Af*FL*E_xb_f_tet_Lce_0_V_0;
E_a_f = a/(1-a)*E_xb_f_Lce_Vce;

E_recovery = (E_xb_f_Lce_Vce+E_a_f )*R;
