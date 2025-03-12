Ks=11.7; % Relative permittivity of Silicon (unitless)
e_0=8.854e-14; % in F/cm
q=1.6e-19; % Charge of an electron in C
k=8.617e-5; % in J/K
T=300; % Temp in Kelvin
ni=1e10; % Intrinic carrier concentration in Si at 300k (cm^-3)
Eg=1.12; % Gap energy for Si in eV

N_A=5e14; % P-side acceptor doping in cm^-3
N_D=1e15; % N-side donor doping in cm^-3
V_A=0; % Bias voltage

Vbi=k*T*log(N_A*N_D/ni^2); % Built in electrostatic potential in V

x_n=sqrt(2*Ks*e_0/q*N_A/(N_D*(N_A+N_D))*(Vbi-V_A)); % n-side depletion region width in cm
x_p=sqrt(2*Ks*e_0/q*N_D/(N_A*(N_A+N_D))*(Vbi-V_A)); % p-side depletion region width in cm

% all in cm
maximum=2e-4;
minimum=-maximum;
dx=1e-7;
x=minimum:dx:maximum;

a=zeros(1, int32((maximum-x_p)/dx));
b=ones(1, int32(x_p/dx))*-N_A*q;
c=ones(1, int32(x_n/dx))*N_D*q;
d=zeros(1, int32((maximum-x_n)/dx+1));
rho=[a b c d]; % charge density in cm^-3

subplot(4,3,1);
plot(x, rho);
title('Plot D')
xlabel('X');
ylabel('Charge Density(C/cm^3)');

E=cumtrapz(rho)/(Ks*e_0)*dx; % multiply by dx since x is spaced dx cm apart
subplot(4,3,2);
plot(x,E);
title('Plot E');
xlabel('X');
ylabel('Electric Field(V/cm)');

V=-cumtrapz(E)*dx; % multiply by dx since x is spaced dx cm apart
subplot(4,3,3);
plot(x,V);
title('Plot F');
xlabel('X');
ylabel('Electrostatic Potential(V)');

E_c=-V+Eg+Vbi; % set bottom of Valence band on n side to 0 energy
E_v=E_c-Eg;
E_i=(E_c+E_v)/2;
% E_f=E_i(int32(length(x)/2))*ones(1, length(x));Ef=Ei at the
% junction-False
E_f=ones(1, length(E_i))*(E_i(length(E_i))+k*T*log(N_D/ni));
subplot(4,3,4);
hold on;
plot(x,E_c, 'DisplayName', 'E_c');
plot(x,E_v, 'DisplayName', 'E_v');
plot(x,E_i, 'DisplayName', 'E_i');
plot(x,E_f, 'DisplayName', 'E_f');
legend('Ec', 'Ev', 'Ei', 'Ef');
title('Plot G');
xlabel('X');
ylabel('E(eV)');
hold off;

p=ni*exp((E_i-E_f)/(k*T));
n=ni^2*ones(1,length(p))./p;
% n=ni*exp((E_f-E_i)/(k*T));
subplot(4,3,5);
hold on;
plot(x,p);
plot(x,n);
hold off;
legend('p', 'n');
title('Plot H');
xlabel('X');
ylabel('Carrier Density(cm^-3)');

%recalculate charge density using carrier densities for a better
%approximation
rho_real=q*(p - n + [-ones(1,int32(length(x)/2))*N_A, ones(1,int32(length(x)/2)-1)*N_D]);
subplot(4, 3, 6);
plot(x, rho_real);
title('Plot I');
xlabel('X');
ylabel('Charge Density(C/cm^3');

% Must recalculate depletion region bounds since bias applies
V_A=.1;
x_n=sqrt(2*Ks*e_0/q*N_A/(N_D*(N_A+N_D))*(Vbi-V_A));
x_p=sqrt(2*Ks*e_0/q*N_D/(N_A*(N_A+N_D))*(Vbi-V_A));

L=.25e-4; % Minority carrier diffusion length in cm
n_p=n(1) + n(1)*(exp(V_A/(k*T))-1)*exp(((x+ones(1, length(x))*(x_p))/L)); % Minority electron concentration in p-side in cm^-3
p_n=p(length(p)) + p(length(p))*(exp(V_A/(k*T))-1)*exp((-(x-ones(1, length(x))*(x_n))/L)); % Minority hole concentration in n-side in cm^-3
bound = -ones(1, length(x))*x_p > x;
n_p=n_p.*bound;
bound = ones(1, length(x))*x_n < x;
p_n=p_n.*bound;

subplot(4,3,7);
hold on;
plot(x, zeros(1, length(x)));
plot(x,n_p);
plot(x,p_n);
hold off;
title('Plot K');
xlabel('X');
ylabel('Minority Carriers(cm^-3)');
legend('V=0V', 'V=0.1V:n_p', 'V=0.1V:p_n');

D=30;
J_ndiff=q*D*gradient(n_p)/dx; % Electron diffusion current on p-side in A/cm^3
J_pdiff=-q*D*gradient(p_n)/dx; % Hole diffusion current on n-side in A/cm^3
% Remove unrealistic points from sharp cuttoff on minotiry carrier concentrations
bound = -ones(1, length(x))*(x_p+dx) > x;
J_ndiff=J_ndiff.*bound;
bound=ones(1, length(x))*(x_n+dx) < x;
J_pdiff=J_pdiff.*bound;
J_diff=J_ndiff+J_pdiff;

subplot(4,3,9);
plot(x,J_diff);
title('Plot L');
xlabel('X');
ylabel('Diffusion Current Density(A/cm^3) V=.1V');

J_drift=ones(1, length(x))*max(J_diff)-J_diff;
subplot(4,3,11);
plot(x,J_drift);
title('Plot M');
xlabel('X');
ylabel('Drift Current Density(A/cm^3) V=.1V');

maximum=5e-4;
minimum=-maximum;
dx=1e-7;
x=minimum:dx:maximum; % Need a wider plot to account for much larger depletion region at negative biasing

% Repeat previous steps for a large negative bias

V_A=-5;
x_n=sqrt(2*Ks*e_0/q*N_A/(N_D*(N_A+N_D))*(Vbi-V_A));
x_p=sqrt(2*Ks*e_0/q*N_D/(N_A*(N_A+N_D))*(Vbi-V_A));

n_p=n(1) + n(1)*(exp(V_A/(k*T))-1)*exp(((x+ones(1, length(x))*(x_p))/L));
p_n=p(length(p)) + p(length(p))*(exp(V_A/(k*T))-1)*exp((-(x-ones(1, length(x))*(x_n))/L));
bound = -ones(1, length(x))*x_p > x;
n_p=n_p.*bound;
bound = ones(1, length(x))*x_n < x;
p_n=p_n.*bound;

subplot(4,3,8);
hold on;
plot(x,n_p);
plot(x,p_n);
hold off;
title('Plot K Part 2');
xlabel('X');
ylabel('Minority Carriers(cm^-3)');
legend('V=-5V:n_p', 'V=-5V:p_n');

J_ndiff=q*D*gradient(n_p)/dx; % Correct for x step size
J_pdiff=-q*D*gradient(p_n)/dx;
bound = -ones(1, length(x))*(x_p+dx) > x;
J_ndiff=J_ndiff.*bound;
bound=ones(1, length(x))*(x_n+dx) < x;
J_pdiff=J_pdiff.*bound;
J_diff=J_ndiff+J_pdiff;

subplot(4,3,10);
plot(x,J_diff);
title('Plot L Part 2');
xlabel('X');
ylabel('Diffusion Current Density(A/cm^3) V=-5V');