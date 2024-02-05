clear all
close all
clc

set(0,'defaultTextInterpreter','latex');
set(0,'DefaultAxesFontSize', 14);

load('tracesbar1.mat');

dt = 1e-3; % Sampling rate = 1000 Hz
w = 0.03; % Web thickness [m]
A_t_vec = pi/4 * (1e-3*[28.80 25.25 21.81]').^2; % Respectively [High Medium Low] pressure [m^2]
C_star_id = [1548 1554 1560]; % [m/s]
rho_p = 1 / (.68/1950 + .18/2700 + .14/920);
D_e = 0.160;
D_i = 0.100;
L   = 0.290;
V_p = pi/4*(D_e^2 - D_i^2)*L;
M_tot = rho_p * V_p;

figure()
hold on
plot(pbar2438(:,1))
p_i = pbar2438(:,1);
p_max = max(p_i);
   
t = 0:dt:((length(p_i))-1)*dt;

idx = p_i >= 0.05*p_max;
idx(1:10) = 0; % Just to rule out initial oscillations (ignition transient)
t_ref = t(idx);
t_A = t_ref(1);
t_G = t_ref(end);

t_ref = t_A:dt:t_G; % i.e. over the Action Time

p_ref = trapz(t_ref, p_i(idx)) / (2*(t_G-t_A)) % reference pressure
yline(p_ref)

idx = p_i > p_ref;
idx(1:10) = 0; % Just to rule out initial oscillations (ignition transient)

p_cut = p_i(idx);
t_eff = t(idx);
t_b = t_eff(end)-t_eff(1); % (effective) Burning time
p_eff = trapz(t_eff, p_cut) / t_b % effective pressure

yline(p_eff)

r_b = (w*1000) / t_b % Burning rate [mm/s]

A_t = A_t_vec(1);
Cstar = 1e+5*A_t * trapz(t_eff, p_cut) / M_tot


 
