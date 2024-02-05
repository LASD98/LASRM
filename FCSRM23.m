clear
close all
clc

DATA = load('tracesbar1.mat');

dt = 1e-3; % Sampling rate = 1000 Hz
w = 0.03; % Web thickness [m]
A_t_vec = pi/4 * (1e-3*[28.80 25.25 21.81]').^2; % [High Medium Low] pressure [m^2]
% C_star_id = [1548.3 1553.7 1559.2]; % [m/s], from CEA, frozen
C_star_id = [1575.7 1579.6 1583.5]; % [m/s], from CEA, equilibrium
rho_p = 1 / (0.68/1950 + 0.18/2700 + 0.14/920);
D_e = 0.16;
D_i = 0.10;
L   = 0.29;
V_p = pi/4*(D_e^2 - D_i^2)*L;
M_tot = rho_p * V_p;

fn = fieldnames(DATA);
for k=1:9

    p = DATA.(fn{k});
    t = 0:dt:(length(p(:,1))-1)*dt;

    subplot(3,3,k)
    hold on

    for i = 1:3
        
        A_t = A_t_vec(i);
        
        p_i = p(:,i);
        pTrace = plot(t, p_i, 'lineWidth', 1);
        xlabel('t [s]')
        ylabel('p [bar]')
        p_max = max(p_i);
       
        % reference pressure
        idx = p_i >= 0.01*p_max;
        t_ref = t(idx);
        t_A = t_ref(1);
        t_G = t_ref(end);
    
        t_ref = t_A:dt:t_G; % i.e. over the Action Time
    
        p_ref = trapz(t_ref, p_i(idx)) / (2*(t_G-t_A));
    
        % effective pressure
        idx = p_i > p_ref;
    
        p_cut = p_i(idx);
        t_eff = t(idx);
        t_b = t_eff(end)-t_eff(1); % (effective) Burning time
        p_eff = trapz(t_eff, p_cut) / t_b;
        
        % burning rate
        r_b = (w*1000) / t_b; % Burning rate [mm/s]
        
        p_mat(i,k) = p_eff;
        r_b_mat(i,k)   = r_b;
        
        % experimental C_star    
        C_star_mat(i,k) = 1e+5*A_t * trapz(t_eff, p_cut) / M_tot;
        
    end
    
end

p_mat
r_b_mat
C_star_mat

C_star_vec = C_star_mat(:);
C_star_mean = mean(C_star_vec)
C_star_id_mean = mean(C_star_id);
C_star_std = std(C_star_vec)
C_star_efficiency = C_star_mean/C_star_id_mean


%% Vieille's law

p_vec = p_mat(:);
r_b_vec = r_b_mat(:);
x_log = log(p_vec);
y_log = log(r_b_vec);

figure(2)
plot(p_vec, r_b_vec, 'o')
% plot(x_log, y_log, 'o')
xlabel('p_e_f_f')
ylabel('r_b')
hold on

% z = log(min(p_vec)):0.01:log(max(p_vec));
% fittingCoeffs = polyfit(x_log,y_log,1);
% Pz = polyval(fittingCoeffs,z);
% plot(z,Pz)

z = min(p_vec):0.01:max(p_vec);
fittingCoeffs = polyfit(p_vec,r_b_vec,1);
Pz = polyval(fittingCoeffs,z);
plot(z,Pz)

[mu_a, sig_a, mu_n, sig_n, ~] = Uncertainty(p_vec, r_b_vec)



%% MONTE-CARLO

for NSRM = 1:3    
A_t = A_t_vec(NSRM);

N = 30;

p_c_vec = [];
t_b_vec = [];

for i=1:N
    a = normrnd(mu_a, sig_a);
    
    for j=1:N
        n = normrnd(mu_n, sig_n);
        
        for k=1:N
            C_star = normrnd(C_star_mean, C_star_std);         
            
            x_c = L/2;
            y_c = w;
            t_b = 0;
            
            while (x_c > 0) && (y_c > 0)
                r_i = (D_e/2-y_c);
                A_b = 2*pi*r_i*(2*x_c(end)) + 2*pi*((D_e/2)^2-r_i^2); %  A_latint + 2*A_base
                p_c = (1e-8*rho_p*mu_a*A_b/A_t*C_star)^(1/(1-mu_n)); % [bar]
                r_b = mu_a*p_c^mu_n;

                t_b = t_b+dt;
                x_c = x_c-r_b*dt/1000;
                y_c = y_c-r_b*dt/1000;
            end

            p_c_vec = [p_c_vec p_c];
            t_b_vec = [t_b_vec t_b];
    
        end
    end
end


% pressure plots

figure()
subplot(1,3,1)
hold on

for l = 1:length(p_c_vec)
    plot(l, mean(p_c_vec(1:l)), '.')
end

xlabel('iterations')
ylabel('\mu(p_c)')
            
subplot(1,3,2)
hold on

for l = 1:length(p_c_vec)
    plot(l, std(p_c_vec(1:l)), '.')
end

xlabel('iterations')
ylabel('\sigma(p_c)')

subplot(1,3,3)
hold on
histfit(p_c_vec)

xlabel('p_c normal distribution')


% burning time plots

figure()
subplot(1,3,1)
hold on

for l = 1:length(t_b_vec)
    plot(l, mean(t_b_vec(1:l)), '.')
end

xlabel('iterations')
ylabel('\mu(t_b)')
            
subplot(1,3,2)
hold on

for l = 1:length(t_b_vec)
    plot(l, std(t_b_vec(1:l)), '.')
end

xlabel('iterations')
ylabel('\sigma(t_b)')

subplot(1,3,3)
hold on
histfit(t_b_vec)

xlabel('t_b normal distribution')
            
hold off            

end











