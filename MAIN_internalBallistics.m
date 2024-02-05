clear all
close all
clc

%set(0,'defaultTextInterpreter','latex');
%set(0,'DefaultAxesFontSize', 14);

DATA = load('tracesbar1.mat');


dt = 1e-3; % Sampling rate = 1000 Hz
w = .030; % Web thickness [m]
A_t_vec = pi/4 * (1e-3*[28.80 25.25 21.81]').^2; % Respectively [High Medium Low] pressure [m^2]
C_star_id = [1548 1554 1560]; %[m/s]
rho_p = 1 / (.68/1950 + .18/2700 + .14/920);
D_e = .160;
D_i = .100;
L   = .290;
V_p = pi/4*(D_e^2 - D_i^2)*L;
M_tot = rho_p * V_p;


h1 = figure;
set(h1, 'Units', 'Normalized', 'OuterPosition', [.0 .0 1 1])
set(h1, 'Name', 'Pressure Traces');

fn = fieldnames(DATA);
for k=1:numel(fn) % Cycle over each 3-traces set
    if( isnumeric(DATA.(fn{k})) )
        p = DATA.(fn{k});
        t = [0:dt:(length(p(:,1))-1)*dt];




subplot(ceil(numel(fn)/3),ceil(numel(fn)/3),k)
hold on
legendInf = {};
legendObjs = [];


for i = 1:size(p,2)
    
    A_t = A_t_vec(i);
    
    p_i = p(:,i);
    pTrace = plot(t, p_i, 'lineWidth', 1);
    xlabel('t [s]')
    ylabel('p [bar]')
    
    p_max = max(p_i);
    
    legendInf{i} = sprintf('$@ \\, \\, p_\\mathrm{max} = %.3f \\, \\, \\mathrm{bar}$', p_max);
    legendObjs(i) = pTrace;
    legend(legendObjs, legendInf, 'AutoUpdate', 'off', 'FontSize', 11, 'Interpreter', 'latex');
   




    % 1. Identify REFERENCE pressure, from the A,G points: trace above 1%
    % of p_max
    idx = p_i >= 0.01*p_max;
    idx(1:10) = 0; % Just to rule out initial oscillations (ignition transient)
    t_ref = t(idx);
    t_A = t_ref(1);
    t_G = t_ref(end);

    t_ref = [t_A:dt:t_G]; % i.e. over the Action Time

    p_ref = trapz(t_ref, p_i(idx)) / (2*(t_G-t_A));

    % 2. Identify the EFFECTIVE pressure, from the B,E points: trace above
    % p_ref
    idx = p_i > p_ref;
    idx(1:10) = 0; % Just to rule out initial oscillations (ignition transient)

    p_cut = p_i(idx);
    t_eff = t(idx);
    t_b = t_eff(end)-t_eff(1); % (effective) Burning time
    p_eff = trapz(t_eff, p_cut) / t_b;
    
    
    % 3. Compute burning rate
    r_b = (w*1000) / t_b; % Burning rate [mm/s]
    
    p_mat(i,k) = p_eff;
    r_b_mat(i,k)   = r_b;
    
    % 4. Compute C_star    
    C_star_mat(i,k) = 1e+5*A_t * trapz(t_eff, p_cut) / M_tot;
    

end
drawnow


    
        end
end


p_mat
r_b_mat
C_star_mat

h1.Renderer = 'Painter';
printFigure(h1, 'pressureTraces_BATES')


h2 = figure;
p_vec = p_mat(:);
r_b_vec = r_b_mat(:);

    x_log = log(p_vec);
    y_log = log(r_b_vec);
    plot(x_log, y_log, 'o', 'MarkerSize', 5)
    xlabel('$\mathrm{log}(p_\mathrm{eff})$')
    ylabel('$\mathrm{log}(r_b)$')
    axis equal
    grid on
    hold on
    
    z = [log(min(p_vec)):.01:log(max(p_vec))];
    fitOrd = 1; % i.e. linear fitting
    fittingCoeffs = polyfit(x_log,y_log,fitOrd)
    Pz = polyval(fittingCoeffs,z);
    plot(z,Pz)
    
    % RMK.:
    % r_b = a*p^n -> log(r_b) = log(a) + n*log(p)

    % Vieille's Law box
    vieilleBoxLocation = [.3 .75 .1 .1];
    vieilleBox = annotation('textbox',vieilleBoxLocation,'String', sprintf('$ r_{b} = %.3fp^{%.3f} \\,\\, \\mathrm{mm/s} $', exp(fittingCoeffs(2)), fittingCoeffs(1)), 'Interpreter', 'latex');
    vieilleBox.FontSize = 16;
    vieilleBox.BackgroundColor = 'w';
    
    mu_a = exp(fittingCoeffs(2));
    mu_n = fittingCoeffs(1);
    r_b_fun = @(a,p,n) a .* p .^n;
    
    [~, sig_a, ~, sig_n, ~] = Uncertainty(p_vec, r_b_vec)
    
    
    mu_C = mean(C_star_mat(:))
    sig_C = std (C_star_mat(:))
    
    
printFigure(h2, 'r_b_vieilleDataFitting')

    
    
    
%% Regression tracking:
%  the upper-right corner will be taken as a refernece

x_c_0 = L/2;
y_c_0 = w;
t0 = 0;

x_c = [x_c_0];
y_c = [y_c_0];
t = [t0];
p_vec = [];

motNum = 1;

a = 1.735; % **********************
n = 0.381; % **********************
C_star = 1550; % **********************



while (x_c(end) > 0) & (y_c(end) > 0)
    
    r_i = (D_e/2-y_c(end));
    A_b = 2*pi*r_i*(2*x_c(end)) + 2*pi*((D_e/2)^2-r_i^2); %  A_latInt + 2*A_base
    p_c = 1e-5*(10^(-5*n-3) * rho_p*a*A_b/A_t_vec(motNum)*C_star)^(1/(1-n)) %[bar] %****************
    
    p_vec = [p_vec; p_c];
    r_b = r_b_fun(a,p_c,n);
    
    t = [t (t(end)+dt)];
    x_c = [x_c (x_c(end)-r_b/1000*dt)];
    y_c = [y_c (y_c(end)-r_b/1000*dt)];    
end


figure
N = length(y_c);
% xlim([-x_c(1) x_c(1)  ])
% ylim([0       y_c(1)])
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
axis equal
hold on

pGeom = 0;
timeStep = 50;

    % Time-box
    timeBoxLocation = [.3 .75 .1 .1];
    timeBox = annotation('textbox',timeBoxLocation,'String', sprintf('$@ p = %.2f \\, \\, \\mathrm{bar}: \\, \\, \\, t = 0 \\,\\, \\mathrm{s}$', p_c), 'Interpreter', 'latex');
    timeBox.FontSize = 40;
    timeBox.BackgroundColor = 'w';
    
for i = [1:timeStep:N]
    if i > timeStep+1
       delete(pGeom)
    end
    delete(timeBox)
   grainGeom = [x_c(i) , 0     ;
                x_c(i) , y_c(i);
                -x_c(i), y_c(i);
                -x_c(i), 0     ];
   if i > 1
       plotStyle = 'r';
   else
       plotStyle = 'b';
   end
   pGeom = plot(grainGeom(:,1), grainGeom(:,2), plotStyle);
   timeBox = annotation('textbox',timeBoxLocation,'String', sprintf('$@ p = %.2f \\, \\, \\mathrm{bar}: \\, \\, \\, t = %.2f \\,\\, \\mathrm{s}$', p_c, t(i)), 'Interpreter', 'latex');

   drawnow
end





%% MONTE CARLO

SRM_names = {'lowP', 'mediumP', 'highP'};

for motNum = 1:size(p,2)    
    % For each motor class we will have a burning time in a different time
    % range
 
    A_t = A_t_vec(motNum);

    
N = 30;

% x_a = [(mu_a - 5 * sig_a) : (sig_a / N) : (mu_a + 5 * sig_a)];
% gauss_a = normpdf(x_a, mu_a, sig_a); % Normal probability density function (pdf).
% 
% x_n = [(mu_n - 5 * sig_n) : (sig_n / N) : (mu_n + 5 * sig_n)];
% gauss_n = normpdf(x_n, mu_n, sig_n);
% 
% 
% x_C = [(mu_C - 5 * sig_C) : (sig_C / N) : (mu_C + 5 * sig_C)];
% gauss_C = normpdf(x_C, mu_C, sig_C);


% figure;
% subplot(1,3,1)
% plot(x_a, gauss_a);
% subplot(1,3,2)
% plot(x_n, gauss_n);
% subplot(1,3,3)
% plot(x_C, gauss_C);


h3 = figure;
set(h3, 'Units', 'Normalized', 'OuterPosition', [.1 .35 .8 .5])
count = 0;
t_b_vec = [];


parfor i=1:N
    i
    a = normrnd(mu_a, sig_a)
    
    for j=1:N
        n = normrnd(mu_n, sig_n)
        
        for k=1:N
            C_star = normrnd(mu_C, sig_C)

         
            
            %count = count+1;
            x_c = L/2;
            y_c = w;
            t_b = 0;
            
            while (x_c > 0) & (y_c > 0)
                r_i = (D_e/2-y_c);
                A_b = 2*pi*r_i*(2*x_c(end)) + 2*pi*((D_e/2)^2-r_i^2); %  A_latInt + 2*A_base
                p_c = (1e-8 * rho_p*mu_a*A_b/A_t*C_star)^(1/(1-mu_n)); % [bar] %****************
                r_b = mu_a*p_c^mu_n;

                t_b = t_b+dt;
                x_c = x_c-r_b/1000*dt;
                y_c = y_c-r_b/1000*dt;
            end
            t_b_vec = [t_b_vec t_b];
%             subplot(1,2,1)
%             if i == 1
%                 xlabel('$\mathrm{iterations}$')
%                 ylabel('$\overline{t_b}$')
%             end
%             hold on
%             plot(count, mean(t_b_vec), 'b.', 'MarkerSize', 5)
%             drawnow
%             subplot(1,2,2)
%             if i == 1
%                 xlabel('$\mathrm{iterations}$')
%                 ylabel('$\sigma_{t_{b}}$')
%             end
%             hold on
%             plot(count, std(t_b_vec), 'b.', 'MarkerSize', 5)
%             drawnow
    
end
end
end



            subplot(1,2,1)
            hold on

            for count = 1:length(t_b_vec)
                plot(count, mean(t_b_vec(1:count)), 'b.', 'MarkerSize', 5)
            end
            
            xlabel('$\mathrm{iterations}$')
            ylabel('$\overline{t_b}$')
            
            drawnow
            
            subplot(1,2,2)
            hold on
            
            for count = 1:length(t_b_vec)
                plot(count, std(t_b_vec(1:count)), 'b.', 'MarkerSize', 5)
            end
            
            xlabel('$\mathrm{iterations}$')
            ylabel('$\sigma_{t_{b}}$')
            
            drawnow
            
            
            


h3.Renderer = 'Painter';
printFigure(h3, sprintf('MC_%s', SRM_names{motNum}))

end
