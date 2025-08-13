%% Plots for physarum's paper
% Daniele Proverbio (daniele.proverbio@unitn.it)
% 15/05/2024


V = -10 : 0.01 : 10;
V_T = 9;
alpha = 0.05;
beta = 5;

r = 10;
%r1=20;

M1 = - (- beta * V + 0.5*(beta - alpha) * (abs(V+V_T) - abs(V - V_T)));
M2 = (V/V_T) .^(2*r + 1);
%M21 = (V/V_T) .^(2*r1 + 1);

figure()
hold on
plot(V, M1, linewidth = 3)
%plot(V, M21, linewidth = 3)

%%
y = -1:0.01:1;
kappa = 800; % 800 : close to ideal; 300: farther from ideal
gamma = 0.5 ; % heterogeneity
u = zeros(size(y)); % arc tensions
epsi = 10e-5; % coefficient

% Piecewise-linear threshold function
u(y > +gamma) = kappa .* (y(y > +gamma) - gamma) + gamma .* epsi;
u(y < -gamma) = kappa .* (y(y < -gamma) + gamma) - gamma .* epsi;
u(y >= -gamma & y <= +gamma) = y(y >= -gamma & y <= +gamma) .* epsi;

V1 = -1 : 0.01 : 1;
beta1 = 800;
alpha1 = 0.00001;
V_T1 = 0.5;
M3 = - (- beta1 * V1 + 0.5*(beta1 - alpha1) * (abs(V1+V_T1) - abs(V1 - V_T1)));


figure()
hold on
plot(y,u)
plot(V1, M3)


%% NB: to plot the slime's development on the network
% I simply used 
% imagesc(obj(:,:,7900))
% colormap winter
% on the obj output of dynamics_systems_simulations