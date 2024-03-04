% 1. Parameters structure 
parameters.n =50;                                     % number of cells
parameters.d = 100;                                   % depth (m)
parameters.Deltaz = parameters.d/parameters.n;        %grid width 
parameters.up = 0.0042*24;                            % velocity current - sinking velocity phytoplankton (m /day) 
parameters.ud = 0.042*24*10;                          % sinking velocity detritus (m /day) 
parameters.Dcoef= 43.2;                               % difussion (m^2/day) 

parameters.m= 0.001*24;                              %loss rate/mortality (1/d) 
parameters.Hi= 30*60*60*24/1000;                     %half saturation constant of light limited growht (mmol photons/m2/d)
parameters.gmax= 1.01;                               % max production rate (1/d)
parameters.kw= 0.19;                                 %background turbidity (1/m) %0.06; 
parameters.kp=  15 * 10^-12/10;                      % specific light attenuation of phytopl(m^2/mmol N)
parameters.Io = 350 *60*60*24/1000;                  %Incident light (mmol photons /m2/d) 
 
parameters.e = 0.8;                                  %nutrient recycling coeff (dimensionless)
parameters.Hn = 2.5;                                 %half saturation ctt of nutrient limited growth (mmol nutrients /m^3) antes 0.025
parameters.Nb = 20;                                  %Nutrient concentration at the bottom (mmol nutrient/m^3

parameters.gamma= 0.02;                              % the zooplankton quadratic grazing parameter (m^3/(mmolN d)) 
parameters.tau= 0.1;                                 % remineralization (d^-1) 

parameters.to= 150*0;                                %Adjusting seasonal light
parameters.B= 325*60*60*24/1000;                     %Adjusting seasonal light
parameters.A= 25*60*60*24/1000 ;                     %Seasonal maximum incident light (mmol photons /m2/d)

% 2. Set up grid
parameters.z = parameters.Deltaz /2 : parameters.Deltaz : (parameters.d - parameters.Deltaz /2);

%3. Set initial conditions
phi0 = zeros(1,parameters.n);
phi0(5) = 1;
N0 = zeros(1,parameters.n)+parameters.Nb;
D0= zeros(1,parameters.n);
 
%4. Solve ODE
[t, PND] = ode45(@(t, PND) derivative(t, PND, parameters), [0 800+365+365+365], [phi0, N0, D0]);



phi = PND(:, 1:parameters.n);
N = PND(:, parameters.n+1: 2*parameters.n);
D= PND(:,2*parameters.n+1: 3*parameters.n);

%5. Plot results
% Plot steady state 
figure;

subplot(1,3,1);
imagesc(t, parameters.z, phi');
xlabel('Time (d)');
ylabel('Depth (m)');
title('Phytoplankton concentration (mmol N/m^3)');
colorbar;

subplot(1,3,2);
imagesc(t, parameters.z, N');
xlabel('Time (d)');
ylabel('Depth (m)');
title('Nutrients concentration (mmol N/m^3)');
colorbar;

subplot(1,3,3);
imagesc(t, parameters.z, D');
xlabel('Time (d)');
ylabel('Depth (m)');
title('Detritus concentration (mmol N/m^3)');
colorbar;

% Plot limiting factors 
figure;
I = integral(t,parameters,phi(end,:));
I_lim= (I./(parameters.Hi+I));
N_lim= (N./(parameters.Hn+N));
g= parameters.gmax*min(I./(parameters.Hi+I), N./(parameters.Hn+N));
plot(I_lim(end,:), parameters.z, 'Color', [0.7 0.7 0.7],'LineWidth', 1.5, 'DisplayName', 'Light intensity');
hold on;
plot(N_lim(end,:), parameters.z, 'Color', [0.9 0.9 0.9], 'LineWidth', 1.5, 'DisplayName', 'Nutrient availability');
hold on;
plot(g(end,:), parameters.z, 'k--','LineWidth', 1, 'DisplayName', 'Phytoplankton growth (mmol N/m^3)');
ax = gca;
ax.YDir = 'reverse';
grid on;
xlabel('Values');
ylabel('Depth (m)');
legend('Light intensity', 'Nutrients availability','Phytoplancton growth (mmol N/m^3)'); 
legend('Location', 'southoutside');

% Plot sensitivity analysis
% SA - Nutrient recycling coefficient
figure;
subplot(1,2,1)
parameters.e_SA = [0.5, 1, 3, 5, 10];
for i = 1:length(parameters.e_SA)
    parameters.e = parameters.e_SA(i);

    % Solve ODE
    [t, PND] = ode45(@(t, PND) derivative(t, PND, parameters), [0 400], [phi0, N0, D0]);

    % Plot results
    phi = PND(:, 1:parameters.n); 
    color = [0 0.4470 0.7410;   
          0.1490 0.5451 0.1490;  
          0.4940 0.4940 0.4940;   
          0.6350 0.0780 0.1840;   
          0.4940 0.1840 0.5560];  
    plot(phi(end,:), parameters.z, 'Color', color(i,:),'LineWidth', 1.5);
    hold on;

end
ax = gca;
ax.YDir = 'reverse';
grid on; 
xlabel('Phytoplankton density (mmol/m^3)');
ylabel('Depth (m)');
legend('ε = 0.5', 'ε = 1','ε = 3', 'ε = 5', 'ε = 10'); 
legend('Location', 'southoutside');


% SA - Mortality
subplot(1,2,2);
parameters.m_SA = [0.001, 0.02, 0.05, 0.08, 0.1];
for i = 1:length(parameters.e_SA)
    parameters.m = parameters.m_SA(i);

    % Solve ODE
    [t, PND] = ode45(@(t, PND) derivative(t, PND, parameters), [0 400], [phi0, N0, D0]);

    % Plot results
    phi = PND(:, 1:parameters.n); 
    color = [0 0.4470 0.7410;   
          0.1490 0.5451 0.1490;   
          0.4940 0.4940 0.4940;   
          0.6350 0.0780 0.1840;   
          0.4940 0.1840 0.5560];  
    plot(phi(end,:), parameters.z, 'Color', color(i,:), 'LineWidth', 1.5);
    hold on;

end
ax = gca;
ax.YDir = 'reverse';
grid on;
xlabel('Phytoplankton density (mmol/m^3)');
ylabel('Depth (m)');
legend('m = 0.002 d^-^1', 'm = 0.01 d^-^1','m = 0.05 d^-^1', 'm = 0.08 d^-^1', 'm = 0.1 d^-^1'); 
legend('Location', 'southoutside');


% Derivative funciton  
function dPNDdt = derivative(t, PND, parameters)
    %Extract paramteters from struct
    n = parameters.n;
    Deltaz = parameters.Deltaz;
    up = parameters.up;
    ud = parameters.ud;
    Dcoef = parameters.Dcoef;
    m = parameters.m;   
    Hi = parameters.Hi;                         
    gmax = parameters.gmax;
    e = parameters.e;
    Hn = parameters.Hn;
    Nb = parameters.Nb;
    gamma = parameters.gamma; 
 
    phi = PND(1:parameters.n);
    N = PND(parameters.n+1:2*parameters.n);
    D = PND(2*parameters.n+1:3*parameters.n);

    %initialize fluxes phytoplankton  
    Jap = zeros (1, n+1);
    Jdp = zeros(1, n+1);

    %calculate fluxes phytoplankton
    for i = 2:n
        Jap(i) = up *phi(i-1);                          %Advective flux
        Jdp(i) = -Dcoef *((phi(i) - phi(i-1))/Deltaz);  %Diffusive flux
    end
    Jp = Jap + Jdp;

    % boundary fluxes phytoplancton
    Jp(1) = 0;
    Jp(n+1) = 0;

    %Initizlize fluxes nutrients
    Jdn = zeros(1, n+1);
    
    %Calcualte fluxes nutrients
    for i = 2:n
        Jdn(i) = -Dcoef *((N(i) - N(i-1))/Deltaz);  %Diffusive flux
    end 

    %boundary fluxes nutrients
    Jdn(1) = 0;
    Jdn(n+1) = -Dcoef * (Nb - N(n))/Deltaz; 

    %initialize fluxes detritus 
    Jad = zeros (1, n+1);
    Jdd = zeros(1, n+1);

    %calculate fluxes detritus
    for i = 2:n
        Jad(i) = ud *D(i-1);                         %Advective flux
        Jdd(i) = -Dcoef *((D(i) - D(i-1))/Deltaz);   %Diffusive flux
    end
    Jd = Jad + Jdd;
    
    %boundary fluxes detritus
    Jd(1)=0;
    Jd(n+1)= ud*D(end);

    %light intensity
    
    

    I = integral(t, parameters, phi);
   
    I_lim= (I./(Hi+I));
    N_lim= (N./(Hn+N));
    g= gmax*min(I_lim, N_lim);
    r=g-m;

    % Calcualte overall derivative 
    
    dphidt=zeros(1,n);
    dNdt=zeros(1,n);
    for i = 1:n
        
        dphidt(i) = r(i) * phi(i) - gamma*phi(i)^2 -(Jp(i+1) - Jp(i))/Deltaz; 

        dNdt(i) = - g(i) * phi(i) + e * D(i) - (Jdn(i+1)-Jdn(i))/Deltaz;   
        
        dDdt(i) =  m*phi(i) + gamma * phi(i)^2  - e*D(i) - (Jd(i+1)-Jd(i))/Deltaz;
    
    end
    dPNDdt = [dphidt, dNdt, dDdt]';
end

%light intensity 
function I = integral(t, parameters, phi)
    n =parameters.n;

    kw = parameters.kw;
    kp = parameters.kp;
    Io = parameters.Io; 
    Deltaz = parameters.Deltaz;

    to = parameters.to;
    B = parameters.B;
    A = parameters.A;

    int_turb = cumsum((kw+kp*phi)*Deltaz);
    Is= A * sin (2*pi/365*(t-to)+B);
    I = Io* exp(-int_turb);
    %I = Is.* exp(-int_turb);
    
     
end