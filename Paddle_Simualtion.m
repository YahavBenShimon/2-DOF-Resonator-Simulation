
% Initialization
clear all;
clc;
tic;

% ================= Dimensional Variables =================

% Scalar Variables
rho = 19300;           % Gold density [kg/m^3]
l = 1e-6;              % Paddle length [m]
g = 1e-6;              % Initial gap [m]
w = 300e-9;            % Paddle width [m]
h = 85e-9;             % Paddle thickness [m]
m = l * w * h * rho;   % Paddle mass [kg]
I = (m/12) * (h^2 + l^2);  % Paddle moment of inertia [kg*m^2]
V_dc = 0.02;           % DC voltage [V]
V_ac = 0.005;          % AC voltage [V]
e = 2.3E-07;           % Eccentricity [m]
S2D = 2.2;             % Static To Dynamic torsion stiffness coeff.
k_T = S2D * 5.842e-14; % Torsion spring constant [N*m]
k_B = 4.18e-1;         % Bending spring constant [N/m]
Q_T = 100;             % Torsion quality factor
Q_B = 50;              % Bending quality factor
c_T = sqrt(I * k_T) / Q_T;  % Torsion damping [N*s]
c_B = sqrt(m * k_B) / Q_B;  % Bending damping [N*s/m]
T = sqrt(I / k_T);     % Natural frequency inverse [s]
F_coff = -1.353e-5;    % Electrical Energy Derivative w.r.t y [F/m]
M_coff = -3.378e-12;   % Electrical Energy Derivative w.r.t y [F/rad]

% Frequency Sweep Settings
frq_start = 1e6;            % Start frequency [Hz]
frq_step = 0.1e6;           % Step frequency [Hz]
frq_stop = 7e6;             % Stop frequency [Hz]
omega_start = frq_start * 2*pi; % Start angular frequency
omega_step = frq_step * 2*pi;   % Step angular frequency
omega_stop = frq_stop * 2*pi;   % Stop angular frequency

% Vector Variables
tspan = linspace(0, 100e-6, 10000); % Time Vector [s]
IC = [0, 0, 0, 0];                 % Initial Conditions (theta, dtheta, y, dy)

% Prepare storage for frequency response
Frq_respone_theta = zeros(length(frq_start:frq_step:frq_stop), 2);
Frq_respone_theta(:, 1) = frq_start:frq_step:frq_stop;
Frq_respone_y = Frq_respone_theta;

% ============== Nondimensional Variables ==============

tspan_non = tspan / T;
omega_start_non = omega_start * T;
omega_step_non = omega_step * T;
omega_stop_non = omega_stop * T;

% ==================== ODE Solver ====================

flag_theta = 0;
flag_y = 0;

for sweep_omega_non = omega_start_non:omega_step_non:omega_stop_non
    % Define the differential equations
    syms y(t) theta(t) Y;
    
    % Bending Equation
    Eq1 = diff(y,2) == -((c_B*T)/m)*diff(y) + ((T^2*k_B*e)/(g*m))*theta - ((T^2*k_B)/m)*y ...
          + (T^2/(g*m))*0.5*F_coff*(V_dc + V_ac*sin(sweep_omega_non*t))^2;
    
    % Angle Equation
    Eq2 = diff(theta,2) == -((T*c_T)/I)*diff(theta) - ((T^2*(k_T+k_B*e^2))/I)*theta + ((T^2*g*k_B*e)/I)*y ...
          + (T^2/I)*0.5*M_coff*(V_dc + V_ac*sin(sweep_omega_non*t))^2;
    
    % Convert ODEs to vector format and solve
    [M,Coff] = odeToVectorField(Eq1, Eq2);
    coupled_ODE = matlabFunction(M, 'Vars', {t, Y});
    [time,v_z] = ode45(coupled_ODE, tspan_non, IC);
    
    % Post-process to find steady-state amplitude
    index = find(abs([omega_start_non:omega_step_non:omega_stop_non] - sweep_omega_non) < 0.0001);
    
    Frq_respone_theta(index, 2) = 0.5*(180/pi)*(max(v_z(0.5*end:end,1)) - min(v_z(0.5*end:end,1)));
    Frq_respone_y(index, 2) = 0.5*g*10^9*(max(v_z(0.5*end:end,3)) - min(v_z(0.5*end:end,3)));
    
    % Determine maximum amplitude resonance for plotting
    if (index == 1)
        flag_theta = Frq_respone_theta(index, 2);
        flag_y = Frq_respone_y(index, 2);
    end
    
    if Frq_respone_theta(index, 2) > flag_theta
        time_theta_max = time * T;
        amp_theta_max = v_z(:, 1) * (180/pi);
        flag_theta = Frq_respone_theta(index, 2);
        Frq_theta = Frq_respone_theta(index, 1);
    end
    
    if Frq_respone_y(index, 2) > flag_y
        time_y_max = time * T;
        amp_y_max = v_z(:, 3) * g * 10^9;
        flag_y = Frq_respone_y(index, 2);
        Frq_y = Frq_respone_y(index, 1);
    end
end

% =================== Plotting Results ===================

% Plot the time response for theta
figure('Name', 'Theta');
plot(time_theta_max, amp_theta_max)
title(['\theta Time Response, f = ', num2str(Frq_theta*10^-6), ' (MHz)'])
xlabel('Time (sec)')
ylabel('\theta (Deg)')
grid off
set(gca,'fontsize',16)

% Plot the time response for y
figure('Name', 'Displacement');
plot(time_y_max, amp_y_max)
title(['y Time Response, f = ', num2str(Frq_y*10^-6), ' (MHz)'])
xlabel('Time (\mus)')
ylabel('y (nm)')
grid off
set(gca,'fontsize',16)

% Plot the frequency response for theta
figure('Name', 'Freq Theta');
plot(Frq_respone_theta(:,1)*10^-6, Frq_respone_theta(:,2))
title('\theta Frequency Response')
xlabel('Frequency (MHz)')
ylabel('\theta (Deg)')
grid off
set(gca,'fontsize',16)

% Plot the frequency response for y
figure('Name', 'Frq y');
plot(Frq_respone_y(:,1)*10^-6, Frq_respone_y(:,2))
title('y Frequency Response')
xlabel('Frequency (MHz)')
ylabel('Amplitude (nm)')
grid off
set(gca,'fontsize',16)

% Display elapsed time
toc
