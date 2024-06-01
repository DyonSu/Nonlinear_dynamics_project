%% Basic variables

ms = 0.001;
mV = 0.001;
msiemens = 0.001;
ufarad = 0.000001;

u = 0;


%% Simulation loop

%Initial values
dt = 0.001*ms; %Microseconds

u_i = -60*mV; % Initial voltage
x_i = [u_i; m_0(u_i); h_0(u_i); s_0(u_i); n_0(u_i)];


simulation_time = 0.5; %seconds
T = round(simulation_time/dt); % The amount of time steps needed
plot_x = zeros(5, T);

percent = round(T/1000); % Timesteps of 0.1% of progress
count = 0; %For tracking progress

%Euler
% for i = 1:100000
% 
%     plot_x(:,i) = x_i;
% 
%     x2 = x_i + dxdt(i*dt, x_i)*dt;
%     
%     x_i = x2;
%     
% end

%rk2
% for i = 1:10000
% 
%     plot_x(:,i) = x_i;
%     
%     % Compute k1
%     k1 = dxdt(i*dt, x_i);
%     
%     % Estimate the state at the midpoint
%     x_mid = x_i + 0.5 * k1 * dt;
%     
%     % Compute k2 using the midpoint state
%     k2 = dxdt(i*dt + 0.5*dt, x_mid);
%     
%     % Update the state
%     x_i = x_i + k2 * dt;
%     
% end

%rk4
for i = 1:T

    plot_x(:, i) = x_i;
    
    % Compute k1
    k1 = dxdt(i * dt, x_i);
    
    % Compute k2
    k2 = dxdt(i * dt + dt / 2, x_i + k1 * dt / 2);
    
    % Compute k3
    k3 = dxdt(i * dt + dt / 2, x_i + k2 * dt / 2);
    
    % Compute k4
    k4 = dxdt(i * dt + dt, x_i + k3 * dt);
    
    % Update state
    x_i = x_i + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6;

    if mod(i, percent)==0
        clc
        count = count+0.1;
        
        disp(count);
    end

    % Check for divergence (optional)
    if abs(x_i(1)) > 100 * mV % Arbitrary threshold for membrane potential
        warning('Divergence detected at iteration %d', i);
        break;
    end
     
end


hold on
plot(plot_x(1,:))
%plot(plot_x(4,:))


%% Complete model



function dxdt = dxdt(t, x)

    % Units
    ms = 0.001; 
    mV = 0.001;
    mA = 0.001;
    pA = 10^-12;
    msiemens = 0.001;

    cm2 = 10^-4; % m2/cm2
    
    % Passive parameters
    Membrane_resistance = 10^4; %Ohms (這裡有刻意降低, 為了讓模型不會發散)
    Membrane_time_constant = 15*ms;  % From: Aradi I and Soltesz I. Modulation of network behaviour by changes in variance in interneuronal properties. J Physiol 538: 227–251, 2002. 
    gl = 1/Membrane_resistance; %Siemens/m2 
    C = Membrane_time_constant/Membrane_resistance; % farads/m2

    % Ion channel parameters
    El = -60 * mV;
    EK = -77* mV;
    ENa = 50 * mV;
    gK = 35 * msiemens/cm2; % 35 mS/cm2 *S/mS *cm2/m2 => Units of S/m2 (i.e. SI units)
    gNa = 200 * msiemens/cm2; % same here
    
 
    

    %  u=x(1), m = x(2), h= x(3), s = x(4), n = x(5)
    dxdt = [-(1/C)*(gNa*x(2)^3*x(3)*x(4)*(x(1)-ENa) + gK*x(5)^4*(x(1)-EK) + gl*(x(1)-El)) + (1/C)*Iex(t);
           -(1/tau_m(x(1)))*(x(2)-m_0(x(1)));
           -(1/tau_h(x(1)))*(x(3)-h_0(x(1)));
           -(1/tau_s(x(1)))*(x(4)-s_0(x(1)));
            alphan(x(1))*(1-x(5))-betan(x(1))*x(5);
           ];

    %m 
    function m_0 = m_0(u)
        u= u/mV;
        m_0 = 1/(1+exp(-(u+21.2)/4.9));
    end

    function tau_m = tau_m(u)
       u=u/mV;
       tau_m = 0.15*ms;
    end


    %h 
    function h_0 = h_0(u)
        u = u/mV;
        h_0 = 1/(1+exp((u+39.7)/7.7));
    end

    function tau_h = tau_h(u)
        u= u/mV;
        tau_h = 20.1*exp(-0.5*((u+61.4)/32.7)^2)*ms;
    end

    %s
    function s_0 = s_0(u)
        u = u/mV;
        s_0 = 1/(1+exp((u+46.1)/5.4));
    end

    function tau_s = tau_s(u)
        u=u/mV;
        tau_s= 1000*106.7*exp(-0.5*((u+52.7)/18.3)^2)*ms;
    end

    %n
    function alphan = alphan(u)
        u= u/mV;
        alphan = 0.02*(u-25)/(1-exp(-(u-25)/9))/ms;
    end

    function betan = betan(u)
        u= u/mV;
        betan = -0.002*(u-25)/(1-exp((u-25)/9))/ms;
    end

    function I = Iex(t)
    % Constant current injected at some interval of time
        if t>=50*ms && t<=250*ms
            I = 0.02*mA;
        else
            I = 0;
        end
    end



end

%m 
function m_0 = m_0(u)
    mV = 0.001;
     u= u/mV;
    m_0 = 1/(1+exp(-(u+21.2)/4.9));
end



%h 
function h_0 = h_0(u)
    mV = 0.001;
    u = u/mV;
    h_0 = 1/(1+exp((u+39.7)/7.7));
end



%s
function s_0 = s_0(u)
    mV = 0.001;
    u = u/mV;
    s_0 = 1/(1+exp((u+46.1)/5.4));
end

function  n_0 = n_0(u)
    n_0 = alphan(u)/(alphan(u)+betan(u));
end

%n
function alphan = alphan(u)
    mV = 0.001;
    ms = 0.001;
    u = u/mV;
    alphan = 0.02*(u-25)/(1-exp(-(u-25)/9))/ms;
end

function betan = betan(u)
    mV = 0.001;
    ms = 0.001;
    u = u/mV;
    betan = -0.002*(u-25)/(1-exp((u-25)/9))/ms;
end









