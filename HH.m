%% Basic variables

ms = 0.001;
mV = 0.001;
msiemens = 0.001;
ufarad = 0.000001;

u = 0;


%% Simulation loop

%Initial values
dt = 0.01*ms; %Microseconds

u_i = -65*mV;
x_i = [u_i; h0(u_i); m0(u_i); n0(u_i)];


plot_x = zeros(4, 100000);

for i = 1:100000

    plot_x(:,i) = x_i;

    x2 = x_i + dxdt(i*dt, x_i)*dt;
    
    x_i = x2;
    
end

hold on
plot(plot_x(1,:))



%% Complete model



function dxdt = dxdt(t, x)
    % Units
    ms = 0.001;
    mV = 0.001;
    msiemens = 0.001;
    ufarad = 0.000001;

    % Ion channel parameters
    El = -65 * mV;
    EK = -77 * mV;
    ENa = 55 * mV;
    gl = 0.3 * msiemens;
    gK = 35 * msiemens;
    gNa = 40 * msiemens;
    C = 1 * ufarad;
    

    %  u=x(1), h = x(2), m = x(3), n = x(4)
    dxdt = [-(1/C)*(gNa*x(3)^3*x(2)*(x(1)-ENa) + gK*x(4)^4*(x(1)-EK) + gl*(x(1)-El)) + (1/C)*Iex(t);
           alphah(x(1))*(1-x(2))-betah(x(1))*x(2);
           alpham(x(1))*(1-x(3))-betam(x(1))*x(3);
           alphan(x(1))*(1-x(4))-betan(x(1))*x(4)];


    %h 
    function alphah = alphah(u)
        u = u/mV;
        alphah = 0.25*exp(-(u+90)/12)/ms;
    end

    function betah = betah(u)
        u= u/mV;
        betah = 0.25*exp((u+62)/6)/exp((u+90)/12)/ms;
    end


    %m 
    function alpham = alpham(u)
        u= u/mV;
        alpham = 0.183*(u+35)/(1-exp(-(u+35)/9))/ms;
    end

    function betam = betam(u)
       u=u/mV;
       betam = -0.124*(u+35)/(1-exp((u+35)/9))/ms;
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
        I = 0.000001;
    
    end

end



function  m0 = m0(u)
    m0 = alpham(u)/(alpham(u)+betam(u));
end

function  n0 = n0(u)
    n0 = alphan(u)/(alphan(u)+betan(u));
end

function  h0 = h0(u)
    h0 = alphah(u)/(alphah(u)+betah(u));
end

%h 
function alphah = alphah(u)
  mV = 0.001;
  alphah = 0.25*exp(-(u/mV+90)/12);
 end

function betah = betah(u)
    mV = 0.001;
    betah = 0.25*exp((u/mV+62)/6)/exp((u/mV+90)/12);
end


%m 
function alpham = alpham(u)
    mV = 0.001;
    alpham = 0.183*(u/mV+35)/(1-exp(-(u/mV+35)/9));
end

function betam = betam(u)
    mV = 0.001;
    betam = -0.124*(u/mV+35)/(1-exp((u/mV+35)/9));
end

%n
function alphan = alphan(u)
    mV = 0.001;
    alphan = 0.02*(u/mV-25)/(1-exp(-(u/mV-25)/9));
end

function betan = betan(u)
    mV = 0.001;
    betan = -0.002*(u/mV-25)/(1-exp((u/mV-25)/9));
end










