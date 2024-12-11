function MamisSSA()
% Simulate a two-state model of gene expression
import Gillespie.*

%% Reaction network:
%   1. transcription:       0       --kR--> mRNA
%   2. translation:         mRNA    --kP--> mRNA + protein
%   3. mRNA decay:          mRNA    --gR--> 0
%   4. protein decay:       protein --gP--> 0

% ODE Function

% Parameters
p.d = 0.81160; %TA -> FD {CHANGES}
p.r = 1/(2.5); %S -> S TA {set} **1/2.5 d-1
p.lambda = 24/30;%TA -> 2TA {set} **1/30 h-1
p.gamma = 1/3.5; %FD -> dies  ** 1/3.5
k = 18; % 
qta = -p.r * k / (-p.d + p.lambda)
qfd = p.d / p.gamma * qta

% ODE Function
function dmdt = NewModelODE(~,m,p)
    s = m(1);
    ta = m(2);
    fd = m(3);
    dsdt = 0;
    dtadt = p.r * s + (p.lambda - p.d) * ta;
    dfddt = p.d * ta - p.gamma * fd;
    dmdt = [dsdt; dtadt; dfddt];
end

% Initial conditions
m0 = [18; 0; 0]; 

% Time span for simulation
tspan = [0 1000];

% Solve ODE
[t, m] = ode45(@(t, m) NewModelODE(t, m,p), tspan, m0);

% Plot results
figure;
plot(t, m(:, 1), '-r', 'LineWidth', 2); hold on;
plot(t, m(:, 2), '-g', 'LineWidth', 2);
plot(t, m(:, 3), '-b', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('pop');
legend('S', 'TA', 'FD');
title('Dynamics of the Model of Original Paper for Ntot=2400');
grid on;





%% Initial state
tspan = [0, 500]; %seconds
x0    = [18, 0, 0];     %mRNA, protein

%% Specify reaction network
pfun = @propensities_3state;
stoich_matrix = [ 0  1  0  %S -> S TA
                  0  1  0  %TA -> 2TA
                  0 -1  1  %TA -> FD
                  0  0 -1]; %FD -> dies

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
figure();
stairs(t,x); set(gca,'XLim',tspan);
title('Gillespie Model of Original Paper for Ntot=2400');
xlabel('time (days)');
ylabel('pop');
legend({'S','TA', 'FD'});

end


function a = propensities_3state(x, p)
% Return reaction propensities given current state x
S    = x(1);
TA = x(2);
FD = x(3);

a = [p.r*S;       
     p.lambda*TA;
     p.d*TA;
     p.gamma*FD;];   %protein decay
end