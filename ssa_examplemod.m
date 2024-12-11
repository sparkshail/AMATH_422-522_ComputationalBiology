function ssa_examplemod()
% Simulate a two-state model of gene expression
import Gillespie.*

%% Reaction network:
%   1. transcription:       0       --kR--> mRNA
%   2. translation:         mRNA    --kP--> mRNA + protein
%   3. mRNA decay:          mRNA    --gR--> 0
%   4. protein decay:       protein --gP--> 0

% ODE Function

% Parameters
p.rs = 1/(10); %S -> 2S, lit based 
p.d = 0.66; %TA -> FD {CHANGES}
p.r = 1/(2.5); %S -> S TA {set} **1/2.5 d-1
p.gs = 1/(10); %S -> dies 
p.lambda = 24/30;%TA -> 2TA {set} **1/30 h-1
p.lambdat = 0.15; %TA -> dies {CHANGES}
p.gamma = 1/3.5; %FD -> dies  ** 1/3.5
k = 18; % {set}
qta = -p.r * k / (-p.d + p.lambda-p.lambdat)
qfd = p.d / p.gamma * qta

% ODE Function
function dmdt = NewModelODE(~, m, p, k)
    s = m(1);
    ta = m(2);
    fd = m(3);
    dsdt = p.rs * s - p.rs / k * s.^2;
    dtadt = p.r * s + (p.lambda - p.d) * ta - p.lambdat * ta;
    dfddt = p.d * ta - p.gamma * fd;
    dmdt = [dsdt; dtadt; dfddt];
end

% Initial conditions
m0 = [18; 0; 0]; 

% Time span for simulation
tspan = [0 1000];

% Solve ODE
[t, m] = ode45(@(t, m) NewModelODE(t, m, p, k), tspan, m0);

% Plot results
figure;
plot(t, m(:, 1), '-r', 'LineWidth', 2); hold on;
plot(t, m(:, 2), '-g', 'LineWidth', 2);
plot(t, m(:, 3), '-b', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('pop');
legend('S', 'TA', 'FD');
title('Dynamics of the Modified Model for Ntot=2400');
grid on;





%% Initial state
tspan = [0, 500]; %seconds
x0    = [18, 0, 0];     %mRNA, protein

%% Specify reaction network
pfun = @propensities_3state;
stoich_matrix = [ 1  0  0   %S -> 2S
                  0  1  0  %S -> S TA
                 -1  0  0  %S -> dies
                  0  1  0  %TA -> 2TA
                  0 -1  0 %TA -> dies
                  0 -1  1  %TA -> FD
                  0  0 -1]; %FD -> dies

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tspan, x0, p);
%[t,x] = firstReactionMethod(stoich_matrix, pfun, tspan, x0, p);

%% Plot time course
figure();
stairs(t,x); set(gca,'XLim',tspan);
xlabel('time (days)');
ylabel('pop');
title('Gillespie Model of Modified model for Ntot=2400');
legend({'S','TA', 'FD'});

end


function a = propensities_3state(x, p)
% Return reaction propensities given current state x
S    = x(1);
TA = x(2);
FD = x(3);

a = [p.rs*S;            %transcription
     p.r*S;       %translation
     p.gs*S;       %mRNA decay
     p.lambda*TA;
     p.lambdat*TA;
     p.d*TA;
     p.gamma*FD;];   %protein decay
end