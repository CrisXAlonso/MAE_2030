% This function first numerically intergrates two spring functions using
% the ODE45 solver. The first function represents a linear spring with an
% attached mass. The second is the same system but now with a linear spring
% of the form -k/m*(x-x0)-c/m*(x-x0)^3. t is the time interval over which
% the integration is carried, xInit is the vector of initial position and
% velocity, c is the nonlinear spring constant, m is the mass attached to
% the spring, x0 is the unstretched length of the spring, and k is the
% spring constant. Secondly it integrates a simple harmonic oscillators
% equation of motion as per problem 2.17 of the Engineering Dynamics
% textbook and plots the position vs time and velocity vs time graphs of
% that simple harmonic osciallator. Lastly, the function plots the
% difference between the analytical and numerical solutions for problem
% 2.17.

function [res1, res2, res3] = cxa2_2(t, xInit1, c, m1, x0, k1, a, m2,...
    k2, xInit2)


t10 = linspace(0, 10, 500);

%Sets the values given in the book or in the PSet as defaults
if ~exist('t','var') || isempty(t)
    t = linspace(0, 20, 200); % in seconds
end
if ~exist('xInit','var') || isempty(xInit1)
    xInit1 = [0.4, 0]; % in meters and meters/second
end
if ~exist('c','var') || isempty(c)
    c = 5;  % in N/m^3
end
if ~exist('m1','var') || isempty(m1)
    m1 = 1; % in kg
end
if ~exist('x0','var') || isempty(x0)
    x0 = 0.25; % in meters
end
if ~exist('k','var') || isempty(k1)
    k1 = 1; % in N/m
end
if ~exist('a','var') || isempty(a)
    a = 1; % in m/s^2
end
if ~exist('m2','var') || isempty(m2)
    m2 = 0.5; % in m
end
if ~exist('k2','var') || isempty(k2)
    k2 = 3; % in N/m
end
if ~exist('xInit2','var') || isempty(xInit2)
    xInit2 = [1, 0.5]; % in m and m/s
end

    % Each of the differential equations are listed below:
    function dx = linearSpringSys(t, x)
        dx = [x(2);...
            -k1/m1*(x(1)-x0)];
    end

    function dx = nonLinearSpringSys(t, x)
        dx = [x(2);...
            -k1/m1*(x(1)-x0)-c/m1*(x(1)-x0)^3];
    end

    function dx = simpleHarmonicOsc(t, x)
        dx = [x(2);...
            a - k2/m2*x(1)];
    end

    % Below is the analytical solution to problem 2.17:
    function x = analyticalSol(t)
        x = (1-m2*a/k2)*cos(sqrt(k2/m2)*t)+0.5*sqrt(m2/k2)*sin(sqrt(k2/m2)...
            *t)+(m2/k2)*a;
    end

    % Below is the function that finds the di

% Solves the differential equations
[~, res1] = ode45(@linearSpringSys, t, xInit1);
[~, res2] = ode45(@nonLinearSpringSys, t, xInit1);
[~, res3] = ode45(@simpleHarmonicOsc, t10, xInit2);

% Plotting:
hold all 

f  = figure(1);
plot(t,res1(:, 1),'--r',t,res2(:, 1),'-b', t, x0, 'k')
legend('Linear Spring','Nonlinear Spring', 'Unstretched Length')
set(f, 'Position', [1, 1, 1200, 800])
xlabel('Time (s)', 'FontSize', 14)
ylabel('Position (m)', 'FontSize', 14)
title('2.13: Linear and Nonlinear Spring Position Comparison', 'FontSize',...
    20)


h = figure(2);
set(h, 'Position', [1, 1, 1100, 400])
s1 = subplot(2, 1, 1);
plot(t10, res3(:, 1), 'k')
title('2.17: Position of Simple Harmonic Oscillator vs Time')
xlabel('Time (s)')
ylabel('Position (m)')
s2 = subplot(2, 1, 2);
plot(t10, res3(:, 2), 'r')
title('2.17: Velocity of Simple Harmonic Oscillator vs Time')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

j = figure(3);
plot(t10, analyticalSol(t10)-res3(:, 1)', 'g')
set(j, 'Position', [1, 1, 1200, 800])
xlabel('Time (s)', 'FontSize', 14)
ylabel('Difference in Position (m)', 'FontSize', 14)
title('2.17: Difference between Analytical and Numerical Solutions',...
    'FontSize', 20)
h1 = text(3.2, 0.002,...
    '(Results are very accurate but begin to deviate as time progresses)');
set(h1, 'FontSize', 12);
end