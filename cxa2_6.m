function [] = cxa2_6()
%This function inegrated the equation of motion of a mass on a cycloid path
%paramterized in terms of a cycloidal angle phi, which is the angle throguh
%which the circle tracing the cycloid path has rotated. It displays the
%time it takes for the mass to reach a point B located at phi=1.59*pi.

termPhi = 1.59*pi; %Creates the phi value at which the mass has reached B.
c = 0.33; %Cycloid constant
g = 9.81; %Gravity constant
tPrev = 0; %Holder variable for time at integration points

init = [0.0000001, sqrt(2*g/c)]; %Setting initial conditions

t = linspace(0, 2*pi, 10000); %Time points to integrate at

%Function that sets the terminal phi value as the stopping point
function [value, isterminal, direction] = events(t,y)
    value = y(1)-termPhi; %Sets stop value at terminal phi
    isterminal = 1; %Makes the condition terminal
    direction = 0; %All zeros
end

%Function that displays the final time at which the mass has reached B.
function status = printOut(t, y, flag)
   if strcmp(flag, 'done') %If done itertating
       disp(tPrev); %Display final time
   else
       tPrev = t(end); %Store latest time;
   end
   status = 0; %Continue running
end  

    %Second order differntial equation of the equation of motion phi dot
    function dp = cycloidDiffEq(t, p)
        dp = [p(2);...
             (g/c-p(2)^2/2)*sin(p(1))/(1-cos(p(1)))];
    end


options = odeset('Events', @events, 'OutputFcn', @printOut); %Setting opts
[~, res1] = ode45(@cycloidDiffEq, t, init, options); %Integrating
end

