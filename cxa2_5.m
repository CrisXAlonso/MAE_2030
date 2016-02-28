function [res1] = cxa2_5(r, l, init)
%This function solves the tetherball differential equation presented in
%Example 4.7 of the Kasdin and Paley dynamics textbook using the provided
%equations of motion and then plots the trajectory of the tetherball where
%r is the radius of the post, l is the length of the tether and init is a
%vecor containing the initial angle and the initial angular velocity of the
%ball

if ~exist('r','var') || isempty(r)
    r = 1; % In meters
end
if ~exist('l','var') || isempty(l)
     l=10; %In meters
end
if ~exist('init','var') || isempty(init)
    init = [0, 0.1]; % In rad and rad/s
end

termTheta = (l-1)/r; %Calculates the theta at which there is not more
                     %tether left, this becomes the stopping condition
tPrev = 0; %Holder variable to display final time
 
%This function sets the terminal event for the ode45, telling it to stop
%once the terminal theta is reached
function [value, isterminal, direction] = events(t,y)
    value = y(1)-termTheta; %Sets stop value at terminal theta
    isterminal = 1; %Makes the condition terminal
    direction = 0; %All zeros
end
 
%This function displays the time at which the tether hits the post by
%storing the previous time on all iterations, and then upon completion
%displays the final time
function status = printOut(t, y, flag)
   if strcmp(flag, 'done') %If done itertating
       disp(tPrev); %Display final time
   else
       tPrev = t(end); %Store latest time;
   end
   status = 0; %Continue running
end   
options = odeset('Events', @events,'OutputFcn', @printOut); %Sets ode options
t = linspace(0, 100, 100000); %Creating timepoint to integrate over
    
%This is the differential equation for the equation of motion in terms of
%theta provided by the Kasdin and Paley textbook
    function do = tetherBallDiffEq(t, o)
        do = [o(2);...
            r*o(2)^2/(l-r*o(1))];
    end
%Solves and stores the resultion in res1
[~, res1] = ode45(@tetherBallDiffEq, t, init, options);

%This is the function for calculating the remaining length of tether after
%a certain amoung of angular displacement
    function d = tetherLength(l, r, theta)
        d = l - r*theta;
    end

%Plotting both the circle and the trajectory of the tetherball:
hold on;
count = 1;
for k = 0:pi/5000:2*pi
    xC(count) = cos(k); %Plots the circle
    yC(count) = sin(k);
    count = count+1;
end

thetaVec = res1(:, 1);                  %Creates vectors to plot in polar
rhoVec = tetherLength(l, r, thetaVec);
%Plotting:
f1 = figure(1);
plot(xC, yC, 'k');
polar(thetaVec, rhoVec);
set(f1, 'Position', [200, 1, 800, 800]);
xlabel('x (m)', 'Fontsize', 16);
ylabel('y (m)','Fontsize', 16);
title('Tetherball Trajectory from Example 4.7', 'Fontsize', 20);
str1 = 'Ball Hits Post at ';
str2 = num2str(tPrev);
str3 = strcat(str1, str2, ' s');
text(-5, 0, str3);


end