function [res1, res2] = cxa2_4()
%This function implements the equations of motion for the simple pendulum
%and then performs the ode45 on them to derive the trajectory of the
%pendulum in both polar and cartesian coordinates. It also analyzes the
%fault in results that occurs when recalculating the length of the pendulum
%arm at each step

g = 9.81; % In m/s
r = 1; %In meters
t = linspace(0, 10, 500); % In seconds
thetaInit = 0.1; %In radians
thetaDotInit = 0; %In radians per second
lInit = 1; % In meters

%Calculating the corresponding initial conditions in the Cartesian
%coordinates from the privided initial conditions in polar coordinates
xInit = lInit/sqrt((tan(thetaInit))^2+1);
yInit = sqrt(lInit^2-xInit^2);

%Creating the initial condition vectors from the ode45 functions from the
%provided and calculated initial conditions
initCart = [xInit 0 yInit 0];
initPol = [thetaInit, thetaDotInit];

    %Function that calculates the length of the pendulum arm from a given x
    %and y coordinate in the Cartesian coordinate frame whose origin lies
    %at the origin of the pendulum arm
    function l = pendulumLength(x, y)
     l = sqrt(x.^2+y.^2); %Pythagorean theorem calculation
    end
    
    %Function that calculates the corresponding theta value of the pendulum
    %given the x and y coordinates of the pendulum mass
    function theta = pendulumAngle(x, y)
     theta = atan2(y,x); %Uses the built-in atan2 function to find the more
                         %precise angle without running into quadrant error
    end
    
    %Function found using the derivative of the above arctan function that
    %calculates the corresponding angular velocity of the mass given the x
    %and y coordinates in a cartesian plane as well as the derivatives of
    %those coordinates at the corresponding time points
    function thDot = thetaDot(x, y, xDot, yDot)
        thDot = (1/(1+(y/x)^2))*((yDot*x-xDot*y)/x^2);
    end

    %Sets up the first-order form for the equations of motion for the
    %pendulum in cartesion coordinates provided by the Kasdin and Paley
    %textbook w(1), w(2), w(3), and w(4) represent x, xDot, y, and yDot
    %respectively
    function dw = simplePendulumCart (t, w)
    dw = [w(2);...
    (-w(1)*w(2)^2-w(1)*w(4)^2+w(1)*w(3)*g)/pendulumLength(w(1), w(3))^2;...
    w(4);...
    (-w(3)*w(2)^2-w(3)*w(4)^2-w(1)^2*g)/pendulumLength(w(1), w(3))^2];
    end

    %Sets up the first-order form for the equations of motion for the
    %pendulum in polar coordinates as I derived to fit the conditions of
    %the new theta taken not from the rest position of the pendulum but
    %instead from an axis parallel to the ground with clockwise rotations
    %being negative and counterclockwise rotations being positive
    function dw = simplePendulumPol (t, w)
        dw = [w(2);...
            -g*cos(w(1))/r];
    end
    
    %Solving the differential equations using ode45
   [~, res1] = ode45(@simplePendulumCart, t, initCart);
   [~, res2] = ode45(@simplePendulumPol, t, initPol);
    
   %Processes the resulting x and y values and assigns the corresponding
   %that value
   thetaVec = zeros(1, length(res1)); %Preassigning space for efficiency
   for n = 1:length(res1) %Iterating over all the results
       theta = pendulumAngle(res1(n,1), res1(n, 3)); %Finding the angle
       if theta > pi/2 %Correcting for being out of atan2's range
        theta = theta - 2*pi; %Moving the theta value back in range
       end
       thetaVec(n) = theta; %Storing the value in the thetaVec
   end
   
   %Processes the x, y, xDot, and yDot values and provides the
   %corresponding thetaDot values
   thDotVec = zeros(1, length(res1)); %Preassigning for efficiency
   for n = 1:length(res1)
       %Uses thetaDot to find the time derivative of theta using the input
       %values at each time point
       thDot = thetaDot(res1(n,1), res1(n,3), res1(n, 2), res1(n, 4));
       thDotVec(n) = thDot; %Stores values in thDotVec
   end
   
   
   %Plotting the required figures
   
   f1 = figure(1);
   plot(t, thetaVec, t, res2(:, 1), 'or');
   set(f1, 'Position', [1,1,1200, 800]);
   title('Pendulum Angle (Theta) as a Function of Time [Cartesian vs Polar]', 'Fontsize', 20);
   xlabel('Time (s)', 'Fontsize', 14);
   ylabel('Pendulum Angle (radians)', 'Fontsize', 14);
   legend('Cartesian', 'Polar');
   
   f2 = figure(2);
   plot(t, thDotVec, t, res2(:, 2), 'or');
   set(f2, 'Position', [1,1,1200, 800]);
   title('Pendulum Angle Rate (ThetaDot) as a Function of Time [Cartesian vs Polar]', 'Fontsize', 20);
   xlabel('Time (s)', 'Fontsize', 14);
   ylabel('Pendulum Angle Rate (radians/second)', 'Fontsize', 14);
   legend('Cartesian', 'Polar');
   
   f3 = figure(3);
   plot(t, pendulumLength(res1(:, 1), res1(:, 3)));
   set(f3, 'Position', [1,1,1200, 800]);
   title('Pendulum Length as a Function of Time', 'Fontsize', 20);
   xlabel('Time (s)', 'Fontsize', 14);
   ylabel('Pendulum Length (meters)', 'Fontsize', 14);
   
end

