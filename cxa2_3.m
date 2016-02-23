function [] = cxa2_3()
% Function that analyses the csv data for the spring-mass damper system and
% then calculates the spring constant (k) and damping coefficient (b)
% values from each ringdown of the data

tolerance = 0.1; %Chosen tolerance for smoothening
springMass = 0.026; %in kg
mass = 0.330; %attached mass in kg
systemMass = mass + (1/3)*springMass;   %in kg
zeroPoint = 0.52; %zero point voltage

%Inputs the ringdown data
data1 = csvread('../Other/data/ringdown1.csv');
data2 = csvread('../Other/data/ringdown2.csv');


%Eliminates all points that have a time value of 0
sz1 = size(data1);
sz2 = size(data2);
rowNum1 = sz1(1);
rowNum2 = sz2(1);
counter1 = 1;
counter2 = 1;
for n=1:rowNum1
    
    if data1(n, 1) ~= 0
        newDat1(counter1, 1) = data1(n, 1);
        newDat1(counter1, 2) = data1(n, 2);
        counter1 = counter1 + 1;
    end
 
end

for n = 1:rowNum2
    if data2(n, 1) ~= 0
        newDat2(counter2, 1) = data2(n, 1);
        newDat2(counter2, 2) = data2(n, 2);
        counter2 = counter2 + 1;
    end
end

%Eliminates all points that don't mark a change higher than the input
%tolerance, so as to filter the data
sz1 = size(newDat1);
sz2 = size(newDat2);
rowNum1 = sz1(1);
rowNum2 = sz2(1);
counter1 = 2;
counter2 = 2;

smoothDat1(1, 1) = (newDat1(1, 1));
smoothDat1(1, 2) = (newDat1(1, 2));
for n = 2:rowNum1 
    if (newDat1(n, 2)-newDat1(n-1, 2)) > tolerance
        smoothDat1(counter1, 1) = (newDat1(n, 1));
        smoothDat1(counter1, 2) = (newDat1(n, 2));
        counter1 = counter1 + 1;
    end
end

smoothDat2(1, 1) = (newDat2(1, 1));
smoothDat2(1, 2) = (newDat2(1, 2));
for n = 2:rowNum2 
    if (newDat2(n, 2)-newDat2(n-1, 2)) > tolerance
        smoothDat2(counter2, 1) = (newDat2(n, 1));
        smoothDat2(counter2, 2) = (newDat2(n, 2));
        counter2 = counter2 + 1;
    end
end

%Smoothens the data using MATLABs built in function
smoothDat1(:,1) = smooth(smoothDat1(:,1));
smoothDat1(:,2) = smooth(smoothDat1(:,2));

%Saving this value for later use in initial condition solving
saveVal = smoothDat1(10, 2);

%Eliminating unnescessary data vectors
clear('counter1', 'counter2', 'data1', 'data2', 'newDat2',...
    'rowNum1', 'rowNum2');

%Eliminates all occurences prior to the release of the mass, this will help
%when determining the maxima of the oscillations to find the logarithmic
%decrement
index = 1;
t = smoothDat1(index, 1);
while t <= 0 
    index = index+1;
    t = smoothDat1(index, 1);
end
smoothDat1 = smoothDat1(index:end, :);



index = 1;
t = smoothDat2(index, 1);
while t <= 0 
    index = index+1;
    t = smoothDat2(index, 1);
end
smoothDat2 = smoothDat2(index:end, :);

%Finds the local peaks of the data
[pks1, locs1] = findpeaks(smoothDat1(:, 2));
[pks2, locs2] = findpeaks(smoothDat2(:, 2));

%Searches through the local peaks in order to find the actual oscillation
%maxima for two consecutive crests

hasSwitched = false;
maxInd1 = zeros(1,2);
maxVals1 = zeros(1,2);
maxTimes1 = zeros(1,2);

[maxVals1(1), maxInd1(1)] = max(pks1); 
maxTimes1(1) = smoothDat1(locs1(maxInd1(1)),1);

n = maxInd1(1);
while ~hasSwitched 
    n = n + 1;
    if pks1(n) < 0.52 
     hasSwitched = true;
    end
end
pks1 = pks1(n+1:end);

[maxVals1(2), maxInd1(2)] = max(pks1);
maxTimes1(2) = smoothDat1(locs1(maxInd1(2)+n),1);

hasSwitched = false;
maxInd2 = zeros(1,2);
maxVals2 = zeros(1,2);
maxTimes2 = zeros(1,2);

[maxVals2(1), maxInd2(1)] = max(pks2); 
maxTimes2(1) = smoothDat2(locs2(maxInd2(1)),1);

n = maxInd2(1);
while ~hasSwitched 
    n = n + 1;
    if pks2(n) < 0.52 
     hasSwitched = true;
    end
end
pks2 = pks2(n+1:end);

[maxVals2(2), maxInd2(2)] = max(pks2);
maxTimes2(2) = smoothDat2(locs2(maxInd2(2)+n),1);

%Cleanup
clear('hasSwitched', 'index', 'locs1', 'locs2', 'maxInd1', 'maxInd2',...
    'n', 'pks1', 'pks2', 'sz1', 'sz2');

%Calculates the damping coefficient using the logarithmic decrement
logDec = zeros(1, 2);
TauD = zeros(1, 2);
b = zeros(1, 2);
logDec(1) = log((maxVals1(1)-zeroPoint)/(maxVals1(2)-zeroPoint));
logDec(2) = log((maxVals2(1)-zeroPoint)/(maxVals2(2)-zeroPoint));
TauD(1) = maxTimes1(2) - maxTimes1(1);
TauD(2) = maxTimes2(2) - maxTimes2(1);
b(1) = 2*systemMass*logDec(1)/TauD(1);
b(2) = 2*systemMass*logDec(2)/TauD(2);
fprintf('\nThe damping coefficient for the first ringdown is %s N/m/s\n\n', b(1));
fprintf('The damping coefficient for the second ringdown is %s N/m/s\n', b(2));

%Calculates the spring constant
k = zeros(1,2);
k(1) = b(1)^2*(1+(2*pi/logDec(1))^2)/(4*systemMass);
k(2) = b(2)^2*(1+(2*pi/logDec(2))^2)/(4*systemMass);
fprintf('\nThe spring constant for the first ringdown is %s N/m\n\n', k(1));
fprintf('The spring constant for the second ringdown is %s N/m\n', k(2));

%Creating the function found analytically for a dampled, driven SHO
function y = forcedDampedSHO(x, c1, c2, k, m, b)
    %To account for the discrepancy between the start of the oscillation
    %and the actual start of the timer
    x = x+0.065;  
    om0 = sqrt(k/m); %Calculating Omega0
    zeta = b/(2*m*om0); %Calculating zeta
    omD = om0*sqrt(1-zeta^2); %Calculating Omega Driving
    y = exp(-zeta*om0*x)*(c1*cos(omD*x)+c2*sin(omD*x))+zeroPoint;
end

%Solving for coefficients through the initial conditions
c1 = saveVal-zeroPoint; %First coefficient
om0 = sqrt(k(1)/systemMass); % Omega naught
zeta = b(1)/(2*systemMass*om0); 
omD = om0*sqrt(1-zeta^2); % Omega drive
c2 = zeta*om0*c1/omD; %Second coefficient

hold on
plot(newDat1(:, 1), newDat1(:, 2), 'r');
for j = newDat1(1,1):0.0005:newDat1(end, 1)
    plot(j, forcedDampedSHO(j, c1, c2, k(1), systemMass, b(1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that analyses the csv data for the louderspeaker system and
% then calculates the spring constant (k) and damping coefficient (b)
% values and the m value of the loudspeaker

deltaM = 0.0093;

data3 = csvread('data', filesep, 'loudspeaker_m.csv');
data4 = csvread('data', filesep, 'loudspeaker_m_plus_dm.csv');

freqM = data3(:, 1);
freqMPlus = data4(:, 1);
datM = data3(:, 2);
datMPlus = data4(:, 2);

%Finding the maximum amplitude for each of the systems in order to find the
%correspoding resonance frequency.
dataInd = zeros(1, 2);
dataMax = zeros(1, 2);
freqVec = zeros(1, 2);

[dataMax(1), dataInd(1)] = max(datM);
[dataMax(2), dataInd(2)] = max(datMPlus);

freqVec(1) = freqM(dataInd(1));
freqVec(2) = freqMPlus(dataInd(2));

freqVec = freqVec * 2 * pi;
%Calculating the original mass of the loudspeaker
finalMass = deltaM/(freqVec(1)/freqVec(2)-1);
fprintf('\nThe original mass of the loudspeaker is  %s kg\n\n', finalMass);
%Calculating the damping coefficient for the loudspeaker
bLS = finalMass*freqVec(1)/(dataMax(1)-4.5);
fprintf('The damping coefficient of the loudspeaker is  %s N/m/s\n\n', bLS);
%Calculating the spring constant for the loudspeaker
kLS = (bLS^2/finalMass)*dataMax(1)^2;
fprintf('The spring constant of the loudspeaker is  %s N/m\n\n', kLS);


end