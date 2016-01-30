%This script plots three points (O, P, and Q) in the three dimensional
%space and then plots two vectors (r_P/O, and r_Q/P) between those points.

% Creates the coordinate vectors for the three points
xCoord = [0 4 -3]; 
yCoord = [0 5 4];
zCoord = [0 6 -2];

hold on;

%Labels the axes
xlabel('X', 'FontSize', 14, 'Color', 'b');
ylabel('Y', 'FontSize', 14, 'Color', 'b');
zlabel('Z', 'FontSize', 14, 'Color', 'b');

%Plots the points and then the two vectors
plot3(xCoord, yCoord, zCoord, 'ok');
plot3(xCoord(1:2), yCoord(1:2), zCoord(1:2), 'k');
plot3(xCoord(2:3), yCoord(2:3), zCoord(2:3), 'k');

% Labels the vecots and then the three points
h1 = text((xCoord(1)+xCoord(2))/2 + 1, (yCoord(1)+yCoord(2))/2, ...
    (zCoord(1)+zCoord(2))/2, 'r_{P/O}');
h2 = text((xCoord(2)+xCoord(3))/2 - 2, (yCoord(2)+yCoord(3))/2, ...
    (zCoord(2)+zCoord(3))/2, 'r_{Q/P}');
h3 = text(xCoord(1) + 1, yCoord(1), zCoord(1), 'O');
h4 = text(xCoord(2) + 1, yCoord(2), zCoord(2), 'P');
h5 = text(xCoord(3) + 1, yCoord(3), zCoord(3), 'Q');

%Sets the font size of the labels
set(h1, 'FontSize', 14);
set(h2, 'FontSize', 14);
set(h3, 'FontSize', 14);
set(h4, 'FontSize', 14);
set(h5, 'FontSize', 14);

view(3);
