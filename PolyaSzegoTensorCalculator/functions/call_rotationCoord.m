%**************************************************************************
% MATLAB implemenetation of the Rotation Matrix of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Compute Rotation Matrix
%
% HISTORY
% A. Amad       11/2018: code implementation
%**************************************************************************

function [Xnew, Ynew] = call_rotationCoord(X, Y, theta)

V = [X; Y];

theta = theta*pi/180;

        R = [cos(theta) -sin(theta); ...
             sin(theta)  cos(theta)];

Vnew = R*V;

Xnew = Vnew(1);
Ynew = Vnew(2);

end
