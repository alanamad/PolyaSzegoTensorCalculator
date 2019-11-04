%**************************************************************************
% MATLAB implemenetation to find the centroid of the object
%**************************************************************************
%  
% DESCRIPTION
% Find the centroid of the object
%
% Input Arguments:
%   -p: vector of coordinates
%   -e: connectivity of object edges
%
% Output Arguments:
%   -x_c: x location of centroid
%   -y_c: y location of centroid
%   -area: area of polygon
%
% HISTORY
% A. Amad       05/2019: code implementation
%**************************************************************************

function [x_c, y_c, area] = centroidPolygon(e,x,y)

    Npanels = length(e(1,:));  % Number of edges

    coordX = [];
    coordY = [];
    for i=1: Npanels
        coordX = [coordX, x(e(1,i)) ]; %#ok<*AGROW>
        coordY = [coordY, y(e(1,i)) ];
    end
    coordX = [coordX, x(e(1,1)) ];
    coordY = [coordY, y(e(1,1)) ];
    x = coordX;
    y = coordY;
    
    det = x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1);
    area = sum(det)/2;
    x_c = (sum((x(2:end)+x(1:end-1)).*det)*1/6)/area;
    y_c = (sum((y(2:end)+y(1:end-1)).*det)*1/6)/area;

    x = x -x_c;
    y = y -y_c;

    disp(['Centroid object: (x_c, y_c)   = (' num2str(x_c),',',num2str(y_c),')' ]);
    
    disp(' ');

    disp('Centring the object at the origem ');

    disp(' ');

    det = x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1);
    area = sum(det)/2;
    x_c2 = (sum((x(2:end)+x(1:end-1)).*det)*1/6)/area;
    y_c2 = (sum((y(2:end)+y(1:end-1)).*det)*1/6)/area;

    disp(['        Object centred at (x_c, y_c)   = (' num2str(x_c2),',',num2str(y_c2),')' ]);

end
