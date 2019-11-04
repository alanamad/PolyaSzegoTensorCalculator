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

function Mr = call_rotationMatrix(M, theta)

theta = theta*pi/180;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Mr = zeros(2,2);

for m=1:2
    for n=1:2
        for i=1:2
            for j=1:2
                Mr(m,n) = Mr(m,n) + R(m,i)*R(n,j)*M(i,j);
            end
        end
    end
end

end
