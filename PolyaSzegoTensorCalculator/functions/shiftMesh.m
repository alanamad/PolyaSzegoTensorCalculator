%**************************************************************************
% MATLAB implemenetation to shift the object to the origem
%**************************************************************************

function p_new = shiftMesh(p, x_c, y_c)
    p_new(1,:) = p(1,:) -x_c;
    p_new(2,:) = p(2,:) -y_c;
end
