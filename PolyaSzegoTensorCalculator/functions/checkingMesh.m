%**************************************************************************
% MATLAB implementation of the checking mesh
%**************************************************************************
%  
% DESCRIPTION
% This program check the quality of the mesh
%
% HISTORY
% A. Amad       06/2019: code implementation
%**************************************************************************


function checkingMesh(p,t)
    %% Checking mesh
    Nnodes = size(p,2);
    x = p(1,:);
    y = p(2,:);

    if length(unique(t))~= Nnodes
        disp(' ')
        disp('Wrong discretisation!')
        disp(' ')
        opts = struct('WindowStyle','modal', 'Interpreter','tex');
        errordlg('More vertices than triangle!','Discretisation Mesh Error',opts)
        return
    end

    det = (x(t(2,:)) - x(t(1,:))) .* (y(t(3,:)) - y(t(1,:))) - (x(t(3,:)) - x(t(1,:))) .* ...
          (y(t(2,:)) - y(t(1,:)));
    if any(det < 0)
        element = find(det < 0);
        disp(' ')
        disp('Element number(s) with negative determinant')
        disp(' ')
        disp(element)
        disp(' ')
        opts = struct('WindowStyle','modal', 'Interpreter','tex');
        errordlg('Element mesh with negative determinant!','Element Mesh Error',opts)
        return
    end
end