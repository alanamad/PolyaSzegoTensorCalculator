%**************************************************************************
% MATLAB implemenetation of the computation of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Main program used to run compute the PS tensor using BEM
% Method options - Layer Potential/Boundary Integral
%                  Finite Element/Interpolated BEM
%
% HISTORY
% A. Amad       10/2018: code implementation
% A. Amad       05/2019: code updating
%**************************************************************************

function PStensor = computePS(data)

method = data.methodOption;
meshOption = data.refType;

switch meshOption
    case 1      % Regular refinement 
        meshOption = 'Reg';
    case 2      % Irregular refinement 
        meshOption = 'Irreg';
    case 3      % Local refinement
        meshOption = 'Local';
end
data.meshOption = meshOption;

cd ..
    cd('problemfiles')
        [mesh, params] = objectData(data);
cd ..


tic 
cd('functions')
    if isequal(method,1) || isequal(method,2)
        switch method
            case 1      % Layer Potential
                methodOption = 'LP';
            case 2      % Boundary Integral
                methodOption = 'BI';
        end
        data.methodOption = methodOption;
        PStensor = main(data, mesh, params);
    elseif isequal(method,3) || isequal(method,4)
        switch method
            case 3      % Boundary Integral
                methodOption = 'FEM';
            case 4      % Boundary Integral
                methodOption = 'IntB';
        end
        data.methodOption = methodOption;
        PStensor = mainFEM_IB(data, mesh, params);
    end
toc
end

