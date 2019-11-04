%**************************************************************************
% MATLAB implementation of the computation of the Polya-Szego tensor 
%**************************************************************************
%  
% DESCRIPTION
% Main program used to compute the Polya-Szego tensor
% Method options - Layer Potential/Boundary Integral
%                  Finite Element/Interpolated BEM
% HISTORY
% A. Amad       10/2018: code implementation
% A. Amad       05/2019: code updating
%**************************************************************************

clc; close all;
format long

%=========================================================================
% Input:
% methodOption  : Layer Potential/Boundary Integral/FEM/Interpolated BEM
% fileName      : Geometry file name
% numberRef     : Number of refinement
% plot          : Figures, tensor as ellipse, checking the object
% dataObject    : Size, contrast, rotation angle
% 
% Output:
% PStensor      : Polya-Szego tensor
%=========================================================================

%-------------------------------------------------------------------------
% Method option
%-------------------------------------------------------------------------
methodOption = 1;  % 1- Layer Potential  (LP)
                   % 2- Boundary Integral (BI)
                   % 3- Finite Element Method (FEM)
                   % 4- Interpolated BEM (IntBEM)

%-------------------------------------------------------------------------
%  Object file name
%-------------------------------------------------------------------------
% Example of object geometry available 
%          fileName = 'circ_r05';
%          fileName = 'ellipse_a2_b1';
%          fileName = 'triangle';
%          fileName = 'Pentagon_NonConvex';

fileName = 'ellipse_a2_b1';

%-------------------------------------------------------------------------
% Object type
% IMPORTANT! If the object is ellipsoidal, the radii must be included
%-------------------------------------------------------------------------
isEllipsoidal = 'No';  % Yes or No

radius_a = 2;       % Circle: a = b 
radius_b = 1;       % Ellipse: a =~ b

%-------------------------------------------------------------------------
% Mesh Refinement Options
%-------------------------------------------------------------------------
refType = 1;  % 1- Uniform refinement 
              % 2- Non-uniform refinement
              % 3- Local refinement

%-------------------------------------------------------------------------
% Local refinement 
%-------------------------------------------------------------------------
vertex=[ ]; % vertices to be refined
r=0.1; % initial radius

%-------------------------------------------------------------------------
% Number of refinement 
%-------------------------------------------------------------------------
numberRefinement = 3;    % Number of refinement 

%-------------------------------------------------------------------------
% Plot OptionsYes
%-------------------------------------------------------------------------
plotFigures = 'Yes'; % Yes or No

plotTensorAsEllipse = 'No';  % Yes or No

plotObject = 'Yes';  % Yes or No

%=========================================================================
%% Input Data
% Size of object
alpha =0.01;

% Contrast
mur = 1.5;

% Rotation angle (degrees)6
theta = 0;
%=========================================================================

%-------------------------------------------------------------------------
% Checking Mesh
%-------------------------------------------------------------------------
checkMesh = 'No'; % Yes or No

%% Compute data

data.methodOption = methodOption;
data.refType = refType;
data.numberRef = numberRefinement;
data.plotFigures = plotFigures;
data.plotTensorAsEllipse = plotTensorAsEllipse;
data.plotObject = plotObject;
data.alpha = alpha;
data.mur = mur;
data.theta = theta;
data.isEllipsoidal = isEllipsoidal;
data.fileName = fileName;
data.checkMesh = checkMesh;
data.vertex = vertex;
data.r = r;

if isequal(isEllipsoidal,'Yes')
    data.radius_a = radius_a;
    data.radius_b = radius_b;
end

cd('functions')
    PStensor = computePS(data);
cd ..

%% Output
disp(['PStensor = ', mat2str(PStensor) ]) ;

