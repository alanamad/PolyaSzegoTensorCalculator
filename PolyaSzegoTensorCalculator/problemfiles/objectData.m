%**************************************************************************
% MATLAB implemenetation of the computation of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Object data
%
% HISTORY
% A. Amad       10/2018: code implementation
%**************************************************************************

function [mesh, params] = objectData(data)

sizeObject = data.alpha;
contrast = data.mur;
theta = data.theta;
isEllipsoidal = data.isEllipsoidal;
method = data.methodOption;
fileName = data.fileName;
meshType = data.meshOption;


%% load geometry
% Geometry generated by pdetool and exported the decomposed geometry g
cd ..
cd('geometry')
    load(sprintf('%s.mat', fileName));  %#ok<LOAD>
cd ..
cd('problemfiles')


rotationAngle = theta;
params.rotationAngle = rotationAngle;

if isequal(isEllipsoidal,'Yes')
    radius_a = data.radius_a;
    radius_b = data.radius_b;
    name_example = sprintf('%s_rot%d_radius_a%d_b%d', fileName, rotationAngle,radius_a,radius_b);
    params.name_example = name_example;
    params.name_obj = 'ellipsoidal';
    params.radius_a = radius_a;
    params.radius_b = radius_b;
else
    name_example = sprintf('%s_rot%d', fileName, rotationAngle);
    params.name_example = name_example;
    params.name_obj = name_example;
end




%% parameters 

% permeability contrast
mur = contrast;
params.mur = mur;

% object size
params.alpha = sizeObject;

% Constant lambda
params.lamb = .5 * (mur+1) / (mur - 1);


%% mesh generation

% shift the object
[p,e,~] = initmesh(g,'hmax',25);

x = p(1,:);
y = p(2,:);

if isequal(method,1) || isequal(method,2)
    [x_c, y_c, ~] = centroidPolygon(e,x,y);

    if isequal(meshType,'Irreg')
        [p,e,t] = initmesh(g,'hmax',0); % irregular mesh
        [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
        [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
        [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh

        e_aux = unique(e(1:2,:));
        it = [];
        for ii = 1: size(e_aux,1)
            for jj = 1:3
                it = [it; find(e_aux(ii)==t(jj,:))'];
            end
        end
        it = unique(it);
        u=zeros(size(p,1));
        [p,e,t] = refinemesh(g,p,e,t,u, it,'regular'); % refinement mesh

    elseif isequal(meshType,'Reg')
        [p,e,t] = initmesh(g,'hmax',25); % regular mesh
        for i=1:2
            [p,e,t] = refinemesh(g,p,e,t); % refinement mesh
        end
        e_aux = unique(e(1:2,:));
        it = [];
        for ii = 1: size(e_aux,1)
            for jj = 1:3
                it = [it; find(e_aux(ii)==t(jj,:))'];
            end
        end
        it = unique(it);
        u=zeros(size(p,1));
        [p,e,t] = refinemesh(g,p,e,t,u, it,'regular'); % refinement mesh
    else
        [p,e,t] = initmesh(g,'hmax',25); % regular mesh
        for i=1:3
            [p,e,t] = refinemesh(g,p,e,t); % refinement mesh
        end
    %     [p,e,t] = initmesh(g,'hmax',0); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh
    %     [p,e,t] = refinemesh(g,p,e,t,'longest'); % irregular mesh

        for i=1:2
            e_aux = unique(e(1:2,:));
            it = [];
            for ii = 1: size(e_aux,1)
                for jj = 1:3
                    it = [it; find(e_aux(ii)==t(jj,:))'];
                end
            end
            it = unique(it);
            u=zeros(size(p,1));
            [p,e,t] = refinemesh(g,p,e,t,u, it,'regular'); % refinement mesh
        end
    end
    
elseif isequal(method,3) || isequal(method,4)
    
    connect_bound_obj = e(6,:)==2;
    e_obj = e(:,connect_bound_obj);

    [x_c, y_c, ~] = centroidPolygon(e_obj,x,y);
    
        if isequal(meshType,'Irreg')
            [p,e,t] = initmesh(g,'hmax',0); % irregular mesh
            % irregular mesh
            [p,e,t] = refinemesh(g,p,e,t,'longest');
            [p,e,t] = refinemesh(g,p,e,t,'longest');

        elseif isequal(meshType,'Reg')
            % regular mesh
            [p,e,t] = initmesh(g,'hmax',25); % regular mesh
            for i=1:2
                [p,e,t] = refinemesh(g,p,e,t); % refinement mesh
            end
        elseif isequal(meshType,'Local')
            % regular mesh
            [p,e,t] = initmesh(g,'hmax',25); % regular mesh
            for i=1:4
                [p,e,t] = refinemesh(g,p,e,t); % refinement mesh
            end
        end

end

if isequal(meshType,'Local')
    % local refinement
    tagvertex = data.vertex;
    if ~isempty(tagvertex)
        r = data.r;

        xc = pdeintrp(p,t,p(1,:)');
        yc = pdeintrp(p,t,p(2,:)');
        [ntag, ~]=size(tagvertex);
        it=[];
        for i=1:ntag
           it=[it;find(sqrt((tagvertex(i,1)-xc).^2+(tagvertex(i,2)-yc).^2)< r)'];%#ok<*AGROW>
        end
        it= unique(it);
        u=zeros(size(p,1));
        [p,e,t] = refinemesh(g,p,e,t,u, it,'regular'); % refinement mesh
    end
end


p0 = p;
e0 = e;
t0 = t;

p = shiftMesh(p, x_c, y_c);

if rotationAngle ~= 0 
    np = size(p,2);
    for i = 1:np
        [p(1,i), p(2,i)] = call_rotationCoord(p(1,i), p(2,i), rotationAngle);
    end
end

mesh.p0 = p0;
mesh.t0 = t0;
mesh.e0 = e0;
mesh.p = p;
mesh.e = e;
mesh.t = t;
mesh.g = g;
mesh.x_c = x_c;
mesh.y_c = y_c;


end