%**************************************************************************
% MATLAB implemenetation of the computation of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Compute the PS tensor using BEM - Layer Potential/Boundary Integral
%
% HISTORY
% A. Amad       10/2018: code implementation
%**************************************************************************

function PStensor = main(data, mesh, params)

%% parameters 

% input data
method = data.methodOption;
nRefinament = data.numberRef;
plotFigures = data.plotFigures;
plotTensor = data.plotTensorAsEllipse;
plotObject = data.plotObject;
meshType = data.meshOption;

% mesh parameter
p = mesh.p; e = mesh.e; t = mesh.t; g = mesh.g;
rotationAngle = params.rotationAngle;

% permeability contrast
mur = params.mur;

% object size
alpha = params.alpha;

lamb = params.lamb;


%% Convergence rate
PS = [];
ndof = [];

m11 = [];
m22 = [];
m12 = [];
m21 = [];

Mex11 = [];
Mex12 = [];
Mex21 = [];
Mex22 = [];

error = [];
errorM12 = [];
errorM21 = [];
errorComparison  = [];
errorM12M21 = [];
relError = [ ];

for nsteps = 0:nRefinament
    disp('----------------');
    disp(['Refinement = ' num2str(nsteps) ]) ;
    
    % mesh data
    x = p(1,:);
    y = p(2,:);

    %% Checking mesh

    checkMesh = data.checkMesh;
    if isequal(checkMesh,'Yes')
        checkingMesh(p,t)
    end

    %% Compute the adjoint double layer potential and the linear system 

    Npanels = size(e,2); % Number of edges
    disp(['Npanels = ' num2str(Npanels) ]) ;
    
    Lamb = lamb * eye(Npanels);
    xm=0.5.*(x(e(1,:))+x(e(2,:))); % x-middle point edge
    ym=0.5.*(y(e(1,:))+y(e(2,:))); % y-middle point edge
    lg=sqrt((x(e(2,:))-x(e(1,:))).^2+(y(e(2,:))-y(e(1,:))).^2); % length edge
    nx=(y(e(2,:))-y(e(1,:)))./lg; % normal x
    ny=(-x(e(2,:))+x(e(1,:)))./lg; % normal y

    if isequal(plotObject,'Yes')
        figure(1)
        set(1,'WindowStyle','docked');
        quiver(xm,ym,nx,ny)
        hold on; pdemesh(p,e,t)
        set(gca,'XColor','k','YColor','k','LineWidth',2); 
        set(gca,'XTick',[],'YTick',[],'Ztick',[]);
        box on;
        cd ..
        cd('results')
            filename = sprintf('Object_step%d', nsteps);
            print(1, '-dtiff', filename);
        cd ..
        cd('functions')
    end
    
    % compute the adjoint double layer potential
    Kstar = zeros(Npanels,Npanels);  % adjoint double layer potential

    for i=1:Npanels
        for j=1:Npanels
            if i==j
                Kstar(i,j)=  0.0;
            else
                dist = sqrt((xm(i)-xm(j))^2 + (ym(i)-ym(j))^2);
                dotxy = (xm(i)-xm(j))*nx(i) + (ym(i)-ym(j))*ny(i);
                Kstar(i,j)= lg(j)*dotxy/(2*pi*dist^2);  % adjoint double layer potential
            end
        end
    end
        
    % compute the linear system
    A= Lamb -Kstar; % lambda*I - K^star

    if isequal(method,'LP')
        normal = [nx; ny];
        phi = A\normal'; % vector solution phi
    elseif isequal(method,'BI')
        rhs_x= ((1/(mur-1))*(-1/2 * eye(Npanels) + Kstar)*nx');
        rhs_y= ((1/(mur-1))*(-1/2 * eye(Npanels) + Kstar)*ny');
        rhs = [rhs_x rhs_y];
        Grad_phi_n = A\rhs; % vector solution phi
    end

    %% Compute the Polya-Szego tensor
    Y =[xm; ym];

    ndof = [ndof; length(Y)]; %#ok<*AGROW>
    if isequal(method,'LP')
        m =zeros(2,2);
        for i=1:2
            for j=1:2
                m(i,j) = alpha^2*dot(phi(:,i),lg.*Y(j,:));
            end
        end
    elseif isequal(method,'BI')
        moment =zeros(2,2);
        area = sum(pdetrg(p,t));
        for i=1:2
            for j=1:2
                moment(i,j) = alpha^2*(mur-1)^2 *dot(Grad_phi_n(:,i),lg.*Y(j,:));
            end
        end
        B = alpha^2*(mur-1)*area* eye(2);
        m = B + moment;
    end
    PS =[PS; Npanels m(1,1) m(1,2) m(2,1) m(2,2)];
    m11 =[m11; m(1,1)];
    m12 =[m12; m(1,2)];
    m21 =[m21; m(2,1)];
    m22 =[m22; m(2,2)];


    disp(['PS = ' mat2str(m) ]) ;

    name_obj = params.name_obj;
    if isequal(name_obj,'ellipsoidal')
        a = params.radius_a; % radius a
        b = params.radius_b; % radius b


        % Exact tensor: Circle or Ellipse
        M_exact = alpha^2*(mur-1)*a*b*pi*[(a+b)/(a+mur*b) 0; 0 (a+b)/(b+mur*a)];

        M_exact = call_rotationMatrix(M_exact, rotationAngle);
        disp(['M_exact = ' mat2str(M_exact) ]) ;


        err = m - M_exact;
        Error = norm(err,'fro');
        error = [error; Error];
        
        er = norm(M_exact,'fro');
        relErr = Error/er;
        relError = [relError; relErr];
        plotData.relError = relError;

        Mex11 = [Mex11; M_exact(1,1)];
        Mex12 = [Mex12; M_exact(1,2)];
        Mex21 = [Mex21; M_exact(2,1)];
        Mex22 = [Mex22; M_exact(2,2)];

        err = m(1,2) - M_exact(1,2);
        Error = norm(err);
        errorM12 = [errorM12; Error];

        err = m(2,1) - M_exact(2,1);
        Error = norm(err);
        errorM21 = [errorM21; Error];
        
        err = m(1,2) - m(2,1);
        Error = norm(err);
        errorM12M21 = [errorM12M21; Error];

        plotData.Mex11 = Mex11;
        plotData.Mex12 = Mex12;
        plotData.Mex21 = Mex21;
        plotData.Mex22 = Mex22;
        plotData.M_exact = M_exact;
        

        plotData.errorM12 = errorM12;
        plotData.errorM21 = errorM21;
        plotData.errorM12M21 = errorM12M21;

    else
        err = m(1,2) - m(2,1);
        Error = norm(err,'fro');
        error = [error; Error];
        errorComparison = [errorComparison; Npanels Error];

    end

    plotData.error = error;

    plotData.m11 = m11;
    plotData.m12 = m12;
    plotData.m21 = m21;
    plotData.m22 = m22;

    %% Plots
    name = params.name_example;
    name_example = sprintf('%s_reftype_%s', name, meshType);

    if isequal(plotFigures,'Yes')
        call_plot(method, ndof, params, data, plotData)
    end

    if and(isequal(name_obj,'ellipsoidal'), isequal(plotTensor,'Yes'))
        plotEllipseSolution(method, mur, m,M_exact, alpha, params, Npanels)
    end


    cd ..
    cd('results')
        fileName  = sprintf('%s_PS_%s.mat', method, name_example);
        save(fileName,'m')

        save(sprintf('%s_PStensor_%s.mat', method, name_example),'PS')

        fileName  = sprintf('%s_Error_M12xM21_%s.mat', method, name_example);
        save(fileName,'errorComparison')
    cd ..
    cd('functions')

    %% mesh refinement 
    
    if nsteps < nRefinament 
        if nsteps == 0
            p = mesh.p0;
            t = mesh.t0;
            e = mesh.e0;
        else
            p = p_aux;
            t = t_aux;
            e = e_aux;
        end

        if isequal(meshType,'Irreg')
            % irregular mesh
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

        elseif isequal(meshType,'Local')
            % local refinement
            tagvertex = data.vertex;
            if ~isempty(tagvertex)
                if nsteps == 0
                    r = data.r;
                else
                    r = r_aux/2;
                end
                r_aux = r;

                for ji=1:2                
                    xc = pdeintrp(p,t,p(1,:)');
                    yc = pdeintrp(p,t,p(2,:)');
                    [ntag, ~]=size(tagvertex);
                    it=[];
                    for i=1:ntag
                       it=[it;find(sqrt((tagvertex(i,1)-xc).^2 +(tagvertex(i,2)-yc).^2)< r)'];
                    end
                    it= unique(it);
                    u=zeros(size(p,1));
                    [p,e,t] = refinemesh(g,p,e,t,u, it,'regular'); % refinement mesh
                end
            else
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
        elseif isequal(meshType,'Reg')
                % regular mesh
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

        e_aux = e;
        t_aux = t;
        p_aux = p;
        x_c = mesh.x_c;
        y_c = mesh.y_c;
        p = shiftMesh(p, x_c, y_c);
        if rotationAngle ~= 0 
            np = size(p,2);
            for i = 1:np
                [p(1,i), p(2,i)] = call_rotationCoord(p(1,i), p(2,i), rotationAngle);
            end
        end
        
    end
    disp('----------------');
end

PStensor =m;

end
