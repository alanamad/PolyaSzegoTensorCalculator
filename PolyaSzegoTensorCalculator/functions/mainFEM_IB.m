%**************************************************************************
% MATLAB implemenetation of the computation of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Compute the PS tensor using Finite Element/Interpolated BEM
%
% HISTORY
% A. Amad       10/2018: code implementation
%**************************************************************************

function PStensor = mainFEM_IB(data, mesh, params)

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

name_obj = params.name_obj;

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
errorM12M21 = [];

relError = [ ];

for nsteps = 0:nRefinament
    disp('----------------');
    disp(['Refinement = ' num2str(nsteps) ]) ;
    
    Nelem = size(t,2);

    disp(' ');
    disp(['Number of elements = ', num2str(Nelem)]);
    disp(' ');
    if isequal(name_obj,'ellipsoidal')
        if isequal(method,'IntB')
            if any(t(4, :) == 1)
                match = (e(6, :) == 2);
                e([1,2,6,7], match) = e([2,1,7,6], match);
            end
        elseif isequal(method,'FEM')
            if any(t(4, :) == 1)
                match = (e(6, :) == 1);
                e([1,2,6,7], match) = e([2,1,7,6], match);
            end
        end
    else
        if any(t(4, :) == 2)
            match = (e(6, :) == 1);
            e([1,2,6,7], match) = e([2,1,7,6], match);
        end
    end

    nt = size(t,2);
    if isequal(name_obj,'ellipsoidal')
        connect_el = find(t(4,:)==1);
        connect_bound_obj = (e(6,:)==1);
    else
        connect_el = find(t(4,:)==2);
        connect_bound_obj = (e(6,:)==2);
    end
    el_complement = setdiff(1:nt,connect_el);
    t_obj = t(:,connect_el);
    t_complement = t(:,el_complement);
    e_obj = e(:,connect_bound_obj);
    Npanels = size(e_obj,2); % Number of edges
    Nnodes = size(p,2);
    
    % mesh data
    x = p(1,:);
    y = p(2,:);

    %% Checking mesh

    checkMesh = data.checkMesh;
    if isequal(checkMesh,'Yes')
        checkingMesh(p,t);
    end
    
    %% Compute the adjoint double layer potential and the linear system 

    disp(['Npanels = ' num2str(Npanels) ]) ;
    
    Lamb = lamb * eye(Npanels);
    if isequal(method,'IntB') || isequal(plotObject,'Yes')
        xm=0.5.*(x(e_obj(1,:))+x(e_obj(2,:))); % x-middle point edge
        ym=0.5.*(y(e_obj(1,:))+y(e_obj(2,:))); % y-middle point edge
        lg=sqrt((x(e_obj(2,:))-x(e_obj(1,:))).^2+(y(e_obj(2,:))-y(e_obj(1,:))).^2); % length edge
        nx=(y(e_obj(2,:))-y(e_obj(1,:)))./lg; % normal x
        ny=(-x(e_obj(2,:))+x(e_obj(1,:)))./lg; % normal y

        if isequal(plotObject,'Yes')
            figure(1)
            set(1,'WindowStyle','docked');
            quiver(xm,ym,nx,ny)
            hold on; pdemesh(p,e_obj,t_obj)
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
    end

    if isequal(method,'IntB')
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
        normal = [nx; ny];
        f = A\normal'; % vector solution phi

        phi = zeros(Nnodes,2);
        for i = 1:2
            for j = 1:Nnodes
                dist     = sqrt((xm - x(j)) .^ 2 + (ym - y(j)) .^2);
                Green    = log(dist) / (2*pi);
                phi(j,i) = sum(f(:, i)' .* lg .* Green ./ (mur-1));
            end
        end
        ndof = [ndof; Npanels];
    end
    
    coef_c =1;
    coef_a = 0;
    fsour = 0;
    [Kc,~,~] = assema(p,t_complement,coef_c,coef_a,fsour); % Stiffness matrix (outside object)
    coef_obj =mur;
    coef_a = 0;
    fsour = 0;
    [Kobj,~,~] = assema(p,t_obj,coef_obj,coef_a,fsour); % Stiffness matrix (object)

    K = Kc+Kobj;

    if isequal(method,'FEM')
        phix = zeros(1,Nnodes);
        phiy = zeros(1,Nnodes);

        phi_sour = [phix; phiy];

        [~,nelem] = size(t_obj);
        for i=1:nelem
            nodes(1,i)=t_obj(1,i);
            nodes(2,i)=t_obj(2,i);
            nodes(3,i)=t_obj(3,i);
        end

        jk=0;
        for nel = 1:nelem   % Begin to assemle by element.
            for j=1:3	    % The coordinates of the nodes in the element.
                jk=jk+1;
                el(jk) = nodes(j,nel);   
            end
        end      
        el = unique(el);
        xx=x(el);
        yy=y(el);
        phi_sour(:,el) = [xx;yy];

        fsour = 0;
        [Ksour,~,~] = assema(p,t_obj,1,0,fsour);

        F = -Ksour*phi_sour';

        connect_bound = (e(6,:)~=2);
        e_bound = e(:,connect_bound);

        [~,npres] = size(e_bound);
        eq = e_bound(1,:);
        K(eq,:) = 0;
        K(:,eq) = 0;
        for i=1:npres  % Imposing boundary conditions
            nod = e_bound(1,i);
            K(nod,nod) = 1;
        end
        F(eq,:) = 0;
        phi = K\F;
        ndof = [ndof; size(phi,1)];
    end


    area = sum(pdetrg(p,t_obj));
    
    moment = alpha^2*(mur-1)^2*(phi'*(K)*phi);
    B = alpha^2*(mur-1)*area* eye(2);
    m = B - moment;

    if isequal(method,'FEM')
        PS =[PS; size(phi,1) m(1,1) m(1,2) m(2,1) m(2,2)];
    elseif isequal(method,'IntB')
        PS =[PS; Npanels m(1,1) m(1,2) m(2,1) m(2,2)];
    end
    m11 =[m11; m(1,1)];
    m12 =[m12; m(1,2)];
    m21 =[m21; m(2,1)];
    m22 =[m22; m(2,2)];


    disp(['PS = ' mat2str(m) ]) ;

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
        errorM12 = [errorM12; Error]; %#ok<*AGROW>

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
            [p,e,t] = refinemesh(g,p,e,t,'longest');
            [p,e,t] = refinemesh(g,p,e,t,'longest');
        elseif isequal(meshType,'Local')
            % local refinement
            tagvertex = data.vertex;
            if ~isempty(tagvertex)
                t_min_area = min(pdetrg(p,t));
                t_max_area = max(pdetrg(p,t));
                r = data.r;
                hstar = 5*sqrt(t_max_area/r + t_min_area);
                r = hstar;
                
                [ntag, ~]=size(tagvertex);
                [~,nelem]=size(t);
                it=[];
                for n=1:nelem
                    for i=1:3
                        xy(i,1)=p(1,t(i,n)); %#ok<*AGROW>
                        xy(i,2)=p(2,t(i,n));            
                    end
                    xc=(xy(1,1)+xy(2,1)+xy(3,1))/3;
                    yc=(xy(1,2)+xy(2,2)+xy(3,2))/3;

                    for i=1:ntag
                       if sqrt((tagvertex(i,1)-xc)^2+(tagvertex(i,2)-yc)^2)< r
                           it=[it;n];
                       end
                    end        
                end
                it= unique(it);
                u=zeros(size(p,1));
                [p,e,t] = refinemesh(g,p,e,t,u, it,'longest'); % refinement mesh
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
                [p,e,t] = refinemesh(g,p,e,t);
        end
        
        p_aux = p;
        e_aux = e;
        t_aux = t;
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
