%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hangingnode.m
% Author: Jay Oswald, Arizona State University 
% E-mail: j-oswald@asu.edu
% Description:
% This code runs a patch test on a quad-tree refined mesh with hanging
% nodes.  The user choses which elements to refine.  Elements must be
% refined such that there is no more than one hanging node per edge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = hangingnode()
clc
    % User settings.
    show_numbering  = 0;
    more_refinement = 0;    
    
    % Material properties (plane strain).
    E = 1.0; v=0.25;
    C = (E/(1+v)/(1-2*v))*[1-v, v, 0; v, 1-v, 0; 0, 0, 0.5-v];
    % Creates a uniform mesh.
    mesh = mesh2d(3,3,1.0,1.0);
    refine = [5,9,2,3];
    for r = refine
        mesh = refine_element(mesh, r);
    end
    % To see more refinement.    
    if more_refinement
        refine = [5,22,5,6,3,30,35,4,2,30,36,41,33];
        for r = refine
            mesh = refine_element(mesh, r);
        end                
    end
    hanging = hanging_nodes(mesh);
    qpts = [-1,-1,1,1;-1,1,-1,1]/sqrt(3);    
    % Number of non-hanging nodes.
    ncn = mesh.nn - size(hanging,1);
    K = spalloc(2*mesh.nn,2*mesh.nn,2*ncn*18);
    f = zeros(2*mesh.nn,1);
    
    % Assemble stiffness matrix.
    for conn = mesh.connect'
        Ke = zeros(8,8);
        xIe = mesh.xI(conn,:);
        for q = qpts
            dN = gradshape(q);
            J = (dN * xIe)';
            dN = J\dN;
            B = [dN(1,:), 0,0,0,0; 0,0,0,0, dN(2,:); dN(2,:), dN(1,:)];
            Ke = Ke + B'*C*B*det(J);
        end                        
        % Scatter matrix condenses out hanging node degrees of freedom.
        L  = spalloc(2*mesh.nn, 8, 16);      
        for i=1:length(conn)
            I = conn(i);
            row = find(I==hanging(:,1));
            if isempty(row) % (Not on a hanging node).
                L(I,i)      = 1.0;
                L(I+mesh.nn,i+4) = 1.0;                
            else
                shared = hanging(row,2:3);
                L(shared,i)           = 1.0/length(shared);
                L(shared+mesh.nn,i+4) = 1.0/length(shared);
            end
        end
        % Put 1 on diagonal for hanging node K matrix rows.
        for i = hanging(:,1)'
            K(i,i) = 1.0;
            K(i+mesh.nn, i+mesh.nn) = 1.0;
        end        
        K = K + L*Ke*L';        
    end    
    
    % Boundary conditions (on exterior nodes).
    fix = find(mesh.xI(:,1)==0 | mesh.xI(:,1)==1 | ...
               mesh.xI(:,2)==0 | mesh.xI(:,2)==1)';
    for n = [fix, fix+mesh.nn]
        K(n,:) = 0.0;
        K(n,n) = 1.0;
    end
    for n = fix
        f(n:mesh.nn:end) = [0.8, 0.2; 0.2, 1.5]*mesh.xI(n,:)';
    end    
    % Solution of linear system.
    u = K\f;    
    % Post processing
    ux = u(1:end/2);
    uy = u(mesh.nn+1:end);        
    for h=hanging'
        ux(h(1)) = mean(ux(h(2:end)));
        uy(h(1)) = mean(uy(h(2:end)));
    end   
    mesh.xI(:,1) = mesh.xI(:,1) + ux;
    mesh.xI(:,2) = mesh.xI(:,2) + uy;
    % Plot displacement.
    clf(); 
    patch('faces', mesh.connect, 'vertices', mesh.xI, ...
          'facecolor', 'interp', 'facevertexcdata', ux, 'edgealpha', 0.5,'linewidth',1.25);
    axis equal; 
    if show_numbering
        % Plot node numbers.
        for i=1:mesh.nn
            text(mesh.xI(i,1), mesh.xI(i,2), sprintf('%d',i));
        end
        % Plot element numbers
        for e=1:mesh.ne
            xc = mean(mesh.xI(mesh.connect(e,:),:));
            text(xc(1),xc(2), sprintf('%d',e));
        end     
    end
    
end

function [N] = shape(xi)
    x = xi(1); y = xi(2);
    N = 0.25*[(1-x).*(1-y), (1+x).*(1-y), (1+x).*(1+y), (1-x).*(1+y)];
end

function [dN] = gradshape(xi)
    x = xi(1); y = xi(2);
    dN = 0.25*[y-1, 1-y, 1+y, -y-1; x-1, -x-1, 1+x,  1-x];
end

% Returns an Nx3 matrix of node indices.  
% The 1st column is the hanging nodes.
% The 2nd & 3rd columns are the edge nodes.
%      x----E1-----------x   
%      |    |            |  (E) edge node
%      |    |            |  (H) hanging node
%      x----H            |
%      |    |            |
%      |    |            |
%      x----E2-----------x
function [h] = hanging_nodes(m)    
    h = [];    
    for e=1:length(m.connect)           
        for edge_local = [1,2,3,4;2,3,4,1]
            edge = m.connect(e,edge_local);
            % Elements connected to each edge.
            e1 = find(any((m.connect == edge(1))'));
            e1(e1==e) = [];
            e2 = find(any((m.connect == edge(2))'));
            e2(e2==e) = [];
            if ~isempty(intersect(e1,e2))
                continue;
            end
            n1 = reshape(m.connect(e1,:),[],1);
            n2 = reshape(m.connect(e2,:),[],1);
            shared = intersect(n1,n2);
            if length(shared) == 2
                midpoint = mean(m.xI(edge,:));
                d1 = norm(m.xI(shared(1),:) - midpoint);
                d2 = norm(m.xI(shared(2),:) - midpoint);
                if d1<d2, h(end+1,:) = [shared(1), edge];
                else      h(end+1,:) = [shared(2), edge];
                end                                
            end            
        end    
    end    
end

% Splits an element into 4, creating new nodes and connectivity.
function [m] = refine_element(m, e)
    % Nodes of element to refine.
    er = m.connect(e,:);
    % Delete element e.
	m.connect(e,:) = [];
    xIe = m.xI(er,:);
    %  o4  a5  o3    o_ = existing nodes
    %  a2  a3  a4    a_ = (potentially) new nodes
    %  o1  a1  o2    
    a = [mean(xIe([1,2],:))
         mean(xIe([1,4],:));
         mean(xIe([1,3],:));
         mean(xIe([2,3],:));
         mean(xIe([3,4],:))];
    % ai are the ids of new nodes.
    ai = zeros(5,1);
    for i=1:length(a) 
        ai(i) = find_node(m.xI, a(i,:));
        if ai(i) < 0
            m.xI(end+1, :) = a(i,:);
            ai(i)          = length(m.xI);
        end        
    end
    m.connect(end+1,:) = [er(1), ai(1), ai(3), ai(2)];
    m.connect(end+1,:) = [ai(1), er(2), ai(4), ai(3)];
    m.connect(end+1,:) = [ai(3), ai(4), er(3), ai(5)];
    m.connect(end+1,:) = [ai(2), ai(3), ai(5), er(4)];
    m.nn = length(m.xI);
    m.ne = length(m.connect);
end

% Returns the node id if found or -1 if not.
function [id] = find_node(xI, x)
    d = sqrt((xI(:,1)-x(1)).^2+(xI(:,2)-x(2)).^2);
    match = find(d < 1e-8);
    if isempty(match), id = -1;
    else               id = match(1);
    end
end

% Makes a simple 2D uniform mesh.
function [m] = mesh2d(ex, ey, lx, ly)       
    m.ne = ex*ey;
    m.nn = (ex+1)*(ey+1);
    m.xI = zeros(m.nn,2);
    m.connect = zeros(m.ne, 4);
    xx = linspace(0,lx,ex+1);
    yy = linspace(0,ly,ey+1);    
    ct = 0;
    for j=1:ey+1
        for i=1:ex+1
            ct = ct + 1;
            m.xI(ct,:) = [xx(i), yy(j)];            
        end
    end        
    ct = 0;
    for j=1:ey
        for i=1:ex
            ct = ct + 1;
            n0 = i + (j-1)*(ex+1);
            m.connect(ct,:) = [n0, n0+1, n0+ex+2, n0+ex+1];
        end
    end
end