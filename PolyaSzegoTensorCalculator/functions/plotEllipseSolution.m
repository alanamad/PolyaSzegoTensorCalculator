%**************************************************************************
% HISTORY
% A. Amad          06/2018: code implementation
%**************************************************************************

function plotEllipseSolution(method, mur, m,M_exact, alpha, params, Npanels)
name_example = params.name_example;

figure(7); clf;
set(7,'WindowStyle','docked');
        [Vect_exac, Lamb_exac] = eig(M_exact);
        Lamb_exac = sort(diag(Lamb_exac));

        vol_exac = (Lamb_exac(1)*Lamb_exac(2)*(mur+1))/((mur-1)*(Lamb_exac(1)+Lamb_exac(2)));
        q_exac = (Lamb_exac(2)-mur*Lamb_exac(1))/(Lamb_exac(1)-mur*Lamb_exac(2));
        a =  sqrt(abs(vol_exac)/(pi*q_exac));
        b =  sqrt((abs(vol_exac)*q_exac)/pi);
        u = Vect_exac(:,1);
        v = Vect_exac(:,2);

        Theta = atan(v(1)/u(1));
        
        ang = Theta;
        ellipse(a/alpha,b/alpha,ang/alpha,0,0,'k',500)

        
        hold on
       [Vect_h, Lamb_h] = eig(m);
        Lamb_h = sort(diag(Lamb_h));

        vol_h = (Lamb_h(1)*Lamb_h(2)*(mur+1))/((mur-1)*(Lamb_h(1)+Lamb_h(2)));
        q_h = (Lamb_h(2)-mur*Lamb_h(1))/(Lamb_h(1)-mur*Lamb_h(2));
        a =  sqrt(abs(vol_h)/(pi*q_h));
        b =  sqrt((abs(vol_h)*q_h)/pi);
        u = Vect_h(:,1);
        v = Vect_h(:,2);


        
        Theta = atan2(u(1),v(1));

        ang = Theta;
        ellipse(a/alpha,b/alpha,ang/alpha,0,0,'r',500)

    leg = legend("$\mathcal{T}^{exact}$", "$\mathcal{T}^{h}$");
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    set(gca,'XColor','k','YColor','k','LineWidth',2); 
    set(gca,'XTick',[],'YTick',[],'Ztick',[]);
    box on;
    

    cd ..
    cd('results')
         savefig(7,sprintf('%s_tensorAsEllipse_%s_Npanels%d.fig', method, name_example,Npanels));

        filename = sprintf('%s_tensorAsEllipse_%s_Npanels%d.pdf', method, name_example,Npanels);
        set(7,'Units','Inches');
        pos = get(7,'Position');
        set(7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(7, filename,'-dpdf','-r0')
    cd ..
    cd('functions')
end
    
