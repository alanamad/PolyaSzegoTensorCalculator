%**************************************************************************
% MATLAB implemenetation of the computation of the Polya-Szego tensor
%**************************************************************************
%  
% DESCRIPTION
% Plots the coefficients of the computed PS tensor
%
% HISTORY
% A. Amad       10/2018: code implementation
%**************************************************************************

function call_plot(method, ndof, params, data, plotData)
    name_obj = params.name_obj;
    meshType = data.meshOption;
    name = params.name_example;
    name_example = sprintf('%s_reftype_%s', name, meshType);

    error = plotData.error;
    
    m11 = plotData.m11;
    m12 = plotData.m12;
    m21 = plotData.m21;
    m22 = plotData.m22;
if isequal(name_obj,'ellipsoidal')

    relError = plotData.relError;
    errorM12 = plotData.errorM12;
    errorM21 = plotData.errorM21;
    errorM12M21 = plotData.errorM12M21;

    Mex11 = plotData.Mex11;
    Mex12 = plotData.Mex12;
    Mex21 = plotData.Mex21;
    Mex22 = plotData.Mex22;
    M_exact = plotData.M_exact; %#ok<*NASGU>
    
    figure(2)
    set(2,'WindowStyle','docked'); 
    loglog(ndof,error,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend("$\|\mathcal{T}^{exact} - \mathcal{T}^{h} \|$");
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('$$ \|\mbox{Error}\|$$','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(21)
    set(21,'WindowStyle','docked');
    loglog(ndof,relError,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend("$\|\mathcal{T}^{exact} - \mathcal{T}^{h} \|/\|\mathcal{T}^{exact} \|$");
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('$$ \mbox{Relative Error}$$','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(22)
    set(22,'WindowStyle','docked');
    loglog(ndof,errorM12,'r*-','LineWidth',1.5,'MarkerSize',5)
    hold on
    loglog(ndof,errorM21,'b--','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\|\mathcal{T}^{exact}_{12} - \mathcal{T}^{h}_{12} \|$', '$\|\mathcal{T}^{exact}_{21} - \mathcal{T}^{h}_{21} \|$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('$$ \|\mbox{Error}\|$$','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(23)
    set(23,'WindowStyle','docked'); 
    loglog(ndof,error,'r*-','LineWidth',1.5,'MarkerSize',5)
    hold on
    loglog(ndof,errorM12M21,'b*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend("$\|\mathcal{T}^{exact} - \mathcal{T}^{h} \|$",'$\|\mathcal{T}^{h}_{12} - \mathcal{T}^{h}_{21}\|$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('$$ \|\mbox{Error}\|$$','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(3)
    set(3,'WindowStyle','docked');
    semilogx(ndof,Mex11,'k-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m11,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{exact}_{11}$', '$\mathcal{T}^{h}_{11}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)


    figure(4)
    set(4,'WindowStyle','docked');
    semilogx(ndof,Mex22,'k-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m22,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{exact}_{22}$', '$\mathcal{T}^{h}_{22}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(5)
    set(5,'WindowStyle','docked');
    semilogx(ndof,Mex12,'k-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m12,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{exact}_{12}$', '$\mathcal{T}^{h}_{12}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(6)
    set(6,'WindowStyle','docked');
    semilogx(ndof,Mex21,'k-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m21,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{exact}_{21}$', '$\mathcal{T}^{h}_{21}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    cd ..
    cd('results')
        fileName  = sprintf('PS_Exact_%s.mat', name_example);
        save(fileName,'M_exact')

        savefig(2,sprintf('%s_errorPS_%s.fig', method, name_example));
        
        filename = sprintf('%s_errorPS_%s.pdf', method, name_example);
        set(2,'Units','Inches');
        pos = get(2,'Position');
        set(2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(2, filename,'-dpdf','-r0')

        savefig(21,sprintf('%s_relErrorPS_%s.fig', method, name_example));
        
        filename = sprintf('%s_relErrorPS_%s.pdf', method, name_example);
        set(21,'Units','Inches');
        pos = get(21,'Position');
        set(21,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(21, filename,'-dpdf','-r0')

        savefig(22,sprintf('%s_error_M12xM21_PS_%s.fig', method, name_example));
        
        filename = sprintf('%s_error_M12xM21_PS_%s.pdf', method, name_example);
        set(22,'Units','Inches');
        pos = get(22,'Position');
        set(22,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(22, filename,'-dpdf','-r0')

        savefig(23,sprintf('%s_errorPS_and_M12xM21_%s.fig', method, name_example));
        
        filename = sprintf('%s_errorPS_and_M12xM21_%s.pdf', method, name_example);
        set(23,'Units','Inches');
        pos = get(23,'Position');
        set(23,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(23, filename,'-dpdf','-r0')

        savefig(3,sprintf('%s_PS_M11_%s.fig', method, name_example));

        filename = sprintf('%s_PS_M11_%s.pdf', method, name_example);
        set(3,'Units','Inches');
        pos = get(3,'Position');
        set(3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(3, filename,'-dpdf','-r0')

        savefig(4,sprintf('%s_PS_M22_%s.fig', method, name_example));

        filename = sprintf('%s_PS_M22_%s.pdf', method, name_example);
        set(4,'Units','Inches');
        pos = get(4,'Position');
        set(4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(4, filename,'-dpdf','-r0')

        savefig(5,sprintf('%s_PS_M12_%s.fig', method, name_example));

        filename = sprintf('%s_PS_M12_%s.pdf', method, name_example);
        set(5,'Units','Inches');
        pos = get(5,'Position');
        set(5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(5, filename,'-dpdf','-r0')

        savefig(6,sprintf('%s_PS_M21_%s.fig', method, name_example));

        filename = sprintf('%s_PS_M21_%s.pdf', method, name_example);
        set(6,'Units','Inches');
        pos = get(6,'Position');
        set(6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(6, filename,'-dpdf','-r0')

    cd ..
    cd('functions')
else
    
    figure(2)
    set(2,'WindowStyle','docked');
    semilogx(ndof,m12,'r*-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m21,'b*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{h}_{12}$', '$\mathcal{T}^{h}_{21}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(3)
    set(3,'WindowStyle','docked');
    loglog(ndof,error,'r*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$|\mathcal{T}^{h}_{12} - \mathcal{T}^{h}_{21}|$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('$$ |\mbox{Error}|$$','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    figure(4)
    set(4,'WindowStyle','docked');
    semilogx(ndof,m11,'r*-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m12,'g*-','LineWidth',1.5,'MarkerSize',5)
    semilogx(ndof,m21,'b*-','LineWidth',1.5,'MarkerSize',5)
    semilogx(ndof,m22,'m*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{h}_{11}$','$\mathcal{T}^{h}_{12}$','$\mathcal{T}^{h}_{21}$','$\mathcal{T}^{h}_{22}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)
    
    figure(5)
    set(5,'WindowStyle','docked');
    semilogx(ndof,m11,'r*-','LineWidth',1.5,'MarkerSize',5)
    hold on
    semilogx(ndof,m22,'b*-','LineWidth',1.5,'MarkerSize',5)
    leg = legend('$\mathcal{T}^{h}_{11}$','$\mathcal{T}^{h}_{22}$');
    set(leg,'Interpreter','latex','Location','best');
    set(leg,'FontSize',13);
    xlabel('$$ \sharp (ndof) $$','Interpreter','latex','FontSize',13)
    ylabel('Tensor coefficient','Interpreter','latex','FontSize',13)
    set(gca,'TickLabelInterpreter','latex')
    set(gca, 'FontSize', 13)

    cd ..
    cd('results')
        savefig(2,sprintf('%s_PS_coeff_M12_and_M21_%s.fig', method, name_example));

        filename = sprintf('%s_PS_coeff_M12_and_M21_%s.pdf', method, name_example);
        set(2,'Units','Inches');
        pos = get(2,'Position');
        set(2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(2, filename,'-dpdf','-r0')

        savefig(3,sprintf('%s_errorPS_M12xM21_%s.fig', method, name_example));

        filename = sprintf('%s_errorPS_M12xM21_%s.pdf', method, name_example);
        set(3,'Units','Inches');
        pos = get(3,'Position');
        set(3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(3, filename,'-dpdf','-r0')

        savefig(4,sprintf('%s_PScoeff_%s.fig', method, name_example));

        filename = sprintf('%s_PScoeff_%s.pdf', method, name_example);
        set(4,'Units','Inches');
        pos = get(4,'Position');
        set(4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(4, filename,'-dpdf','-r0')

        savefig(5,sprintf('%s_PS_coeff_M11_and_M22_%s.fig', method, name_example));

        filename = sprintf('%s_PS_coeff_M11_and_M22_%s.pdf', method, name_example);
        set(5,'Units','Inches');
        pos = get(5,'Position');
        set(5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(5, filename,'-dpdf','-r0')

    cd ..
    cd('functions')
end



end