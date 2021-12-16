function plot_cst_gainlag_phase(params)

    lambda_sweep = 2:2:10;
    gain = -1.2;
    lag = 0.1;
    if nargin>0
        assignParams(who,params)
    end
    
    x_init = 0;
    u_init = 10*rand-5;
    
    figure('defaultaxesfontsize',18)
    for lambdanum = 1:length(lambda_sweep)
        lambda = lambda_sweep(lambdanum);
        [~,x,u] = integrate_gainlag_cst(lambda,x_init,u_init,gain,lag);
        
        subplot(1,length(lambda_sweep),lambdanum)
        plot_cst_phase(x,u)
        title(sprintf('\\lambda = %f',lambda))
    end
    suptitle(sprintf('Gain: %f, lag: %f s',gain,lag))
end

function plot_cst_phase(x,u)
    plot([-60 60],[0 0],'-k','linewidth',1)
    hold on
    plot([0 0],[-60 60],'-k','linewidth',1)
    plot([-60 60],[60 -60],'-k','linewidth',1)
    plot([-50 -50],[-60 60],'-r','linewidth',1)
    patch(...
        [0 0 -60 60 0],...
        [-60 60 60 -60 -60],...
        [0.8 0.8 0.8],'edgecolor','none')
    if max(abs(x))>50
        plot([-50 -50],[-60 60],'-r','linewidth',1)
        plot([50 50],[-60 60],'-r','linewidth',1)
    else
        plot([-50 -50],[-60 60],'-g','linewidth',1)
        plot([50 50],[-60 60],'-g','linewidth',1)
    end
    scatter(x,u,[],viridis(length(u)),'filled')
    axis equal
    set(gca,'box','off','tickdir','out','xlim',[-60 60],'ylim',[-60 60])
    xlabel('Cursor position (mm)')
    ylabel('Hand position (mm)')
end