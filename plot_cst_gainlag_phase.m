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
    
function dxdt = cst_update(lambda,x,u)
    dxdt = lambda*(x + u);
end

function [t,x,u] = integrate_gainlag_cst(lambda,x_init,u_init,gain,lag)
    % forward integrate with Euler
    t_init = 0;
    t_step = 0.01;
    t_final = 6;
    lag_bins = floor(lag/t_step);
    t = (t_init:t_step:t_final)';
    [x,u] = deal(nan(size(t)));
    x(1) = x_init;
    u(1) = u_init;
    for idx = 2:length(t)
        prev_idx = idx-1;
        if t(idx-1)>lag
%             u(idx) = gain*x(idx-lag_bins)+0.05*randn;
            u(idx) = gain*x(idx-lag_bins);
        else
%             u(idx) = u(prev_idx)+0.05*randn;
            u(idx) = u(prev_idx);
        end
        
        dx = t_step*cst_update(lambda,x(prev_idx),u(prev_idx));
        
        x(idx) = x(prev_idx)+dx;
    end
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