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

function dxdt = cst_update(lambda,x,u)
    dxdt = lambda*(x + u);
end