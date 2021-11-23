function handle = make_interactive_td_plot(td,params,plot_params)
% wrapper function to make an interactive plot with trial data
% Inputs:
%   td - trial_data structure to iterate over
%   params - params for overall plotting
%       .shape - 2x1 int array specifying subplot shaping
%       .colormap - colormap to use for plot
%       .suptitle_fcn - function handle to generate supertitle of plotting window
%           <string> = suptitle_fcn(trial)
%       .gcf_params - struct of params to pass to `set(gcf,...)` for figure
%   plot_params - struct array for individual plot paramerters
%       .plot_loc - integer specifying which subplot this plot should be
%       .plot_fcn - function handle to use when plotting. Must be of form
%           plot_fcn(trial)
%       .xlabel - label for x axis of plot
%       .ylabel - label for y axis of plot
%       .title - title of plot
%       .gca_params - struct of params to pass to `set(gca,...)` for individual plot

% TODO: find some smart way to specify how axes should be linked if desired

    handle = figure;
    if isfield(params,'gcf_params')
        set(gcf,params.gcf_params)
    end
    colormap(params.colormap)
    trialnum = 1;
    while trialnum<=length(td)
        clf
        
        for plotnum = 1:length(plot_params)
            subplot(params.shape(1),params.shape(2),plot_params(plotnum).plot_loc);
            plot_params(plotnum).plot_fcn(td(trialnum))
            xlabel(plot_params(plotnum).xlabel)
            ylabel(plot_params(plotnum).ylabel)
            set(gca,plot_params(plotnum).gca_params)
        end
        
        sgtitle(params.suptitle_fcn(td(trialnum)))

        % set up navigation keys
        while true
            if waitforbuttonpress==1
                charpressed = get(gcf,'CurrentCharacter');
                if charpressed == 'h'
                    trialnum = max(trialnum-1,1);
                    break;
                elseif charpressed == 'l'
                    trialnum = min(trialnum+1,length(td));
                    break;
                elseif charpressed == 'q'
                    return
                end
            end
        end
    end