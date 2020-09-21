%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';

%% Loop through files to check out temporal smoothness at different lambdas
filetic = tic;
for filenum = 1:length(filenames)
    td = load(fullfile(dataroot,'library',filenames{filenum}));
    td = td.trial_data;
    
    % get rid of unsorted neurons
    bad_units = td(1).M1_unit_guide(:,2)<=1;
    % for file 31 only! remove neurons 8,2 and 64,2
    if contains(td(1).date_time,'2018/06/26')
        corr_units = [8 2;64 2];
        bad_units = bad_units | ismember(td(1).M1_unit_guide,corr_units,'rows');
    end
    for trialnum = 1:length(td)
        td(trialnum).M1_spikes = td(trialnum).M1_spikes(:,~bad_units);
        td(trialnum).M1_unit_guide = td(trialnum).M1_unit_guide(~bad_units,:);
    end
    
    if isempty(td(1).M1_unit_guide)
        fprintf('Skipping file %d because no spike data...\n',filenum)
        continue
    end
    
    [~,td_cst] = getTDidx(td,'task','CST','result','R');
    [~,td_co] = getTDidx(td,'task','CO','result','R');
    if ~isempty(td_cst) && ~isempty(td_co)
        td_cst = trimTD(td_cst,{'idx_goCueTime',0},{'idx_goCueTime',floor(6/td_cst(1).bin_size)});

        td_cst = removeBadNeurons(td_cst,struct(...
            'arrays',{{'M1'}},...
            'min_fr',0.1,...
            'calc_fr',true,...
            'do_shunt_check',true));
        
        
        % check different lambdas for different correlations
        lambdas = unique([td_cst.lambda]);
        covars_cst = cell(1,length(lambdas));
        for lambdanum = 1:length(lambdas)
            covars_cst{lambdanum} = pairwiseCorr(td_cst,struct(...
                'signals','M1_spikes',...
                'method',@cov,...
                'zero_diagonal',false,...
                'trial_idx',getTDidx(td_cst,'lambda',lambdas(lambdanum))));
        end
        
        cst_prs = cellfun(@(x) trace(x)^2/trace(x^2),covars_cst);
    
        % now CO
        td_co = trimTD(td_co,{'idx_goCueTime',0},{'idx_goCueTime',floor(0.5/td_co(1).bin_size)});
        
        td_co = removeBadNeurons(td_co,struct(...
            'arrays',{{'M1'}},...
            'min_fr',0.1,...
            'calc_fr',true,...
            'do_shunt_check',true));
        
        % check different lambdas for different correlations
%         ot_locs = unique(vertcat(td_co.ot_location),'rows');
%         covars_co = cell(1,size(ot_locs,1));
%         for otnum = 1:size(ot_locs,1)
%             covars_co{otnum} = pairwiseCorr(td_co,struct(...
%                 'signals','M1_spikes',...
%                 'method',@cov,...
%                 'zero_diagonal',false,...
%                 'trial_idx',getTDidx(td_co,'ot_location',ot_locs(otnum,:))));
%         end
        covar_co = pairwiseCorr(td_cst,struct(...
            'signals','M1_spikes',...
            'method',@cov,...
            'zero_diagonal',false));
        co_pr = trace(covar_co)^2/trace(covar_co^2);
        
        figure
        plot([1.5 4.5],[co_pr co_pr],'--k','linewidth',2)
        hold on
        plot(lambdas,cst_prs,'-k','linewidth',2)
        text(1.5,co_pr,'Center-out PR')
        xlabel('\lambda')
        ylabel('Participation ratio')
        title(sprintf('%s %s',td(1).monkey,td(1).date_time))
        set(gca,'box','off','tickdir','out')
        
        fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
    else
        fprintf('Incomplete dataset for file %d\n',filenum)
    end
end