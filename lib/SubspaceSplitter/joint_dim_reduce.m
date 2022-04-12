function [td_out,w_new] = joint_dim_reduce(td_cell,params)
    % find joint space across multiple datasets

    num_dims = 10;
    signals = '';
    combine_before_projection = false; %whether to put data together and renormalize before projection into joint space
    assignParams(who,params);

    % soft normalize and dim reduce separately
    w_cell = cell(size(td_cell));
    for tdnum = 1:length(td_cell)
        td_temp = softNormalize(td_cell{tdnum},struct('signals',signals));
        [~,pca_temp] = dimReduce(td_temp,struct(...
            'algorithm','pca',...
            'signals',signals,...
            'num_dims',num_dims));
    
        % put spaces together and orthogonalize
        w_cell{tdnum} = pca_temp.w(:,1:num_dims);
    end
    [w_new,~,~] = svd(cat(2,w_cell{:}),0);

    % reproject data
    td_out = cell(size(td_cell));
    if combine_before_projection
        td_full = cat(2,td_cell{:});
        td_full = softNormalize(td_full,struct('signals',signals));
        td_full = centerSignals(td_full,struct('signals',signals));
        for trialnum = 1:length(td_full)
            td_full(trialnum).M1_pca_joint = getSig(td_full(trialnum),signals) * w_new;
        end
    
        % split data again
        td_lens = cellfun(@length,td_cell);
        td_idx = [0 cumsum(td_lens)];
        for tdnum = 1:length(td_cell)
            td_out{tdnum} = td_full((td_idx(tdnum)+1):td_idx(tdnum+1));
        end
    else
        for tdnum = 1:length(td_cell)
            td_temp = softNormalize(td_cell{tdnum},struct('signals',signals));
            td_temp = centerSignals(td_temp,struct('signals',signals));
            for trialnum = 1:length(td_temp)
                td_temp(trialnum).M1_pca_joint = getSig(td_temp(trialnum),signals) * w_new;
            end
            td_out{tdnum}=td_temp;
        end
    end
end