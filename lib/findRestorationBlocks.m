function td_cst = findRestorationBlocks(td_cst)
% finds start and end indices of blocks where cursor is heading back to center

    for trialnum = 1:length(td_cst)
        % for restoration, cursor velocity and cursor position have opposite sign
        is_restoring = (sign((td_cst(trialnum).hand_pos(:,1) + td_cst(trialnum).cursor_pos(:,1)).*td_cst(trialnum).cursor_pos(:,1)) < 0);
        idx_restoreStart = find(diff(is_restoring)>0)+1;
        idx_restoreEnd = find(diff(is_restoring)<0)+1;

        % append or prepend if we start or end in restoration zone
        if is_restoring(1)>0
            idx_restoreStart = [1; idx_restoreStart];
        end
        if is_restoring(end)>0
            idx_restoreEnd = [idx_restoreEnd; length(is_restoring)];
        end

        % make sure there's the same number of starts as ends and ends come after starts
        assert(length(idx_restoreStart)==length(idx_restoreEnd) && all(idx_restoreStart<idx_restoreEnd),'Something went wrong with the restore index finding');

        td_cst(trialnum).idx_restoreStart = idx_restoreStart';
        td_cst(trialnum).idx_restoreEnd = idx_restoreEnd';
        td_cst(trialnum).is_restoring = is_restoring;
    end