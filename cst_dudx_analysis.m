%% Set up
    dataroot = '/data/raeed/project-data/smile/cst-gainlag';
    if ispc
        dataroot = 'C:\Users\Raeed\data\project-data\smile\cst-gainlag';
    end

    file_info = dir(fullfile(dataroot,'library','*COCST*'));
    filenames = horzcat({file_info.name})';

%% Loop through file to extract gains and lags for each lambda into a huge table
filetic = tic;
for filenum = 31%1:length(filenames)
    td = load(fullfile(dataroot,'library',filenames{filenum}));
    td = td.trial_data;

    [~,td_cst] = getTDidx(td,'task','CST','result','R');
    if ~isempty(td_cst)
        td_cst = binTD(td_cst,50);
        td_cst = getDifferential(td_cst,struct('signals','hand_pos','alias','hand_vel'));
        td_cst = getDifferential(td_cst,struct('signals','cursor_pos','alias','cursor_vel'));
    
        td_cst = trimTD(td_cst,{'idx_goCueTime',floor(0.1/td_cst(1).bin_size)},{'idx_rewardTime',floor(0.1/td_cst(1).bin_size)});

        for trialnum = 1:length(td_cst)
            td_cst(trialnum).dudx = td_cst(trialnum).hand_vel(:,1)./td_cst(trialnum).cursor_vel(:,1);
        end
%         fprintf('Finished file %d of %d at time %f\n',filenum,length(filenames),toc(filetic))
    else
        fprintf('No CST trials for file %d\n',filenum)
    end
end

%%
td_bin = binTD(td_cst,'average');
scatter(1:length(td_bin),cat(2,td_bin.dudx),[],cat(2,td_bin.lambda),'filled')

%%
plot(cat(1,td_cst.dudx),'k')