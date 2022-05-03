
dtstr = '20120528';
E = loadEyeData(dtstr);

%% find avg pupil measurements per trial, and plot smoothed per block

nm = 'Left Pupil';
minTime = 7; % end of freeze period
maxTime = 20; % end of early phase of trial (arbitrary definition)

ceyeind = find(ismember(E.dataNames, 'Left Eye X'));
cind = find(ismember(E.dataNames, nm));
pupil_data = nan(numel(E.Data), 2);
for ii = 1:numel(E.Data)
    eyedata = E.Data(ii).TrialData.analogData;
    pups = eyedata(:,cind);
    xs = 1:size(pups,1); xs = floor(xs/45)+1;
    ix = (xs >= minTime) & (xs <= maxTime);
    
    pupil_data(ii,1) = E.Data(ii).trialTime;
    pupil_data(ii,2) = nanmedian(pups(ix,:),1);
end
mu = nanmean(pupil_data(:,2));
sd = nanstd(pupil_data(:,2));
pupil_data(:,2) = (pupil_data(:,2)-mu)./sd; % z-score

% plot
figure; set(gcf, 'color', 'w'); hold on; set(gca, 'FontSize', 14);
plot(pupil_data(:,1), pupil_data(:,2), '.-');
xlabel('time'); ylabel('pupil (a.u.)'); axis tight;
title(E.datestr);
