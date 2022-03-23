function [Subspaces,Projs,out] = SubspaceSplitterKaoOrthGreedy(X1,X2,alphaNullSpace,var_cutoff,do_plot)
% SubspaceSplitterKaoOrth performs an iterative rotation of the N-d space
% containing data X1 and X2 to sort dimensions into those with variance
% unique to one condition, or shared by both. 
% 
% Inputs
%   X1             : Data matrix for condition 1 (size M x N: M = samples, N = dimensions)
%   X2             : Data matrix for condition 2 (size(X2,2) must equal size (X1,2))
%   alphaNullSpace : Value dictating percentage of variance allowed to
%                    exist in opposite unique space 
% 
%                     e.g. alphaNullSpace = 0.01 allows 1% of condition 2 
%                     variance to exist in a condition 1 unique axis
%
%                     NOTE: this is really a "best case" constraint. After
%                     identifying all possible unique axes, the function
%                     will return shared axes that violate this constraint
%
%   var_cutoff     : Total % variance cutoff used for defining discrete
%                    uniqeu subspaces. 
%
%                     e.g. var_cutoff = 0.05 means that the unique spaces
%                     will explain a maximum of 5% of the opposite
%                     condition variance. 
%
%   do_plot (optional) : Flag to return summary plot (default = true) 

%%
% Add toolboxes from:
% Jiang, X., Saggar, H., Ryu, S. I., Shenoy, K. V., & Kao, J. C. (2020). 
% Structure in Neural Activity during Observed and Executed Movements Is 
% Shared at the Neural Population Level, Not in Single Neurons. Cell Reports, 32(6), 108006.
% https://doi.org/10.1016/j.celrep.2020.108006
curdir = cd; 
addpath([curdir, '\KaoHelpers']); 

if nargin < 5 || isempty(do_plot)
    do_plot = true; 
end

assert(size(X1,2)==size(X2,2),'X1 and X2 must have the same dimensionality'); 

X1 = X1(~any(isnan(X1),2),:); 
X2 = X2(~any(isnan(X2),2),:); 

v1 = sum(var(X1)); % total overt variance
v2 = sum(var(X2)); % total covert variance

%%
d = size(X1,2); 
tol_cutoff = 2e-5;%1e-5; % 4e-5;
[nullers,pots] = deal(cell(1,d-1)); % Initialize
data1 = X1; % set initial data matrix 1
data2 = X2; % set initial data matrix 2
iterbar = ['[' repmat('-',1,d-1) ']']; iterbar(2) = char(29); % progress bar
cond_tracker = NaN(1,d-1); 
for i = 1:(d-1) % loop through dimensions

    clc; fprintf('Identifying orthogonal spaces:\n%s\n',iterbar); iterbar(i+2) = char(29); 
    

    flag = -1; alphanull = alphaNullSpace; % initialize params
    while flag < 0 % try using defined alphaNullSpace. If it fails, increase and try again
        [Q1,flag] = exclusive_subspace(cov(data1),cov(data2),1,alphanull,tol_cutoff); 
        alphanull = alphanull*1.1; % increase by 10%
    end

    flag = -1; alphanull = alphaNullSpace; 
    while flag < 0
        [Q2,flag] = exclusive_subspace(cov(data2),cov(data1),1,alphanull,tol_cutoff); 
        alphanull = alphanull*1.1; 
    end

    cross_v1 = var(data2*Q1)./v2; % condition 2 variance explained by axis 1
    cross_v2 = var(data1*Q2)./v1; % condition 1 variance explained by axis 2 
    
    if sum(cross_v1) < sum(cross_v2)
        Q = Q1; 
        cond_tracker(i) = 1; 
    else
        Q = Q2; 
        cond_tracker(i) = 2; 
    end
    
    [U,~,~] = svd(Q); % SVD to find full space
    nuller = U(:,2:end); % trailing eigenvectors correspond to null dimensions
    data1 = data1*nuller; % constrain data to new, lower-d subspace
    data2 = data2*nuller; 
        
    nullers{i} = nuller; % Keep track of (reduced space) null spaces
    pots{i} = U(:,1); % Save (reduced space) eigenvector 
end

%% Compile reduced space eigenvectors into full space
Qf = NaN(d,d-1); % Initialize 
Qf(:,1) = pots{1}; 
for i = 2:length(pots) % Loop through dimensions (first dimension is already in full space)
    % calculate matrix for re-projecting reduced space vectors into full dimensional space
    nullerpre = eye(size(nullers{1},1)); 
    for j = 1:(i-1)
        nullerpre = nullerpre*nullers{j}; 
    end
    Qf(:,i) = nullerpre*pots{i}; % Re-project
end    
% do svd to find the remaining last dimension
[U,~,~] = svd(Qf); 
Qf = [Qf, U(:,end)]; 

% Re-order Qf so dimensions are ordered cond 1 --> shared --> cond 2
Qord = [Qf(:,cond_tracker==1) fliplr(Qf(:,cond_tracker==2))]; 

%% Split into Unique and Shared
varexs = [var(X1*Qord)./v1; var(X2*Qord)./v2]*100; % variance explained

% use provided var_cutoff to split
cumvar = [cumsum(varexs(1,end:-1:1)); cumsum(varexs(2,:))];  
inds{1} = cumvar(2,:)<(var_cutoff*100); 
inds{2} = fliplr(cumvar(1,:)<(var_cutoff*100)); 
inds{3} = ~inds{1}&~inds{2}; 

Subspaces = cellfun(@(x) Qord(:,x),inds,'uni',0); 

% order the factors within each subspace from high->low variance
[~,ord1] = sortrows(var(X1*Subspaces{1})','descend'); 
[~,ord2] = sortrows(var(X2*Subspaces{2})','descend'); 
[~,ord3] = sortrows(sum([var(X1*Subspaces{3})./sum(var(X1)); var(X2*Subspaces{3})./sum(var(X2))])','descend'); 

% Subspaces: {Unique1, Unique2, Shared}; 
Subspaces{1} = Subspaces{1}(:,ord1); 
Subspaces{2} = Subspaces{2}(:,ord2); 
Subspaces{3} = Subspaces{3}(:,ord3); 

%% Compile Outputs
out.exp_var.cond1.original = var(X1)./v1*100; 
out.exp_var.cond1.sorted = var(X1*Qord)./v1*100; 
out.exp_var.cond1.split = cellfun(@(x) var(X1*x)./v1*100,Subspaces,'uni',0); 

out.exp_var.cond2.original = var(X2)./v2*100; 
out.exp_var.cond2.sorted = var(X2*Qord)./v2*100; 
out.exp_var.cond2.split = cellfun(@(x) var(X2*x)./v2*100,Subspaces,'uni',0);

out.subspace_inds.unique1 = inds{1}; 
out.subspace_inds.unique2 = inds{2}; 
out.subspace_inds.shared = inds{3}; 

out.axes.sorted = Qord;  

Projs.cond1.sorted = X1*Qord; 
Projs.cond1.subspaces = cellfun(@(x) X1*x,Subspaces,'uni',0); 

Projs.cond1.sorted = X2*Qord; 
Projs.cond1.subspaces = cellfun(@(x) X2*x,Subspaces,'uni',0); 

%% Do plotting if requested
if do_plot
    
    figure('Color','w','Position',[450 219 1039 735]); hold on; 
    
    subplot(3,3,[1 2]); hold on; 
    bar([out.exp_var.cond1.original; out.exp_var.cond2.original]','EdgeColor','none');
    legend({'cond 1','cond 2'},'Location','Best','box','off','FontSize',14); 
    set(gca,'TickDir','out','FontSize',12); title('Original','FontSize',16); 
    ylabel('% variance'); xlim([0 d+1]); 
    
    subplot(3,3,3); hold on; 
    plot(out.exp_var.cond1.original,out.exp_var.cond2.original,'.','MarkerSize',15,'Color',[.5 .5 .5]); 
    axis square; set(gca,'TickDir','out','FontSize',12); xlabel('% var cond 1'); ylabel('% var cond 2'); 
    xl = xlim; yl = ylim; xlim([min([xl yl]),max([xl yl])]); ylim([min([xl yl]),max([xl yl])]); 
    
    subplot(3,3,[4 5]); hold on; 
    bar([out.exp_var.cond1.sorted; out.exp_var.cond2.sorted]','EdgeColor','none'); 
    set(gca,'TickDir','out','FontSize',12); xlabel('Dimension');  title('Rotated','FontSize',16); 
    yl = ylim; ylim(yl); 
    plot((sum(inds{1})+.5)*[1 1],yl,':','Color',[.5 .5 .5],'LineWidth',1.5); 
    plot((sum(inds{1}+inds{3})+.5)*[1 1],yl,':','Color',[.5 .5 .5],'LineWidth',1.5); 
    xlim([0 d+1]); 
    
    subplot(3,3,6); hold on; 
    clrs = [lines(2); .5 .5 .5]; 
    for i = 1:3
        plot(out.exp_var.cond1.split{i},out.exp_var.cond2.split{i},'.','MarkerSize',15,'Color',clrs(i,:)); 
    end
    axis square; set(gca,'TickDir','out','FontSize',12); xlabel('% var cond 1'); ylabel('% var cond 2'); 
    xl = xlim; yl = ylim; xlim([min([xl yl]),max([xl yl])]); ylim([min([xl yl]),max([xl yl])]); 
    
    subplot(3,3,[7 8]); hold on; 
    ev_sub1 = cellfun(@sum,out.exp_var.cond1.split([1 3 2])); 
    ev_sub2 = cellfun(@sum,out.exp_var.cond2.split([1 3 2])); 
    bar([ev_sub1;ev_sub2]','EdgeColor','none'); 
    set(gca,'TickDir','out','FontSize',12,'Xtick',1:3);
    xticklabels({'unique 1','shared','unique 2'}); 
    ylabel('% variance');
end

end