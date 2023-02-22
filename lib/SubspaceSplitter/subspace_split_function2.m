function [Q,varexp] = subspace_split_function2(C1,C2,variance_cutoff,do_plot)
options.verbosity = 0; 
options.maxiter = 1000; 
warning('off', 'manopt:getHessian:approx');
if nargin < 4
    do_plot = false; 
end
if nargin < 3 || isempty(variance_cutoff)
    variance_cutoff = 'auto'; 
end
if size(C1,1) > size(C1,2) % data, not covariance
    AO = C1; 
    BO = C2; 
    C1 = cov(AO(~any(isnan(AO),2),:)); 
    C2 = cov(BO(~any(isnan(BO),2),:)); 
end

assert(all(size(C1)==size(C2)),'C1 and C2 must be the same size'); 

problem = []; 
problem.M = rotationsfactory(size(C1,2)); 
problem.cost = @(Q) sharedcostfun(Q,C1,C2); 
problem.egrad = @(Q) sharedgradfun(Q,C1,C2); 

% Qrot = trustregions(problem,[],options); 
[Qinit,~] = pcacov(C1./trace(C1) + C2./trace(C2)); 

C1O = C1; C2O = C2; 
C1 = Qinit'*C1O*Qinit; 
C2 = Qinit'*C2O*Qinit; 

problem = []; 
problem.cost = @(Q) uniquecostfun(Q,C1,C2);
problem.egrad = @(Q) uniquegradfun(Q,C1,C2); 

ds = size(C1,2); 
[Qs,maxleakage_rats] = deal(cell(ds,1)); 
[costs] = deal(NaN(ds,1)); 
[withinv,withinperv,crossv,crossperv,maxleaks] = deal(NaN(ds,2)); 
I = eye(size(C1)); 
for i = 1:ds
    clc; fprintf('shared subspace: %dD\n',ds-i); 
    
    problem.M = stiefelfactory(size(C1,2),i); 
     
    Qs{i} = trustregions(problem,I(:,1:i),options); 

    costs(i) = uniquecostfun(Qs{i},C1,C2); 
    
    vQ1 = diag(Qs{i}'*C1*Qs{i})/trace(C1); 
    vQ2 = diag(Qs{i}'*C2*Qs{i})/trace(C2); 

    dims1 = vQ1>vQ2; 
    dims2 = vQ2>vQ1; 
    
    [~,lat_sub1_c1] = pcacov(Qs{i}(:,dims1)'*C1*Qs{i}(:,dims1)); 
    [~,lat_sub1_c2] = pcacov(Qs{i}(:,dims1)'*C2*Qs{i}(:,dims1)); 
    [~,lat_sub2_c1] = pcacov(Qs{i}(:,dims2)'*C1*Qs{i}(:,dims2)); 
    [~,lat_sub2_c2] = pcacov(Qs{i}(:,dims2)'*C2*Qs{i}(:,dims2)); 
    
    if isempty(lat_sub1_c1); lat_sub1_c1 = NaN; end
    if isempty(lat_sub1_c2); lat_sub1_c2 = NaN; end
    if isempty(lat_sub2_c1); lat_sub2_c1 = NaN; end
    if isempty(lat_sub2_c2); lat_sub2_c2 = NaN; end
    
    v_sub1_c1 = lat_sub1_c1./trace(C1); % percent var explained (C1 in subspace 1) 
    v_sub1_c2 = lat_sub1_c2./trace(C2); % percent var explained (C2 in subspace 1) 
    v_sub2_c1 = lat_sub2_c1./trace(C1); % percent var explained (C1 in subspace 2) 
    v_sub2_c2 = lat_sub2_c2./trace(C2); % percent var explained (C2 in subspace 2) 
    
    maxleakage_ratio_sub1 = v_sub1_c2(1)./v_sub1_c1(end); % (largest leakage) / (smallest within)
    maxleakage_ratio_sub2 = v_sub2_c1(1)./v_sub2_c2(end); 
    
    maxleakage_rats{i} = [maxleakage_ratio_sub1, maxleakage_ratio_sub2]; 
    
    maxleaks(i,:) = [v_sub1_c2(1), v_sub2_c1(1)]; 
    crossv(i,:) = [sum(v_sub1_c2), sum(v_sub2_c1)]; 
    withinv(i,:) = [sum(v_sub1_c1), sum(v_sub2_c2)]; 
    crossperv(i,:) = [mean(v_sub1_c2), mean(v_sub2_c1)]; 
    withinperv(i,:) = [mean(v_sub1_c1), mean(v_sub2_c2)]; 
end
% sharevar = [1,1; 1 - nansum(cat(3,fliplr(withinv(1:end-1,:)),crossv(1:end-1,:)),3)];
% if isnumeric(variance_cutoff)
% %     unique_d = find(all(crossv < variance_cutoff,2),1,'last'); 
%     unique_d = find(all(maxleaks < variance_cutoff,2),1,'last'); 
% else
% %     unique_d = size(C1,2) - elbow_dimensionality(flipud(nansum(crossv,2))); 
% %     unique_d = find(max(withinperv,[],2)==min(max(withinperv,[],2)),1,'first'); 
%     
%     share_d = find(any(flipud(diff(sharevar))>0,2),1,'first'); % # dimensions with strictly increasing shared variance
%     if isempty(share_d); share_d = 0; end
%     unique_d = size(C1,2)-share_d; 
%     
% end
if do_plot
    figure('Position',[146 68 1690 922]); hold on; clrs = lines(10); 
end
d = size(C1,1); 
fs = factor(d); 
spd = [d/fs(end) fs(end)]; 
[V1,V2] = deal(zeros(3,d)); 
[V1c,V2c] = deal(cell(3,d)); 
sph = cell(d,1); 
for i = 1:d
    
    unique_d  = i; 
    
    vQ1 = diag(Qs{unique_d}'*C1*Qs{unique_d})/trace(C1);
    vQ2 = diag(Qs{unique_d}'*C2*Qs{unique_d})/trace(C2); 

    dims_unique1 = vQ1>vQ2; 
    dims_unique2 = vQ2>vQ1; 

    [U,~,~] = svd(Qs{unique_d}); 
    ss_shared = U(:,(unique_d+1):end); 

    [pc_sss] = pcacov(ss_shared'*(C1./trace(C1)+C2./trace(C2))*ss_shared); 
    [pc_ss1] = pcacov(Qs{unique_d}(:,dims_unique1)'*C1*Qs{unique_d}(:,dims_unique1)); 
    [pc_ss2] = pcacov(Qs{unique_d}(:,dims_unique2)'*C2*Qs{unique_d}(:,dims_unique2)); 

    ss_unique1 = Qs{unique_d}(:,dims_unique1)*pc_ss1; 
    ss_unique2 = Qs{unique_d}(:,dims_unique2)*pc_ss2; 
    ss_shared = ss_shared*pc_sss; 

    Q.unique1 = Qinit*ss_unique1; 
    Q.unique2 = Qinit*ss_unique2; 
    Q.shared = Qinit*ss_shared; 

    varexp.unique1_C1 = diag(ss_unique1'*C1*ss_unique1)/trace(C1); 
    varexp.unique1_C2 = diag(ss_unique1'*C2*ss_unique1)/trace(C2); 

    varexp.unique2_C1 = diag(ss_unique2'*C1*ss_unique2)/trace(C1); 
    varexp.unique2_C2 = diag(ss_unique2'*C2*ss_unique2)/trace(C2); 

    varexp.shared_C1 = diag(ss_shared'*C1*ss_shared)/trace(C1); 
    varexp.shared_C2 = diag(ss_shared'*C2*ss_shared)/trace(C2); 

    V1(:,i) = [sum(100*varexp.unique1_C1);sum(100*varexp.unique2_C1);sum(100*varexp.shared_C1)];  
    V2(:,i) = [sum(100*varexp.unique1_C2);sum(100*varexp.unique2_C2);sum(100*varexp.shared_C2)];  
    
    V1c(:,i) = {100*varexp.unique1_C1;100*varexp.unique2_C1;100*varexp.shared_C1};  
    V2c(:,i) = {100*varexp.unique1_C2;100*varexp.unique2_C2;100*varexp.shared_C2};  
    
    if do_plot
        sph{i} = subplot(spd(1),spd(2),i); hold on; 
        plot(100*varexp.unique1_C1,100*varexp.unique1_C2,'.','MarkerSize',20); 
        plot(100*varexp.unique2_C1,100*varexp.unique2_C2,'.','MarkerSize',20); 
        plot(100*varexp.shared_C1,100*varexp.shared_C2,'.','MarkerSize',20,'Color',clrs(4,:)); 

        xl = xlim; yl = ylim; zl = [min([xl yl]), max([xl yl])]; xlim(zl); ylim(zl); 
        axis square; 
        if i == 1
            title(sprintf('%d unique dim',i),'FontWeight','normal'); 
        else
            title(sprintf('%d',i),'FontWeight','normal'); 
        end
    end
end
% if do_plot; make_pretty; end

[Vu,Vsh] = deal(cell(1,d)); 
for i = 1:d
    Vu{i} = [V1c{1,i}, V2c{1,i}; V2c{2,i}, V1c{2,i}]; 
    Vsh{i} = [V1c{3,i}, V2c{3,i}]; 
end

 
atans = cellfun(@(x) atan2(x(:,2),x(:,1)),Vu,'uni',0); 
vspread = cellfun(@(x) max(x) - mean(x),atans);

oppo_vex1 = cellfun(@(x) sum(x),V1c(2,:)); 
oppo_vex2 = cellfun(@(x) sum(x),V2c(1,:)); 
oppo_vex = cellfun(@(x) sum(x(:,2)),Vu); 
% w_vex = cellfun(@(x) sum(x(:,1)),Vu); 
% sh_vex = cellfun(@(x) sum(x(:)),Vsh); 

% shared_d = elbow_dimensionality(fliplr(vspread)); 
shared_d = elbow_dimensionality(fliplr(oppo_vex)); 
unique_d = d-shared_d; 
    
if do_plot
    set(sph{unique_d},'Color',[.8 .8 .8]); 
    sph{unique_d}.Title.FontWeight = 'bold'; 
    
    figure('Position',[567 544 700 350]); hold on; 
    
    plot(oppo_vex1,'LineWidth',.5); 
    plot(oppo_vex2,'LineWidth',.5); 
    plot(oppo_vex/2,'k.-'); 
    
    plot([unique_d unique_d],[0 oppo_vex(unique_d)/2],'Color','g','LineWidth',2); 
    plot(unique_d,oppo_vex(unique_d)/2,'.','Color','g','MarkerSize',30); 
    xlabel('# Unique Dimensions'); 
    ylabel('Average opposite variance'); 
%     make_pretty; 
end 


vQ1 = diag(Qs{unique_d}'*C1*Qs{unique_d})/trace(C1);
vQ2 = diag(Qs{unique_d}'*C2*Qs{unique_d})/trace(C2); 

dims_unique1 = vQ1>vQ2; 
dims_unique2 = vQ2>vQ1; 

[U,~,~] = svd(Qs{unique_d}); 
ss_shared = U(:,(unique_d+1):end); 

[pc_sss] = pcacov(ss_shared'*(C1./trace(C1)+C2./trace(C2))*ss_shared); 
[pc_ss1] = pcacov(Qs{unique_d}(:,dims_unique1)'*C1*Qs{unique_d}(:,dims_unique1)); 
[pc_ss2] = pcacov(Qs{unique_d}(:,dims_unique2)'*C2*Qs{unique_d}(:,dims_unique2)); 

ss_unique1 = Qs{unique_d}(:,dims_unique1)*pc_ss1; 
ss_unique2 = Qs{unique_d}(:,dims_unique2)*pc_ss2; 
ss_shared = ss_shared*pc_sss; 

Q.unique1 = Qinit*ss_unique1; 
Q.unique2 = Qinit*ss_unique2; 
Q.shared = Qinit*ss_shared; 

varexp.unique1_C1 = diag(ss_unique1'*C1*ss_unique1)/trace(C1); 
varexp.unique1_C2 = diag(ss_unique1'*C2*ss_unique1)/trace(C2); 

varexp.unique2_C1 = diag(ss_unique2'*C1*ss_unique2)/trace(C1); 
varexp.unique2_C2 = diag(ss_unique2'*C2*ss_unique2)/trace(C2); 

varexp.shared_C1 = diag(ss_shared'*C1*ss_shared)/trace(C1); 
varexp.shared_C2 = diag(ss_shared'*C2*ss_shared)/trace(C2); 

%  plot(100*varexp.unique1_C1,100*varexp.unique1_C2,'.','MarkerSize',20); 
%     plot(100*varexp.unique2_C1,100*varexp.unique2_C2,'.','MarkerSize',20); 
%     plot(100*varexp.shared_C1,100*varexp.shared_C2,'.','MarkerSize',20,'Color',clrs(4,:)); 
% end; make_pretty
%%
if do_plot
    clrs = lines(4); 
    figure('Position',[370 507 1138 433]); hold on; subplot(1,2,1); hold on; 
    plot(100*varexp.unique1_C1,100*varexp.unique1_C2,'.','MarkerSize',20); 
    plot(100*varexp.unique2_C1,100*varexp.unique2_C2,'.','MarkerSize',20); 
    plot(100*varexp.shared_C1,100*varexp.shared_C2,'.','MarkerSize',20,'Color',clrs(4,:)); 
    xlabel('C1 %var explained'); 
    ylabel('C2 %var explained'); 
    axis square; 
    xl = xlim; yl = ylim; xyl = [0, max([xl yl])]; xlim(xyl); ylim(xyl); 
    plot([0 xyl(2)],[0 xyl(2)],'k:'); 
    
    lax = legend({'unique1','unique2','shared'},'location','best'); 
    
    s2 = subplot(2,2,2); hold on; 
    v1 = [varexp.unique1_C1,varexp.unique1_C2]; 
    v2 = [varexp.unique2_C1,varexp.unique2_C2]; 
    vs = [varexp.shared_C1,varexp.shared_C2]; 
    bh = bar([v1;[NaN,NaN];vs;[NaN,NaN];v2],'EdgeColor','none','BarWidth',1); 
    set(gca,'XTick',[]);%1:size(Qs{unique_d},1)); 
    bh(1).ShowBaseLine = 'off'; 
    yl = ylim; 
    plot([.5 size(v1,1)+.5],-0.05*range(yl)*[1 1],'Color',clrs(1,:),'LineWidth',2); 
    plot([.5 size(vs,1)+.5]+size(v1,1)+1,-0.05*range(yl)*[1 1],'Color',clrs(4,:),'LineWidth',2); 
    plot([.5 size(v2,1)+.5]+size([v1;vs],1)+2,-0.05*range(yl)*[1 1],'Color',clrs(2,:),'LineWidth',2); 
    text(mean([.5 size(v1,1)+.5]),-0.05*range(yl),'unique 1','VerticalAlignment','top','HorizontalAlignment','center','FontName','Segoe UI','FontSize',14,'Color',clrs(1,:)); 
    text(mean([.5 size(vs,1)+.5]+size(v1,1)+1),-0.05*range(yl),'shared','VerticalAlignment','top','HorizontalAlignment','center','FontName','Segoe UI','FontSize',14,'Color',clrs(4,:)); 
    text(mean([.5 size(v2,1)+.5]+size([v1;vs],1)+2),-0.05*range(yl),'unique 2','VerticalAlignment','top','HorizontalAlignment','center','FontName','Segoe UI','FontSize',14,'Color',clrs(2,:)); 
    xlim([0 size(Qs{1},1)+3]); 
    ylabel('VAF'); 
    
    subplot(2,2,4); hold on; 
    vex_c1 = 100*[sum(varexp.unique1_C1), sum(varexp.shared_C1), sum(varexp.unique2_C1)]';
    vex_c2 = 100*[sum(varexp.unique1_C2), sum(varexp.shared_C2), sum(varexp.unique2_C2)]'; 
    bar([1 2 3],[vex_c1, vex_c2],'EdgeColor','none','FaceAlpha',0.5); 
    set(gca,'XTick',[1 2 3]); 
    xticklabels({'Unique1','Shared','Unique2'}); 
    ylabel('VAF'); 
%     make_pretty; 
    
    set(lax,'box','on'); 
    s2.XAxis.Visible = 'off'; 
end
clc; fprintf('Done\n'); 

end

%% Helper functions
function qcost = uniquecostfun(Q,C1,C2)
    vQ1 = diag(Q'*C1*Q)/trace(C1); 
    vQ2 = diag(Q'*C2*Q)/trace(C2); 

%     qcost = sum(1-cos(4*atan2(vQ2,vQ1))); 
    
    qcost = sum((0.5-0.5*cos(4*atan2(vQ2,vQ1)))-(vQ1-vQ2).^2); 
    
end

function grad = uniquegradfun(Q,A,B)
%     t0 = 1/trace(B); 
%     t1 = diag(Q'*B*Q); 
%     t2 = trace(A); 
%     t3 = diag(Q'*A*Q); 
%     t4 = (1/t2)*t3; 
%     t5 = t0*t1./t4; 
%     t6 = sin(4*atan2(t0*t1,t4));
%     t7 = 1./(t5.^2 + 1); 
%     t8 = 4*t6.*t7./t4; 
%     t9 = t0*t1.*t6.*t7./(1/(t2.^2)*t3.*t3); 
% 
%     grad = t0*B*Q*diag(t8) + 1/trace(B')*B'*Q*diag(t8)-(4/t2*A*Q*diag(t9)+4/trace(A')*A'*Q*diag(t9)); 

    t0 = trace(B); 
    t1 = 1./t0; 
    t2 = diag(Q'*B*Q); 
    t3 = t1*t2; 
    t4 = trace(A'); 
    t5 = diag(Q'*A*Q); 
    t6 = 1./t4*t5; 
    t7 = t3./t6;
    t8 = trace(A); 
    t9 = (1./t8)*t5; 
    t10 = t3./t9; 
    t11 = -1; 
    t12 = 2./t8; 
    t13 = A*Q; 
    t14 = sin(4*atan2(t3,t6)); 
    t15 = (t7.^2 + 1).^t11;
    t16 = sin(4*atan2(t3,t9)); 
    t17 = (t10.^2 + 1).^t11; 
    t18 = t5.*t5; 
    t19 = B*Q; 
    t20 = t6-t3; 
    t21 = trace(B'); 
    t22 = B'*Q; 
    t23 = diag(t9-t3); 

    grad = -(t12*t13*diag(t20)-(t1*t19*diag(2*t14.*t15./t6)+(1/t21)*t22*...
           diag(2*t16.*t17./t9)-(t12*t13*diag(t1*t2.*t14.*t15./((1/(t4.^2))*t18))+(4/t8)*...
           t13*diag(t1*t2.*.5.*t16.*t17./(1/(t8.^2)*t18))))+t12*t13*t23-(2/t0*t19*diag(t20)+2/t21*t22*t23)); 
end

function qcost = sharedcostfun(Q,C1,C2)
    vQ1 = diag(Q'*C1*Q)/trace(C1); 
    vQ2 = diag(Q'*C2*Q)/trace(C2); 

    qcost = sum(cos(4*atan2(vQ2,vQ1)));     
end

function grad = sharedgradfun(Q,A,B)
    t0 = 1/trace(B); 
    t1 = diag(Q'*B*Q); 
    t2 = trace(A); 
    t3 = diag(Q'*A*Q); 
    t4 = (1/t2)*t3; 
    t5 = t0*t1./t4; 
    t6 = sin(4*atan2(t0*t1,t4));
    t7 = 1./(t5.^2 + 1); 
    t8 = 4*t6.*t7./t4; 
    t9 = t0*t1.*t6.*t7./(1/(t2.^2)*t3.*t3); 

    grad = -(t0*B*Q*diag(t8) + 1/trace(B')*B'*Q*diag(t8)-(4/t2*A*Q*diag(t9)+4/trace(A')*A'*Q*diag(t9))); 

end

function ndims = elbow_dimensionality(scree)

    N = length(scree); 
    d2_line = @(x,y) norm(cross([1 scree(1) 0] - [length(scree) scree(end) 0],[x y 0] -...
        [length(scree) scree(end) 0])) / norm([1 scree(1) 0] - [length(scree) scree(end) 0]); 

    ds = NaN(N,1); 
    for k = 1:N
        ds(k) = d2_line(k,scree(k)); 
    end

    ndims = find(ds==max(ds),1,'first'); 
end