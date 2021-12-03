function pa_dim = parallel_analysis_dimensionality(X,params)
% PARALLEL_ANALYSIS uses shuffling to determine a null distribution for the eigenvalue spectrum,
% allowing us to set a noise floor on the spectrum to help in deciding a dimensionality for data
% Note: this assumes that X has observations on different rows and features (e.g. neurons) on cols

    num_shuffles = 1000;
    do_plot = false;
    if nargin>1
        assignParams(who,params)
    end
    
    X = X-mean(X,1);
    
    % First, get eigenvalue spectrum of data covariance
    true_eig_spectrum = get_covar_eigs(X);
    
    shuffle_eig_spectra = zeros(num_shuffles,size(X,2));
    for shufflenum = 1:num_shuffles
        X_shuffle = X;
        % shuffle each col independently
        for colnum = 1:size(X,2)
            shuffle_idx = randperm(size(X,1));
            X_shuffle(:,colnum) = X(shuffle_idx,colnum);
        end
        
        shuffle_eig_spectra(shufflenum,:) = get_covar_eigs(X_shuffle);
    end
    
    boundary_eig_spectrum = prctile(shuffle_eig_spectra,95);
    
    if do_plot
        figure('defaultaxesfontsize',10)
        patch([1:size(X,2) fliplr(1:size(X,2))],[boundary_eig_spectrum zeros(1,size(X,2))],[0.8 0.8 0.8],'edgecolor','none')
        hold on
        plot(true_eig_spectrum,'-k')
        set(gca,'box','off','tickdir','out','xlim',[0 length(true_eig_spectrum)])
        ylabel('Eigenvalue')
        xlabel('Eig num')
    end
    
    pa_dim = sum(true_eig_spectrum>boundary_eig_spectrum);
end

function eig_spectrum = get_covar_eigs(X)
% get eigenvalue spectrum of signal covariance

    N_T = size(X,1);
    C = (1/N_T)*(X'*X);
    eig_spectrum = sort(eig(C)','descend');
end