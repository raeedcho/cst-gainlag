function pr = calc_participation_ratio(X)
% calculates the particpation ratio for a given signal matrix X, where rows are observations and columns are features

    X = X-mean(X,1);

    N_T = size(X,1);
    C = (1/N_T)*(X'*X);
    pr = (trace(C)^2)/trace(C^2);
%     eigs = eig(C);
%     pr = (sum(eigs)^2)/(sum(eigs.^2));