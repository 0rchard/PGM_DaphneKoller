function [P loglikelihood] = LearnCPDsGivenGraph(dataset, G, labels)
%
% Inputs:
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% G: graph parameterization as explained in PA description
% labels: N x 2 true class labels for the examples. labels(i,j)=1 if the 
%         the ith example belongs to class j and 0 elsewhere        
%
% Outputs:
% P: struct array parameters (explained in PA description)
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
K = size(labels,2);

loglikelihood = 0;
P.c = zeros(1,K);

% estimate parameters
% fill in P.c, MLE for class probabilities
% fill in P.clg for each body part and each class
% choose the right parameterization based on G(i,1)
% compute the likelihood - you may want to use ComputeLogLikelihood.m
% you just implemented.
%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_C = sum(labels, 1);
P.c = count_C / sum(count_C);

P.clg = repmat(struct('mu_y', [], 'sigma_y', [], 'mu_x', [], 'sigma_x', [], 'mu_angle', [], 'sigma_angle', [], 'theta', []), 1, 10);

for i=1:10
    if G(i,1) == 0
        P.clg(i).mu_y = zeros(1, K);
        P.clg(i).sigma_y = zeros(1, K);
        P.clg(i).mu_x = zeros(1, K);
        P.clg(i).sigma_x = zeros(1, K);
        P.clg(i).mu_angle = zeros(1, K);
        P.clg(i).sigma_angle = zeros(1, K);
    else
        P.clg(i).sigma_y = zeros(1, K);
        P.clg(i).sigma_x = zeros(1, K);
        P.clg(i).sigma_angle = zeros(1, K);
        P.clg(i).theta = zeros(K, 12);
    end
end


for i=1:10
    for k=1:K
        if G(i,1) == 0
            U = squeeze(dataset(:,i,:));
            U_y = labels(:,k) .* U(:,1);
            U_y = U_y(find(U_y~=0));
            [P.clg(i).mu_y(k), P.clg(i).sigma_y(k)] = FitGaussianParameters(U_y);
            U_x = labels(:,k) .* U(:,2);
            U_x = U_x(find(U_x~=0));
            [P.clg(i).mu_x(k), P.clg(i).sigma_x(k)] = FitGaussianParameters(U_x);
            U_angle = labels(:,k) .* U(:,3);
            U_angle = U_angle(find(U_angle~=0));
            [P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k)] = FitGaussianParameters(U_angle);
        else
            par = G(i,2);
            U_y = squeeze(dataset(:,par,:));            
            X_y = squeeze(dataset(:,i,1));
            U_y = U_y(find(labels(:,k) == 1),:);
            X_y = X_y(find(labels(:,k) == 1));
            [theta_y sigma_y] = FitLinearGaussianParameters(X_y, U_y);
            P.clg(i).theta(k,1:4) = [theta_y(4), theta_y(1), theta_y(2), theta_y(3)];
            P.clg(i).sigma_y(k) = sigma_y;
            %disp(sprintf('POS:%d, K:%d, theta:[%f,%f,%f,%f]', i, k, theta_y(1), theta_y(2), theta_y(3), theta_y(4)));
            
            U_x = squeeze(dataset(:,par,:));
            X_x = squeeze(dataset(:,i,2));
            U_x = U_x(find(labels(:,k) == 1),:);
            X_x = X_x(find(labels(:,k) == 1));
            [theta_x sigma_x] = FitLinearGaussianParameters(X_x, U_x);
            P.clg(i).theta(k,5:8) = [theta_x(4), theta_x(1), theta_x(2), theta_x(3)];
            P.clg(i).sigma_x(k) = sigma_x;
            
            U_angle = squeeze(dataset(:,par,:));
            X_angle = squeeze(dataset(:,i,3));
            U_angle = U_angle(find(labels(:,k) == 1),:);
            X_angle = X_angle(find(labels(:,k) == 1));
            [theta_angle sigma_angle] = FitLinearGaussianParameters(X_angle, U_angle);
            P.clg(i).theta(k,9:12) = [theta_angle(4), theta_angle(1), theta_angle(2), theta_angle(3)];
            P.clg(i).sigma_angle(k) = sigma_angle;
        end
    end
end

loglikelihood = ComputeLogLikelihood(P, G, dataset);
fprintf('log likelihood: %f\n', loglikelihood);

