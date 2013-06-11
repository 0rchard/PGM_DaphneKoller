function [P G loglikelihood] = LearnGraphAndCPDs(dataset, labels)

% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% labels: N x 2 true class labels for the examples. labels(i,j)=1 if the
%         the ith example belongs to class j
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
K = size(labels,2);

G = zeros(10,2,K); % graph structures to learn
% initialization
for k=1:K
    G(2:end,:,k) = ones(9,2);
end

% estimate graph structure for each class
for k=1:K
    % fill in G(:,:,k)
    % use ConvertAtoG to convert a maximum spanning tree to a graph G
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%
    t_data = dataset(find(labels(:,k)==1),:,:);
    [A W] = LearnGraphStructure(t_data);
    G(:,:,k) = ConvertAtoG(A);
end

% estimate parameters

P.c = zeros(1,K);
% compute P.c

% the following code can be copied from LearnCPDsGivenGraph.m
% with little or no modification
%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_C = sum(labels, 1);
P.c = count_C / sum(count_C);

P.clg = repmat(struct('mu_y', [], 'sigma_y', [], 'mu_x', [], 'sigma_x', [], 'mu_angle', [], 'sigma_angle', [], 'theta', []), 1, 10);

for i=1:10
    for k=1:K
        if G(i,1,k) == 0
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
            par = G(i,2,k);
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

loglikelihood = 0.0;
logPC = log(P.c);
for n=1:N % samples
    sum_k = 0.0;
    for k=1:K % classes
        t = logPC(k);
        for  i=1:10
            % body part i only has the class variable as its parents
            if G(i,1,k) == 0
                t_1 = lognormpdf(dataset(n,i,1), P.clg(i).mu_y(k), P.clg(i).sigma_y(k));
                t_2 = lognormpdf(dataset(n,i,2), P.clg(i).mu_x(k), P.clg(i).sigma_x(k));
                t_3 = lognormpdf(dataset(n,i,3), P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k));
                t = t + t_1 + t_2 + t_3;
                %disp(sprintf('|SAMPLE-%d| class %d, data %d, %f',n, k, i, t_1+t_2+t_3));
            % body part i has, besides class variable, another parent
            % G(i,2)
            else
                par = G(i, 2, k);
                v = [1, dataset(n,par,1), dataset(n,par,2), dataset(n,par,3)];
                t_1 = lognormpdf(dataset(n,i,1), P.clg(i).theta(k,1:4) * v', P.clg(i).sigma_y(k));
                t_2 = lognormpdf(dataset(n,i,2), P.clg(i).theta(k,5:8) * v', P.clg(i).sigma_x(k));
                t_3 = lognormpdf(dataset(n,i,3), P.clg(i).theta(k,9:12) * v', P.clg(i).sigma_angle(k));
                t = t + t_1 + t_2 + t_3;
                %disp(sprintf('|SAMPLE-%d| class %d, data %d, %f',n, k, i, t_1+t_2+t_3));
            end
            % disp(sprintf('|SAMPLE-%d| class %d, data %d, %f',n, k, i, t));
        end
        sum_k = sum_k + exp(t);        
    end    
    loglikelihood = loglikelihood + log(sum_k);
end

fprintf('log likelihood: %f\n', loglikelihood);