% File: EM_cluster.m
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [P loglikelihood ClassProb] = EM_cluster(poseData, G, InitialClassProb, maxIter)

% INPUTS
% poseData: N x 10 x 3 matrix, where N is number of poses;
%   poseData(i,:,:) yields the 10x3 matrix for pose i.
% G: graph parameterization as explained in PA8
% InitialClassProb: N x K, initial allocation of the N poses to the K
%   classes. InitialClassProb(i,j) is the probability that example i belongs
%   to class j
% maxIter: max number of iterations to run EM

% OUTPUTS
% P: structure holding the learned parameters as described in the PA
% loglikelihood: #(iterations run) x 1 vector of loglikelihoods stored for
%   each iteration
% ClassProb: N x K, conditional class probability of the N examples to the
%   K classes in the final iteration. ClassProb(i,j) is the probability that
%   example i belongs to class j

% Initialize variables
N = size(poseData, 1);
K = size(InitialClassProb, 2);

ClassProb = InitialClassProb;

loglikelihood = zeros(maxIter,1);

P.c = [];
P.clg.sigma_x = [];
P.clg.sigma_y = [];
P.clg.sigma_angle = [];

% EM algorithm
for iter=1:maxIter

    % M-STEP to estimate parameters for Gaussians
    %
    % Fill in P.c with the estimates for prior class probabilities
    % Fill in P.clg for each body part and each class
    % Make sure to choose the right parameterization based on G(i,1)
    %
    % Hint: This part should be similar to your work from PA8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P.c = zeros(1,K);
    P.c = sum(ClassProb, 1) / N;
    P.clg = repmat(struct('mu_y', [], 'sigma_y', [], 'mu_x', [], 'sigma_x', [], 'mu_angle', [], 'sigma_angle', [], 'theta', []), 1, 10);

    for i=1:10
        for k=1:K
            if G(i,1) == 0
                U = squeeze(poseData(:,i,:));
                [P.clg(i).mu_y(k), P.clg(i).sigma_y(k)] = FitG(U(:,1), ClassProb(:,k));
                [P.clg(i).mu_x(k), P.clg(i).sigma_x(k)] = FitG(U(:,2), ClassProb(:,k));
                [P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k)] = FitG(U(:,3), ClassProb(:,k));
            else
                par = G(i,2);
                U_y = squeeze(poseData(:,par,:));
                X_y = squeeze(poseData(:,i,1));                
                [theta_y sigma_y] = FitLG(X_y, U_y, ClassProb(:,k));
                P.clg(i).theta(k,1:4) = [theta_y(4), theta_y(1), theta_y(2), theta_y(3)];
                P.clg(i).sigma_y(k) = sigma_y;

                U_x = squeeze(poseData(:,par,:));
                X_x = squeeze(poseData(:,i,2));                
                [theta_x sigma_x] = FitLG(X_x, U_x, ClassProb(:,k));
                P.clg(i).theta(k,5:8) = [theta_x(4), theta_x(1), theta_x(2), theta_x(3)];
                P.clg(i).sigma_x(k) = sigma_x;

                U_angle = squeeze(poseData(:,par,:));
                X_angle = squeeze(poseData(:,i,3));
                [theta_angle sigma_angle] = FitLG(X_angle, U_angle, ClassProb(:,k));
                P.clg(i).theta(k,9:12) = [theta_angle(4), theta_angle(1), theta_angle(2), theta_angle(3)];
                P.clg(i).sigma_angle(k) = sigma_angle;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % E-STEP to re-estimate ClassProb using the new parameters
    %
    % Update ClassProb with the new conditional class probabilities.
    % Recall that ClassProb(i,j) is the probability that example i belongs to
    % class j.
    %
    % You should compute everything in log space, and only convert to
    % probability space at the end.
    %
    % Tip: To make things faster, try to reduce the number of calls to
    % lognormpdf, and inline the function (i.e., copy the lognormpdf code
    % into this file)
    %
    % Hint: You should use the logsumexp() function here to do
    % probability normalization in log space to avoid numerical issues

    ClassProb = zeros(N,K);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loglikelihood(iter) = 0.0;
    logPC = log(P.c);
    temp = 0.0 * zeros(N, K);
    for n=1:N
        for k=1:K % classes
            t = logPC(k);
            for  i=1:10
                % body part i only has the class variable as its parents
                if G(i,1) == 0
                    t_1 = lognormpdf(poseData(n,i,1), P.clg(i).mu_y(k), P.clg(i).sigma_y(k));
                    t_2 = lognormpdf(poseData(n,i,2), P.clg(i).mu_x(k), P.clg(i).sigma_x(k));
                    t_3 = lognormpdf(poseData(n,i,3), P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k));
                    t = t + t_1 + t_2 + t_3;
                    % body part i has, besides class variable, another parent
                    % G(i,2)
                else
                    par = G(i, 2);
                    v = [1, poseData(n,par,1), poseData(n,par,2), poseData(n,par,3)];
                    t_1 = lognormpdf(poseData(n,i,1), P.clg(i).theta(k,1:4) * v', P.clg(i).sigma_y(k));
                    t_2 = lognormpdf(poseData(n,i,2), P.clg(i).theta(k,5:8) * v', P.clg(i).sigma_x(k));
                    t_3 = lognormpdf(poseData(n,i,3), P.clg(i).theta(k,9:12) * v', P.clg(i).sigma_angle(k));
                    t = t + t_1 + t_2 + t_3;
                end
            end
            ClassProb(n,k) = exp(t);     
            temp(n,k) = t;
        end
        ClassProb(n,:) = NormalizeProb(ClassProb(n,:));
    end
    loglikelihood(iter) = sum(logsumexp(temp));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute log likelihood of dataset for this iteration
    % Hint: You should use the logsumexp() function here
    %loglikelihood(iter) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Print out loglikelihood
    disp(sprintf('EM iteration %d: log likelihood: %f', ...
        iter, loglikelihood(iter)));
    if exist('OCTAVE_VERSION')
        fflush(stdout);
    end

    % Check for overfitting: when loglikelihood decreases
    if iter > 1
        if loglikelihood(iter) < loglikelihood(iter-1)
            break;
        end
    end
end

% Remove iterations if we exited early
loglikelihood = loglikelihood(1:iter);
