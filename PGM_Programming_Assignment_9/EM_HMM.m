% File: EM_HMM.m
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [P loglikelihood ClassProb PairProb] = EM_HMM(actionData, poseData, G, InitialClassProb, InitialPairProb, maxIter)
format long g;
% INPUTS
% actionData: structure holding the actions as described in the PA
% poseData: N x 10 x 3 matrix, where N is number of poses in all actions
% G: graph parameterization as explained in PA description
% InitialClassProb: N x K matrix, initial allocation of the N poses to the K
%   states. InitialClassProb(i,j) is the probability that example i belongs
%   to state j.
%   This is described in more detail in the PA.
% InitialPairProb: V x K^2 matrix, where V is the total number of pose
%   transitions in all HMM action models, and K is the number of states.
%   This is described in more detail in the PA.
% maxIter: max number of iterations to run EM

% OUTPUTS
% P: structure holding the learned parameters as described in the PA
% loglikelihood: #(iterations run) x 1 vector of loglikelihoods stored for
%   each iteration
% ClassProb: N x K matrix of the conditional class probability of the N examples to the
%   K states in the final iteration. ClassProb(i,j) is the probability that
%   example i belongs to state j. This is described in more detail in the PA.
% PairProb: V x K^2 matrix, where V is the total number of pose transitions
%   in all HMM action models, and K is the number of states. This is
%   described in more detail in the PA.

% Initialize variables
N = size(poseData, 1);
K = size(InitialClassProb, 2);
L = size(actionData, 2); % number of actions
V = size(InitialPairProb, 1);

ClassProb = InitialClassProb;
PairProb = InitialPairProb;

loglikelihood = zeros(maxIter,1);

P.c = [];
P.clg.sigma_x = [];
P.clg.sigma_y = [];
P.clg.sigma_angle = [];

% EM algorithm
for iter=1:maxIter

    % M-STEP to estimate parameters for Gaussians
    % Fill in P.c, the initial state prior probability (NOT the class probability as in PA8 and EM_cluster.m)
    % Fill in P.clg for each body part and each class
    % Make sure to choose the right parameterization based on G(i,1)
    % Hint: This part should be similar to your work from PA8 and EM_cluster.m

    P.c = zeros(1,K);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial state prior
    for i=1:L
        P.c = P.c + ClassProb(actionData(i).marg_ind(1),:);
    end
    P.c = P.c / sum(P.c);

    % Emission CPD
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

    % M-STEP to estimate parameters for transition matrix
    % Fill in P.transMatrix, the transition matrix for states
    % P.transMatrix(i,j) is the probability of transitioning from state i to state j
    P.transMatrix = zeros(K,K);

    % Add Dirichlet prior based on size of poseData to avoid 0 probabilities
    P.transMatrix = P.transMatrix + size(PairProb,1) * .05;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:L
        idx = actionData(n).pair_ind;
        for i=1:length(idx)
            P.transMatrix = P.transMatrix + reshape(PairProb(idx(i), :), K, K);
        end
    end
    for r=1:K
        P.transMatrix(r,:) = P.transMatrix(r,:) / sum(P.transMatrix(r,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % E-STEP preparation: compute the emission model factors (emission probabilities) in log space for each
    % of the poses in all actions = log( P(Pose | State) )
    % Hint: This part should be similar to (but NOT the same as) your code in EM_cluster.m

    logEmissionProb = zeros(N,K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:N
        sum_k = 0;
        for k=1:K % classes
            t = 0;
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
            %disp(sprintf('|DEBUG|n:%d, k:%d, t:%f', n, k, t));
            logEmissionProb(n,k) = t;
            sum_k = sum_k + exp(t);
        end
        %ClassProb(n,:) = NormalizeProb(ClassProb(n,:));
        %loglikelihood(iter) = loglikelihood(iter) + log(sum_k);
    end
    loglikelihood(iter) = sum(logsumexp(logEmissionProb));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % E-STEP to compute expected sufficient statistics
    % ClassProb contains the conditional class probabilities for each pose in all actions
    % PairProb contains the expected sufficient statistics for the transition CPDs (pairwise transition probabilities)
    % Also compute log likelihood of dataset for this iteration
    % You should do inference and compute everything in log space, only converting to probability space at the end
    % Hint: You should use the logsumexp() function here to do probability normalization in log space to avoid numerical issues

    %ClassProb = zeros(N,K);
    PairProb = zeros(V,K^2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:L
        singletonFactors = repmat(struct('var', 0, 'card', 0, 'val', []), 1 + length(actionData(i).marg_ind), 1);
        LSF = length(singletonFactors);
        for j=1:LSF
            singletonFactors(j).var = j;
            singletonFactors(j).card = K;
        end
        singletonFactors(1).val = log(P.c);
        for j=2:LSF
            idx = actionData(i).marg_ind(j-1);
            singletonFactors(j).val = log(ClassProb(idx, :));
        end

        pairFactors = repmat(struct('var', 0, 'card', 0, 'val', []), length(actionData(i).pair_ind), 1);
        LPF = length(pairFactors);
        for j=1:LPF
            pairFactors(j).var = [j j+1];
            pairFactors(j).card = [K K];
            pairFactors(j).val = log(reshape(P.transMatrix, 1, K^2));
        end
        Factors = [singletonFactors; pairFactors];
        
        [sFactors, PCalibrated] = ComputeExactMarginalsHMM(Factors);
        
        for j=2:length(sFactors)
            idx = actionData(i).marg_ind(j-1);
            ClassProb(idx,:) = exp(sFactors(j).val - logsumexp(sFactors(j).val));
        end

        for j=1:LPF
            idx = actionData(i).pair_ind(j);
            PairProb(idx,:) = exp(PCalibrated.cliqueList(j).val-logsumexp(PCalibrated.cliqueList(j).val));
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Print out loglikelihood
    disp(sprintf('EM iteration %d: log likelihood: %f', ...
        iter, loglikelihood(iter)));
    if exist('OCTAVE_VERSION')
        fflush(stdout);
    end

    % Check for overfitting by decreasing loglikelihood
    if iter > 1
        if loglikelihood(iter) < loglikelihood(iter-1)
            break;
        end
    end

end

% Remove iterations if we exited early
loglikelihood = loglikelihood(1:iter);