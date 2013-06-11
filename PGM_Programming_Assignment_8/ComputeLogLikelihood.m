function loglikelihood = ComputeLogLikelihood(P, G, dataset)
% returns the (natural) log-likelihood of data given the model and graph structure
%
% Inputs:
% P: struct array parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description)
%
%    NOTICE that G could be either 10x2 (same graph shared by all classes)
%    or 10x2x2 (each class has its own graph). your code should compute
%    the log-likelihood using the right graph.
%
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% 
% Output:
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset,1); % number of examples
K = length(P.c); % number of classes

loglikelihood = 0;
% You should compute the log likelihood of data as in eq. (12) and (13)
% in the PA description
% Hint: Use lognormpdf instead of log(normpdf) to prevent underflow.
%       You may use log(sum(exp(logProb))) to do addition in the original
%       space, sum(Prob).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logPC = log(P.c);
for n=1:N % samples
    sum_k = 0.0;
    for k=1:K % classes
        t = logPC(k);
        for  i=1:10
            % body part i only has the class variable as its parents
            if G(i,1) == 0
                t_1 = lognormpdf(dataset(n,i,1), P.clg(i).mu_y(k), P.clg(i).sigma_y(k));
                t_2 = lognormpdf(dataset(n,i,2), P.clg(i).mu_x(k), P.clg(i).sigma_x(k));
                t_3 = lognormpdf(dataset(n,i,3), P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k));
                t = t + t_1 + t_2 + t_3;
                %disp(sprintf('|SAMPLE-%d| class %d, data %d, %f',n, k, i, t_1+t_2+t_3));
            % body part i has, besides class variable, another parent
            % G(i,2)
            else
                par = G(i, 2);
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