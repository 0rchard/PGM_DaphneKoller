function accuracy = ClassifyDataset(dataset, labels, P, G)
% returns the accuracy of the model P and graph G on the dataset 
%
% Inputs:
% dataset: N x 10 x 3, N test instances represented by 10 parts
% labels:  N x 2 true class labels for the instances.
%          labels(i,j)=1 if the ith instance belongs to class j 
% P: struct array model parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description) 
%
% Outputs:
% accuracy: fraction of correctly classified instances (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset,1); % number of examples
K = length(P.c); % number of classes
accuracy = 0.0;
correct = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logPC = log(P.c);
for n=1:N
    prob = zeros(1,K);
    for k=1:K % classes
        t = logPC(k);
        for  i=1:10
            % body part i only has the class variable as its parents
            if G(i,1) == 0
                t_1 = lognormpdf(dataset(n,i,1), P.clg(i).mu_y(k), P.clg(i).sigma_y(k));
                t_2 = lognormpdf(dataset(n,i,2), P.clg(i).mu_x(k), P.clg(i).sigma_x(k));
                t_3 = lognormpdf(dataset(n,i,3), P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k));
                t = t + t_1 + t_2 + t_3;             
            % body part i has, besides class variable, another parent
            % G(i,2)
            else
                par = G(i, 2);
                v = [1, dataset(n,par,1), dataset(n,par,2), dataset(n,par,3)];
                t_1 = lognormpdf(dataset(n,i,1), P.clg(i).theta(k,1:4) * v', P.clg(i).sigma_y(k));
                t_2 = lognormpdf(dataset(n,i,2), P.clg(i).theta(k,5:8) * v', P.clg(i).sigma_x(k));
                t_3 = lognormpdf(dataset(n,i,3), P.clg(i).theta(k,9:12) * v', P.clg(i).sigma_angle(k));
                t = t + t_1 + t_2 + t_3;      
            end            
        end
        prob(k) = exp(t);
    end
    [m idx] = max(prob);    
    if labels(n, idx) == 1
        correct = correct + 1;
    end
end

accuracy = correct / N;


fprintf('Accuracy: %.2f\n', accuracy);