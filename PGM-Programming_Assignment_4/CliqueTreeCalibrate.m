%CLIQUETREECALIBRATE Performs sum-product or max-product algorithm for 
%clique tree calibration.

%   P = CLIQUETREECALIBRATE(P, isMax) calibrates a given clique tree, P 
%   according to the value of isMax flag. If isMax is 1, it uses max-sum
%   message passing, otherwise uses sum-product. This function 
%   returns the clique tree where the .val for each clique in .cliqueList
%   is set to the final calibrated potentials.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function P = CliqueTreeCalibrate(P, isMax)


% Number of cliques in the tree.
N = length(P.cliqueList);

% Setting up the messages that will be passed.
% MESSAGES(i,j) represents the message going from clique i to clique j. 
MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We have split the coding part for this function in two chunks with
% specific comments. This will make implementation much easier.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
% While there are ready cliques to pass messages between, keep passing
% messages. Use GetNextCliques to find cliques to pass messages between.
% Once you have clique i that is ready to send message to clique
% j, compute the message and put it in MESSAGES(i,j).
% Remember that you only need an upward pass and a downward pass.
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug = 0;
if debug == 1
    disp(sprintf('|DEBUG|MESSAGE PASSING'));
    disp(sprintf('--------------------------------------------------------------------'));
end

if isMax == 1
    for i=1:N
        P.cliqueList(i).val = log(P.cliqueList(i).val);
    end
end

while 1
    [i,j] = GetNextCliques(P, MESSAGES);
    if [i,j] == [0 0]
        break;
    end
    if debug == 1        
        disp(sprintf('|NextClique| %d->%d',i,j));
    end
    % incoming message to i except j    
    factor_i = P.cliqueList(i);
    
    for k=1:N
        if i~=k && j~=k && ~isempty(MESSAGES(k,i).var)
            if isMax == 1
                factor_i = FactorSum(factor_i, MESSAGES(k,i));
            else
                factor_i = FactorProduct(factor_i, MESSAGES(k,i));
            end
        end
    end
    % variables to be marginalized
    var_e = setdiff(P.cliqueList(i).var, P.cliqueList(j).var);
    if isMax == 1
        factor_mar = FactorMaxMarginalization(factor_i, var_e);
        MESSAGES(i,j) = factor_mar;
    else
        factor_mar = FactorMarginalization(factor_i, var_e);
        % normalize after factor marginalization
        factor_norm = FactorNormalization(factor_mar);    
        MESSAGES(i,j) = factor_norm;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
% Now the clique tree has been calibrated. 
% Compute the final potentials for the cliques and place them in P.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    % incoming message to i
    for k=1:N
        if ~isempty(MESSAGES(k,i).var)
            if isMax == 1
                P.cliqueList(i) = FactorSum(P.cliqueList(i), MESSAGES(k,i));
            else
                P.cliqueList(i) = FactorProduct(P.cliqueList(i), MESSAGES(k,i));
            end
        end
    end
end

return
