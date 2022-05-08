function [PM] = fnc_computePathMarkovityMatrix(MI,SP)

%
% MI: stack of mutual information matrices
% SP: k-shortest path
% PM: k path markovity matrix (symmetric, binary)
%
% Author
% Alessandra Griffa
% University of Geneva
% May 2022
%

% Number of nodes (brain regions)
nn = size(MI,1);
% Number of considered shortest paths between region pairs (k-shortest paths)
nk = length(SP{1,2});
% Number of subjects
ns = size(MI,3);

% Initialize output
PM = zeros(nn,nn,nk,ns);

% Loop over subjects
for s = 1:ns
    
    thisMI = MI(:,:,s);
    
    % Loop over k (k-shortest paths)
    for k = 1:nk

        temp = zeros(nn);

        % Loop over node pairs
        for i = 1:nn-1
            for j = i+1:nn

                % Node sequence
                nodeSeq = SP{i,j}{k};

                % If path length = 1 hop, continue
                if length(nodeSeq) <= 2
                    continue
                end

                % Assess whether path is Markovian (i.e., DPI respected in both directions)
                temp(i,j) = fnc_is_markovian(nodeSeq, thisMI);

            end
        end

        PM(:,:,k,s) = temp + temp';

    end
end
