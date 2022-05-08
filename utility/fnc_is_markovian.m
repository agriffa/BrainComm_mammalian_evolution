function val = fnc_is_markovian(nodeSeq, MI)

%
% Determine whether a nodes' sequence corresponds to a Markovian path, 
% i.e., the mutual information values with respect to the source node are
% monotonically non-increaseing along the paths. This condition must be 
% respected in both directions of the nodes' sequence. A Markovian path
% supports relay communication between a region pair.
%
% nodeSeq: nodes sequence (required: length(nodeSeq > 2))
% MI: mutual information matrix
% val: 0: the nodes' sequence does not correspond to a Markovian path
%      1: the nodes' sequence corresponds to a Markovian path
%
% Author
% Alessandra Griffa
% University of Geneva
% May 2022
%


% MI values along nodes' sequences, forward and backward directions
MI_fwd = MI(nodeSeq(1), nodeSeq(2:1:end));
MI_bwd = MI(nodeSeq(end), nodeSeq(end-1:-1:1));

% Test whether MI values are monotonically non-increasing
val1 = issorted(MI_fwd, 'descend');
val2 = issorted(MI_bwd, 'descend');

% Output 
val = min([val1, val2]);
