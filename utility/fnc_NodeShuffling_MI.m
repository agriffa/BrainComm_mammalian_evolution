function [MIrand] = fnc_NodeShuffling_MI(MI, nrand, data_dir)
%
% Generate null model
% MIrand: time series shuffling across brain regions
%
% Alessandra Griffa
% alessandra.griffa@gmail.com
% University of Geneva
% May 2022
%


rng('default');

% Initialize output
MIrand = zeros([size(MI), nrand]);

% Number of nodes (brain regions)
nn = size(MI,1);

% Number of subjects
ns = size(MI,3);
    
% Loop over randomizations
for r = 1:nrand

    disp(['  ... randomization ' num2str(r) ' of ' num2str(nrand)]);

    % Loop over subjects
    for s = 1:ns
        
        ii = randperm(nn);
        MIrand(:,:,s,r) = MI(ii,ii,s);
        
    end
end

% Save randomized MI data
out_filename = fullfile(data_dir,strcat('MIrand_',num2str(nrand),'nrand_',num2str(nn),'nodes_',num2str(ns),'subj_',date,'.mat'));
save(out_filename,'MIrand','-v7.3');
