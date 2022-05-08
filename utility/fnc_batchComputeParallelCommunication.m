function [PM,PCS] = fnc_batchComputeParallelCommunication(data_dir, mi_prefix, weis) 
%
% Compute Path Markovity matrices
%
% Input
% data_dir: data directory (e.g., 'DATA/hHCP')
% mi_prefix: prefix of mutual information mat file ('MI_' or 'MIrand_')
% weis: cell array to select the weight that were used for the k-shortest
%       path computation (e.g., {'mm'}) 
%
% Output (output variables are saved to data_dir)
% PM: [nnode X nnode X kmax X ns]
%     one binary symmetric matrix per subject, per k value
%     entry(i,j)=1: the multi-step path connecting regions i and j is Markovian / supports realy communication
%                   (i.e., the DPI is repsected in both directions) 
%     entry(i,j)=0: the multi-step path connecting regions i and j is NOT Markovian / does NOT support relay communication,
%                   or reiongs i and j are connected through a single-step path
% PCS: [nnode X nnode X ns]
%       values are comprised between 0 and nk;
%       for a single subject, the PCS matrix is obtained by summing
%       his/her nk path markovity (PM) matrices obtained by considergin the
%       mutual information values along the the nk shortest paths
%       (k-shortestest paths) connecting region pairs, according to the
%       data processing inequality (DPI)
%
% Author
% Alessandra Griffa
% University of Geneva
% May 2022
%



% Loop over SC weigths (e.g., {'mm'})
for b = 1:length(weis)
    wei = weis{b};

    % Load data
    disp(wei);

    % MI
    temp = dir(fullfile(data_dir,[mi_prefix,'*.mat']));
    if length(temp) ~= 1
        warning(['WRONG MI mat FILE IN ' data_dir]);
        return
    end
    valMI = load(fullfile(data_dir,temp.name));
    if isfield(valMI,'MIrand')
        valMI.MI = valMI.MIrand;
    end
    
    % KSP
    temp = dir(fullfile(data_dir,['KSP_*SC',wei,'*.mat']));
    if length(temp) ~= 1
        warning(['WRONG SP mat FILE IN ' data_dir]);
        return
    end
    valKSP = load(fullfile(data_dir,temp.name));

    % Number of nodes (brain regions)
    nn = size(valMI.MI,1);
    % Number of subjects
    ns = size(valMI.MI,3);
    % Number of randomizations
    nrand = size(valMI.MI,4);
    % Number of considered shortest paths (k-shortest paths)
    nk = length(valKSP.SP{1,2});

    % Path Markovity and parallel communication score matrices
    PM = nan(nn, nn, nk, ns, nrand);
    PCS = nan(nn, nn, ns, nrand);
    % Loop over randomization
    for r = 1:nrand
        disp([' > loop ' num2str(r) ' of ' num2str(nrand)]);
        thisMI = valMI.MI(:,:,:,r);
        thisPM = fnc_computePathMarkovityMatrix(thisMI,valKSP.SP);
        PCS(:,:,:,r) = squeeze(sum(thisPM,3));
        PM(:,:,:,:,r) = thisPM;
    end
    PM = squeeze(PM);
    PCS = squeeze(PCS);

    % Save outputs
    if nrand > 1
        filename = ['PMrand_',num2str(nn),'nodes_','_SC',wei,'_',date,'.mat'];
        save(fullfile(data_dir,filename),'PM','-v7.3');
        filename = ['PCSrand_',num2str(nn),'nodes_','_SC',wei,'_',date,'.mat'];
        save(fullfile(data_dir,filename),'PCS','-v7.3');
    else
        filename = ['PM_',num2str(nn),'nodes_',num2str(ns),'subj_SC',wei,'_',date,'.mat'];
        save(fullfile(data_dir,filename),'PM','subj','-v7.3');
        filename = ['PCS_',num2str(nn),'nodes_',num2str(ns),'subj_SC',wei,'_',date,'.mat'];
        save(fullfile(data_dir,filename),'PCS','subj','-v7.3');
    end

end
           
 
        
        

