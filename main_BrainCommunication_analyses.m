%
% This code reproduce the main results of the article 'The evolution of 
% information transmission in mammalian brain networks' - Alessandra 
% Griffa, Mathieu Mach, Julien Dedelley, Daniel Gutierrez-Barragan, 
% Alessandro Gozzi, Gilles Allali, Joanes Grandjean, Dimitri Van De Ville,
% Enrico Amico.
%
%
% Alessandra Griffa
% University of Geneva
% May 2022
%


clear
close all
clc

%
% -> CHANGE MATLAB PATH TO BrainComm_mammalian_evolution LOCAL REPOSITORY
%
% -> SET DATA DIRECTORY
%    Data to reproduce the results of Griffa et al. (202x) can be
%    downloadad at ZENODO REPOSITORY
data_dir = '/Users/alli/Desktop/PROJECT_IT/DATASET_BrainComm_mammalian_evolution';



%% Preliminary settings
% Add directories to Matlab path
gitdir = pwd;
addpath(genpath(fullfile(gitdir,'utility')));         

% Set plotkit directory (for cortical surface renderings)
plotkit_dir = fullfile(gitdir, 'utility', 'surface_plotkit');

% Structural connectivity weigth 
wei = 'mm';              % Euclidean distance (in millimiters)

% LOAD MULTIPLE DATASETS
% -------------------------------------------------------------------------
% NOTE: the first letter of 'spec_string' and 'spec' entries indicates the
% species (h = human; q = macaque; m = mouse)
spec_string = {'hHCP','qTVB','mAD3'};
spec = {'h','q','m'};                           
        
% Number of subecjts for each dataset
% hHCP: up to 100; qTVB: 9; mAD3: 10
pii{1} = 1:100;     % h HCP 
pii{2} = 1:9;       % q TVB
pii{3} = 1:10;      % m AD3

% k (k-shortest path - number of considered shortest path per region pair)
nk = 5;     

% Colors
hcol = [25, 118, 210]./255;     % human: blue
qcol = [0, 137, 123]./255;      % macaque: green
mcol = [216, 67, 21]./255;      % mouse: orange
rcol = [207, 216, 220]./255;    % random: grey
hmap = 'Blues9';                % human colormap
qmap = 'Greens9';               % macaque colormap
mmap = 'Reds9';                 % mouse colormap



%% LOAD parcellation information
% Human
load(fullfile(data_dir, 'parcellations', 'hinfo100_with_rsn.mat'),'hinfo');

% Macaque
load(fullfile(data_dir, 'parcellations', 'qinfo82_with_rsn.mat'),'qinfo');

% Mouse
load(fullfile(data_dir, 'parcellations', 'minfo78_with_rsn.mat'),'minfo');



%% LOAD data 
%
% nn = number of nodes (brain regions)
% ns = number of subejcts
% nk = number of considered shortest paths (k-shortest paths)
% nrand = number of randomizations
%
% STRUCTURAL CONNECTIVITY INFORMATION
% SP: k-shortest path node sequence
% SPhops: [nn X nn X nk] k-shortest path 'binary' length (number of hops) 
% SPlen: [nn X nn X nk] k-shortest path 'weighted' length (sum of connection-length weights)
% SC: [nn X nn] group-representiative structural connectivity matrix (conneciton-length weights -> Euclidean distance (mm))
% 
% FUNCTIONAL INFORMATION
% MI: [nn X nn X ns] individual mutual information matricies (symmetric)
% ntp: number of time points of original time series
% TR: temporal resolution of original time series (ms)
% subj: experimental subjects' IDs
%
% RANDOMIZED FUNCTIONAL DATA
% MIrand: [nn X nn X ns X nrand] randomized mutual information matricies (symmetric)
%         A null model was defined by randomly shuffling the raw fMRI time series across
%         brain regions while preserving the original structural connectivity information.
%         Note that this randomization preserves the statistical properties of both
%         the original functional and structural data, since we are merely rearranging spatially 
%         fMRI time series across the brain network.
% subj: experimental subjects' IDs
%

% Loop over datasets
for s = 1:length(spec_string)
    
    disp('');
    disp('');
    disp(['    ' spec_string{s} ' - dataset ' num2str(s) ' of ' num2str(length(spec_string)) ' - LOAD DATA']);
    
    % Assign color and parcellation info according to species
    if spec_string{s}(1) == 'h'
        eval([spec{s},'col = hcol;']); 
        eval([spec{s},'map = hmap;']); 
        eval([spec{s},'info = hinfo;']);
    elseif spec_string{s}(1) == 'q'
        eval([spec{s},'col = qcol;']); 
        eval([spec{s},'map = qmap;']);         
        eval([spec{s},'info = qinfo;']);
    elseif spec_string{s}(1) == 'm'
        eval([spec{s},'col = mcol;']); 
        eval([spec{s},'map = mmap;']);
        if contains(spec_string{s},'GG')
            eval([spec{s},'info = m66info;']);
        else
            eval([spec{s},'info = minfo;']);
        end
    else
        warning('Wrong species string definition');
        continue
    end
    
    this_data_dir = fullfile(data_dir, 'DATA', spec_string{s});
    
    % LOAD STRUCTURAL CONNECTIVITY INFORMATION
    filename = dir(fullfile(this_data_dir,['KSP_*',wei,'*.mat']));
    if length(filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - KSP FILE']);
        continue
    else
        a = load(fullfile(this_data_dir,filename.name));
        eval([spec{s},'SP = a.SP; ',spec{s},'SPhops = a.SPhops; ',spec{s},'SPlen = a.SPlen;']); 
        eval([spec{s},'SC = a.SC' wei ';']);
        eval([spec{s},'SC(isinf(',spec{s},'SC)) = 0;']);
    end
    
    % LOAD FUNCTIONAL INFORMATION
    filename = dir(fullfile(this_data_dir,'MI_*.mat'));
    if length(filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - MI FILE']);
        continue
    else
        a = load(fullfile(this_data_dir,filename.name));
        eval([spec{s},'MI = a.MI(:,:,pii{s});']);
    end
    
    % LOAD RANDOMIZED FUNCTIONAL DATA
    filename = dir(fullfile(this_data_dir,'MIrand_*.mat'));
    if length(filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - MI FILE']);
        continue
    else
        a = load(fullfile(this_data_dir,filename.name));
        eval([spec{s},'MIrand = a.MIrand(:,:,:,pii{s});']);
    end    
    
    clear a filename this_data_dir
    
end



%% VISUALIZE INPUT DATA
close all
% Loop over datasets
for s = 1:length(spec)   
    
    eval(['thisSC = ',spec{s},'SC;']);
    eval(['thisSPhops = ',spec{s},'SPhops;']);
    eval(['thisInfo = ',spec{s},'info;']);
    nn = size(thisSC,1);
    d(s) = sum(thisSC(:) > 0) / (nn * (nn-1));
    
    figure, ax(1) = subplot(2,2,1);
    imagesc(thisSC(thisInfo.rsnOrder,thisInfo.rsnOrder)), colormap(ax(1), jet), colorbar, axis equal tight, fnc_setfig;
    fnc_addYeoLines(thisInfo, [0, 0, 0], 2), title([spec_string{s} ' SC' wei ' (length), d = ' num2str(d(s),2)]);
    
    % Percentage of 2-step and multi-step paths as a function of k
    eval(['thisK = ',spec{s},'SPhops;']);
    val_1step = zeros(nk,1);
    val_2step = zeros(nk,1);
    val_MULTIstep = zeros(nk,1);
    for k = 1:nk
        temp = squeeze(thisK(:,:,k));
        a = zeros(size(temp)); a(temp > 2) = 1;
        b = zeros(size(temp)); b(temp == 2) = 1;
        c = zeros(size(temp)); c(temp == 1) = 1;
        val_MULTIstep(k) = sum(a(:)) / (nn * (nn-1));
        val_2step(k) = sum(b(:)) / (nn * (nn-1));
        val_1step(k) = sum(c(:)) / (nn * (nn-1));
    end
    subplot(2,2,2), bar([val_2step, val_MULTIstep, val_1step],'stacked'), xlabel('k'), ylabel('paths percentage'), ylim([0 1]);
    legend('2-step paths','longer paths','direct','Location','northwest','fontsize',20), title(spec_string{s}), fnc_setfig; 
    
    % MI
    eval(['thisMI = ',spec{s},'MI;']);
    val = mean(thisMI,3);
    minval = min(val(:));
    maxval = max(val(:));
    ax(2) = subplot(2,2,3);
    imagesc(val(thisInfo.rsnOrder,thisInfo.rsnOrder)), colormap(ax(2),othercolor('Purples8')), colorbar, axis equal tight, fnc_setfig;
    fnc_addYeoLines(thisInfo, [0, 0, 0], 2), caxis([minval, maxval]), title('group-average MI');
    
    % MI null model
    eval(['thisMI = ',spec{s},'MIrand;']);
    val = mean(mean(thisMI,4),3);
    ax(3) = subplot(2,2,4);
    imagesc(val(thisInfo.rsnOrder,thisInfo.rsnOrder)), colormap(ax(3), othercolor('Purples8')), colorbar, axis equal tight, fnc_setfig;
    fnc_addYeoLines(thisInfo, [0, 0, 0], 2), caxis([minval, maxval]), title([spec_string{s} ' null model MI (average)']);    
     
end



%% Compute (or load) PARALLEL COMMUNICATION INFORMATION
%
% nn = number of nodes (brain regions)
% ns = number of subejcts
% nk = number of considered shortest paths (k-shortest paths)
%
% INDIVIDUAL PATH MARKOVITY MATRICES
% PM: [nn X nn X nk X ns]
%     For each subject, nk binary and symmetric communication matrices 
%     which indicate whether the k-shortest path between a region pair is 
%     or is not a relay (markovian) communication pathway
%
% INDIVIDUAL PARALLEL COMMUNICATION SCORE (PCS) MATRICES
% PM: [nn X nn X ns]
%

% Loop over datasets
for s = 1:length(spec_string)
    
    disp('');
    disp('');
    disp(['    ' spec_string{s} ' - dataset ' num2str(s) ' of ' num2str(length(spec_string)) ' - PARALLEL COMMUNICATION INFORMAITON']);

    this_data_dir = fullfile(data_dir, 'DATA', spec_string{s});
    
    % LOAD PATH MARKOVITY MATRICES AND PARALLEL COMMUNICATION MATRICES
    filenamePM = dir(fullfile(this_data_dir,['PM_*',wei,'*.mat']));
    filenamePCS = dir(fullfile(this_data_dir,['PCS_*',wei,'*.mat']));
    filenamePMrand = dir(fullfile(this_data_dir,['PMrand_*',wei,'*.mat']));
    filenamePCSrandavg = dir(fullfile(this_data_dir,['PCSrandavg_*',wei,'*.mat']));
    
    % DATA
    % If mat data does not exist, compute parallel communicaiton information
    if ( length(filenamePM) ~= 1 ) || ( length(filenamePCS) ~= 1 )
        [thisPM, thisPCS] = fnc_batchComputeParallelCommunication(this_data_dir, 'MI_', {wei});
        eval([spec{s},'PM = thisPM(:,:,:,pii{s});']);
        eval([spec{s},'PCS = thisPCS(:,:,pii{s});']);
    % else, load mat data
    else
        a = load(fullfile(this_data_dir,filenamePM.name));
        eval([spec{s},'PM = a.PM(:,:,:,pii{s});']);
        a = load(fullfile(this_data_dir,filenamePCS.name));
        eval([spec{s},'PCS = a.PCS(:,:,pii{s});']);
    end
    
    % RANDOMIZED DATA - individual path markovity and parallel
    % communication score matrices
    % If mat data does not exist, compute parallel communicaiton information
    if length(filenamePMrand) ~= 1
        [thisPM, ~] = fnc_batchComputeParallelCommunication(this_data_dir, 'MIrand_', {wei});
        eval([spec{s},'PMrand = thisPM(:,:,:,pii{s});']);
    % else, load mat data
    else
        a = load(fullfile(this_data_dir,filenamePMrand.name));
        eval([spec{s},'PMrand = a.PM(:,:,:,pii{s},:);']);
    end
    
    % RANDOMIZED DATA - null model for group-average parallel communication score matrices
    % This null model was obtained by repating the randomization
    % procedure 3000 time over all subjects; for each randomization (over
    % all subejcts), a group-average randomized PCS matrix was stored.
    % NOTE: only 10 human subjects (instead of 100) were used for the null
    % model and related analyses because the null distribution's standard 
    % deviation (and therefore the z-scores) are affected by the number of 
    % subjects, which we kept as similar as possible across species for
    % these analyses (i.e., 10 humans, 9 macaques, 10 mice).
    if length(filenamePCSrandavg) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - PCSrandavg file !']);
        continue
    else
        a = load(fullfile(this_data_dir,filenamePCSrandavg.name));
        eval([spec{s},'PCSravg = a.PCSravg;']);
    end
    
    clear a filenamePM filenamePCS this_data_dir
    
end



%% FIGURE 2a
clc
for s = 1:length(spec)

    disp(['    ' spec_string{s} ' - dataset ' num2str(s) ' of ' num2str(length(spec_string)) ' - Fig.2a']);
    
    eval(['PM = ', spec{s}, 'PM;']);
    eval(['SPhops = ', spec{s}, 'SPhops;']);
    eval(['PMrand = ', spec{s}, 'PMrand;']);
    nn = size(PM,1);
    nk = size(PM,3);
    ns = size(PM,4);
    nrand = size(PMrand,5);
    eval(['thiscol = ',spec{s},'col;']);
    
    % Boxplot whole-brain density of relay communication channels for k = 1 to 5
    % (one point per subject or per randomization)
    % Real data
    temp = squeeze(sum(PM,1));
    temp = squeeze(sum(temp,1));

    % Randomized data
    temprand = squeeze(sum(PMrand,1));
    temprand = squeeze(sum(temprand,1));

    % Normalize with respect to number of polysynaptic structural paths
    for k = 1:nk
        temp(k,:) = temp(k,:) ./ sum(sum(SPhops(:,:,k) > 1));
        temprand(k,:,:) = temprand(k,:,:) ./ sum(sum(SPhops(:,:,k) > 1));
    end
    
    % Build boxplots
    A = nan(ns*nrand, nk * 2);
    p = nan(nk,1);
    disp(spec_string{s});
    for k = 1:nk
        A(1:ns, 2*(k-1)+1) = temp(k,:);
        tempr = squeeze(temprand(k,:,:));
        A(1:ns*nrand, 2*(k-1)+2) = tempr(:);
        p(k) = ranksum(temp(k,:),tempr(:));
        disp(['Average (std) density, k = ' num2str(k) ': ' num2str(mean(temp(k,:)))...
            ' (' num2str(std(temp(k,:))) ')']);
        disp(['    comparison with randomized data: p = ' num2str(p(k))]);
    end
    val = temp([2,3,4,5],:);
    disp(['Average (std) density, k = 2,3,4,5: ' num2str(mean(val(:)))...
            ' (' num2str(std(val(:))) ')']);
        
    figure('Name',spec_string{s});
    notBoxPlot_dotcolor(A(:,1:2:nk*2),[1:4:nk*4],0.4,'line',thiscol), hold on;
    notBoxPlot_dotcolor(A(:,2:2:nk*2),[2:4:nk*4],0.4,'line',rcol), hold on;
    xticks([1.5:4:nk*4]), xticklabels({'1st s.p.','2nd s.p.','3rd s.p.','4th s.p.','5th s.p.'}); 
    xlim([0 nk*4]), ylim([0 0.6]), ylabel('density of relay communication paths'), grid on;
    title([num2str(ns) ' ' spec_string{s} ' participants']); fnc_setfig;
    
    disp(' ');

end



%% FIGURE 2b
clc
perc5_human = 0.130; 
perc95_human = 3.023; 
medianPCS = nan(length(spec),1);
perc5PCS = nan(length(spec),1);
perc95PCS = nan(length(spec),1);
for s = 1:length(spec)

    disp(['Species ' num2str(s) ' of ' num2str(length(spec))]);
    
    % PCS [nn X nn X ns]
    eval(['PCS = ', spec{s}, 'PCS;']);
    nn = size(PCS,1);
    ns = size(PCS,3);
    ii = find(triu(ones(nn),1));
    
    % Group-average PCS
    PCSavg = nanmean(PCS,3);
    disp(spec_string{s});
    disp(['  median PCS = ' num2str(median(PCSavg(ii)))]);
    disp(['  mean PCS = ' num2str(mean(PCSavg(ii)))]);
    disp(['  5-percentile PCS = ' num2str(prctile(PCSavg(ii),5))]);
    disp(['  95-percentile PCS = ' num2str(prctile(PCSavg(ii),95))]);
    disp(['  minimum PCS = ' num2str(min(PCSavg(ii)))]);
    disp(['  maximum PCS = ' num2str(max(PCSavg(ii)))]);
    medianPCS(s) = median(PCSavg(ii));
    perc5PCS(s) = prctile(PCSavg(ii),5);
    perc95PCS(s) = prctile(PCSavg(ii),95);
    
    % Plot
    figure, subplot(1,2,1), imagesc(PCSavg), axis equal tight;
    colormap(othercolor('YlOrBr9')), colorbar, caxis([perc5_human perc95_human]), fnc_setfig;
    title({['parallel communication score - ' spec_string{s}],['median = ' num2str(median(PCSavg(find(triu(ones(nn),1)))),3)]});
    subplot(1,2,2), histogram(PCSavg(ii),[0:0.25:5],'normalization','probability','facealpha',0.6,'edgecolor','k','facecolor','w','linewidth',2), hold on, fnc_setfig;
    set(gcf,'position',[38 586 897 391]);
    
end

% Two-sample Kolmogorov-Smirnov test between PCS histograms
hval = nanmean(hPCS,3);
hval = hval(find(triu(ones(size(hPCS,1)),1)));
qval = nanmean(qPCS,3);
qval = qval(find(triu(ones(size(qPCS,1)),1)));
mval = nanmean(mPCS,3);
mval = mval(find(triu(ones(size(mPCS,1)),1)));
[h,p,ks2stat] = kstest2(hval,qval)
[h,p,ks2stat] = kstest2(hval,mval)
[h,p,ks2stat] = kstest2(qval,mval)



%% FIGUREs S6, S7, S8
clc
perc5_human = 0.130; 
perc95_human = 3.023; 
medianPCS_Z = nan(length(spec),1);
perc5PCS_Z = nan(length(spec),1);
perc95PCS_Z = nan(length(spec),1);
for s = 1:length(spec)

    disp(['Species ' num2str(s) ' of ' num2str(length(spec))]);
    
    % PCS [nn X nn X ns]
    eval(['PCS = ', spec{s}, 'PCS;']);
    eval(['PCSravg = ', spec{s}, 'PCSravg;']);
    eval(['rsninfo = ', spec{s}, 'info;']);
    nn = size(PCS,1);
    ns = size(PCS,3);
    ii = find(triu(ones(nn),1));
    
    % Group-average PCS
    PCSavg = nanmean(PCS,3);
    PCSavg10 = nanmean(PCS(:,:,1:min(size(PCS,3),10)),3);
    
    % Plot real and randomized data
    figure, subplot(1,3,1), imagesc(PCSavg), axis equal tight;
    colormap(othercolor('YlOrBr9')), colorbar, caxis([perc5_human perc95_human]), fnc_setfig;
    title({['parallel communication score - ' spec_string{s}],['median = ' num2str(median(PCSavg(find(triu(ones(nn),1)))),3)]});
    subplot(1,3,3), histogram(PCSavg(ii),[0:0.25:5],'normalization','probability','facealpha',0.6), hold on, fnc_setfig;
    thisMean = nanmean(PCSravg,3);
    subplot(1,3,2), imagesc(thisMean), axis equal tight;
    colormap(othercolor('YlOrBr9')), colorbar, caxis([perc5_human perc95_human]), fnc_setfig;
    title({['parallel communication score - ' spec_string{s}],['median = ' num2str(median(PCSravg(:)),3)]});
    subplot(1,3,3), histogram(thisMean(ii),[0:0.25:5],'normalization','probability','facealpha',0.3), xlim([0 4]), ylabel('probability'); 
    set(gcf,'position',[977 710 1523 416]);
    
    % Z-scores and p-values for group-average PCS scores
    thisMean = mean(PCSravg,3);
    thisStd = std(PCSravg,[],3);
    PCSz = ( PCSavg10 - thisMean ) ./ thisStd;
    maskZ = zeros(nn); maskZ(PCSz > 1.96) = 1;
    pval = nan(nn);
    for i = 1:nn-1
        for j = i+1:nn
            pval(i,j) = sum( squeeze(PCSravg(i,j,:)) >= PCSavg10(i,j) ) / size(PCSravg,3);
            pval(i,j) = max(pval(i,j), 1/(size(PCSravg,3)+1));
        end
    end
    adj_p = mafdr(pval(ii),'BHFDR',true);
    maskFDR = zeros(nn); maskFDR(ii(adj_p < 0.05)) = 1; maskFDR = maskFDR + maskFDR';
    eval([spec{s} 'maskZ = maskZ;']);
    eval([spec{s} 'maskFDR = maskFDR;']);
    
    % Plot
    figure, ax(1) = subplot(2,3,1); imagesc(PCSz(rsninfo.rsnOrder,rsninfo.rsnOrder)), axis equal tight, axis off;
    title([spec_string{s} ': PCS z-scores']), colormap(ax(1),jet), colorbar, fnc_addYeoLines(rsninfo, [0, 0, 0], 2), fnc_setfig;
    temp = PCSavg(rsninfo.rsnOrder,rsninfo.rsnOrder) .* maskZ(rsninfo.rsnOrder,rsninfo.rsnOrder);
    ax(2) = subplot(2,3,2); imagesc(temp), axis equal tight, axis off, title('PCS z>1.96');
    colormap(ax(2),othercolor('YlOrBr9')), colorbar, caxis([perc5_human perc95_human]), fnc_addYeoLines(rsninfo, [0, 0, 0], 2), fnc_setfig;
    subplot(2,3,3), histogram(PCSavg(maskZ==1),[0:0.25:5],'normalization','probability','facealpha',0.6,...
        'edgecolor','k','facecolor','w','linewidth',2), title(['median = ' num2str(median(PCSavg(maskZ==1)))]);
    fnc_setfig, ylim([0 0.18]), xlim([0 4.1]);
    
    disp(['  median PCS (z>1.96) = ' num2str(median(PCSavg(maskZ==1)))]);
    disp(['  5-percentile PCS (z>1.96) = ' num2str(prctile(PCSavg(maskZ==1),5))]);
    disp(['  95-percentile PCS (z>1.96) = ' num2str(prctile(PCSavg(maskZ==1),95))]);
    disp(['  percentage z>1.96 = ' num2str(sum(maskZ(:) * 100 / (nn*(nn-1))))]);
    disp(' ');   
    
    % Plot
    temp = PCSavg(rsninfo.rsnOrder,rsninfo.rsnOrder) .* maskFDR(rsninfo.rsnOrder,rsninfo.rsnOrder);
    ax(3) = subplot(2,3,5); imagesc(temp), axis equal tight, axis off, title('PCS FDR-corrected p<.05');
    colormap(ax(3),othercolor('YlOrBr9')), colorbar, caxis([perc5_human perc95_human]), fnc_addYeoLines(rsninfo, [0, 0, 0], 2), fnc_setfig;
    ax(4) = subplot(2,3,6); histogram(PCSavg(maskFDR==1),[0:0.25:5],'normalization','probability','facealpha',0.6,...
        'edgecolor','k','facecolor','w','linewidth',2), title(['median = ' num2str(median(PCSavg(maskFDR==1)))]);
    fnc_setfig, ylim([0 0.18]), xlim([0 5.1]);

    disp(['  median PCS (FDR) = ' num2str(median(PCSavg(maskFDR==1)))]);
    disp(['  5-percentile PCS (FDR) = ' num2str(prctile(PCSavg(maskFDR==1),5))]);
    disp(['  95-percentile PCS (FDR) = ' num2str(prctile(PCSavg(maskFDR==1),95))]);
    disp(['  percentage FDR = ' num2str(sum(maskFDR(:) * 100 / (nn*(nn-1))))]);
    disp(' ');

end

% Two-sample Kolmogorov-Smirnov test between PCS histograms
hval = nanmean(hPCS,3); temp = triu(hmaskFDR,1);
hval = hval(temp == 1);
qval = nanmean(qPCS,3); temp = triu(qmaskFDR,1);
qval = qval(temp == 1);
mval = nanmean(mPCS,3); temp = triu(mmaskFDR,1);
mval = mval(temp == 1);
[h,p,ks2stat] = kstest2(hval,qval)
[h,p,ks2stat] = kstest2(hval,mval)
[h,p,ks2stat] = kstest2(qval,mval)



%% FIGURES 3a, 3b
clc
close all
flagNullModel = 0;  % 0 : original data; 1 : z<1.96; 2 : FDR-correction 
for s = 1:length(spec)

    disp(['Species ' num2str(s) ' of 3 - ' spec_string{s}]);
    
    % PCS [nn X nn X ns]
    eval(['PCS = ', spec{s}, 'PCS;']);
    eval(['thisInfo = ', spec{s}, 'info;']);
    nn = size(PCS,1);
    
    PCSavg = nanmean(PCS,3);
    if flagNullModel == 1
        eval(['mask = ', spec{s}, 'maskZ;']);
        PCSavg = PCSavg .* mask;
    elseif flagNullModel == 2
        eval(['mask = ', spec{s}, 'maskFDR;']);
        PCSavg = PCSavg .* mask;
    end
    PCSnode = sum(PCSavg) ./ (nn-1);
    minval = prctile(PCSnode,5);   
    maxval = prctile(PCSnode,95);
    
    % Cortical plot (human and macaque data)
    % ---------------------------------------------------------------------
    if spec{s} == 'h'
        path_surf_rh = fullfile(plotkit_dir, 'data', 'fsaverage6', 'surf', 'rh.pial'); % inflated / pial
        path_surf_lh = fullfile(plotkit_dir, 'data', 'fsaverage6', 'surf', 'lh.pial'); % inflated / pial
        path_annot_rh = fullfile(plotkit_dir, 'data', 'fsaverage6', 'label', 'rh.Schaefer2018_100Parcels_7Networks_order.annot'); % SET PARCELLATION
        path_annot_lh = fullfile(plotkit_dir, 'data', 'fsaverage6', 'label', 'lh.Schaefer2018_100Parcels_7Networks_order.annot'); % SET PARCELLATION  
        this_cm = othercolor('YlOrBr9');
        CM = squeeze(mapsc2rgb(PCSnode, this_cm, minval, maxval));
        [k1,k2] = colorsurf_2hemi_5perspectives_Schaefer(path_surf_rh, path_annot_rh, ...
                path_surf_lh, path_annot_lh, CM, thisInfo.Label, 5, 'gouraud');
    end
    if spec{s} == 'q'
        this_cm = othercolor('YlOrBr9');
        CM = squeeze(mapsc2rgb(PCSnode, this_cm, minval, maxval));
        [k3,k4,k5,k6,k7] = colorsurf_2hemi_5perspectives_macaqueRM82(thisInfo, CM, thisInfo.Id, 'gouraud', 5); % nview=5 ~ very slow !!
    end
    
    
    % Average per RSN
    % ---------------------------------------------------------------------
    nrsn = length(thisInfo.rsn);
    PCSnode_rsn = zeros(nrsn,1);
    for i = 1:nrsn
        PCSnode_rsn(i) = mean(PCSnode(thisInfo.rsnROIs == i));
    end
    
    [~,bb] = sort(PCSnode_rsn,'descend');
    this_cm = othercolor('YlOrBr9');
    CM = squeeze(mapsc2rgb(PCSnode_rsn, this_cm, minval, maxval));
    figure, hold on;
    for i = 1:nrsn
        bar(i,PCSnode_rsn(bb(i)),'facecolor','w','edgecolor','k','linewidth',2,'facecolor',CM(bb(i),:),'barwidth',1);
    end
    title([spec_string{s}, ' PCS']);
    xticks(1:nrsn), xticklabels(thisInfo.rsn(bb)), xtickangle(45);
    ylim([min(PCSnode_rsn) * 0.9, max(PCSnode_rsn) * 1.1]), fnc_setfig;
    
    
    % Connections within and between RSNs
    % ---------------------------------------------------------------------
    PCSrsn = nan(nrsn);
    for i = 1:nrsn
        for j = 1:nrsn
            ii = find(thisInfo.rsnROIs==i);
            jj = find(thisInfo.rsnROIs==j);
            temp = PCSavg(ii,jj);
            if i==j
                temp = temp(find(triu(ones(size(temp)),1)));
            end
            PCSrsn(i,j) = mean(temp(:));
        end
    end
    figure, imagesc(PCSrsn), colormap(othercolor('YlOrBr9')), colorbar;
    xticks(1:nrsn), xticklabels(thisInfo.rsn), xtickangle(45), title([spec_string{s}, ' PCS']);;
    yticks(1:nrsn), yticklabels(thisInfo.rsn), axis equal tight, fnc_setfig;
    
end



%% FIGURE 3c
clc
for s = 1:length(spec)

    disp(['Species ' num2str(s) ' of ' num2str(length(spec))]);
    
    % PCS [nn X nn X ns]
    eval(['PCS = ', spec{s}, 'PCS;']);
    eval(['thisInfo = ', spec{s}, 'info;']);
    ns = size(PCS,3);
    nrsn = length(thisInfo.rsn);

    % Loop over subjects
    % Individual average PCS:
    %   1. within and between unimodal/multimodal systems
    %   2. within and between transmodal systems
    %   3. between unimodal/multimodal and transmodal systems
    PCSuti = nan(ns,3);
    for i = 1:ns
        val1 = []; % unimodal/multimodal
        val2 = []; % transmodal
        val3 = []; % inter
        rsn = zeros(nrsn);
        PCSsubj = PCS(:,:,i);
        % Loop over systems (RSNs)
        % Intra-system PCS
        for a = 1:nrsn
            ii = find(thisInfo.rsnROIs == a);
            val = PCSsubj(ii,ii);
            val = val(triu(ones(size(val,1)),1) == 1);
            if thisInfo.unimodal(a)
                val1 = [val1; val(:)];
            else
                val2 = [val2; val(:)];
            end
            if a == nrsn
                continue
            end
            % Inter-system PCS
            for b = a+1:nrsn
                jj = find(thisInfo.rsnROIs == b);
                val = PCSsubj(ii,jj);
                if thisInfo.unimodal(a)
                    if thisInfo.unimodal(b)
                        val1 = [val1; val(:)];
                    else
                        val3 = [val3; val(:)];
                    end
                else
                    if thisInfo.unimodal(b)
                        val3 = [val3; val(:)];
                    else
                        val2 = [val2; val(:)];
                    end
                end
            end
        end
        
        PCSuti(i,1) = mean(val1);
        PCSuti(i,2) = mean(val2);
        PCSuti(i,3) = mean(val3);
    end
    
    % Boxplots
    figure;
    thiscol = [217, 136, 128; 169, 50, 38; 123, 36, 28] ./ 255;
    al_goodplot(PCSuti(:,1),1,0.5,thiscol(1,:));
    al_goodplot(PCSuti(:,2),2,0.5,thiscol(2,:));
    al_goodplot(PCSuti(:,3),3,0.5,thiscol(3,:))
    title([spec_string{s} ' - individual PCS'])
    xticks([1 2 3]), xticklabels({'uni','trans','inter'}), ylabel('individual PCS')
    ylim([0.5 2.5]), fnc_setfig;

    % Statistical tests - within species
	spec_string{s}
    disp(['median unimodal = ' num2str(median(PCSuti(:,1)))]);
    disp(['median transmodal = ' num2str(median(PCSuti(:,2)))]);
    disp(['median inter = ' num2str(median(PCSuti(:,3)))]);
    p = kruskalwallis(PCSuti,[],'off');
    disp(['p-value Kruskas-Wallis test = ' num2str(p)]);
     
    % Store values
    eval([spec{s} 'PCSuti = PCSuti;']);
    
end



%% FIGURE S2
% Between-species statistics

% Select unimodal (ii=1) / transmodal (ii=2) / cross-modal (ii=3)
ii = 1;
val = nan(100,3);
val(1:size(hPCSuti,1),1) = hPCSuti(:,ii);
val(1:size(qPCSuti,1),2) = qPCSuti(:,ii);
val(1:size(mPCSuti,1),3) = mPCSuti(:,ii);
p = kruskalwallis(val,[],'off');
% Boxplots
figure;
thiscol = [hcol; qcol; mcol];
al_goodplot(val(1:size(hPCSuti,1),1),1,0.5,thiscol(1,:));
al_goodplot(val(1:size(qPCSuti,1),2),2,0.5,thiscol(2,:));
al_goodplot(val(1:size(mPCSuti,1),3),3,0.5,thiscol(3,:));
xticks([1 2 3]), xticklabels({'humans','macaques','mice'}), ylabel('individual PCSs');
title(['Unimodal region pairs (p = ' num2str(p) ')']);
ylim([0.5 2.5]), fnc_setfig;

ii = 2;
val = nan(100,3);
val(1:size(hPCSuti,1),1) = hPCSuti(:,ii);
val(1:size(qPCSuti,1),2) = qPCSuti(:,ii);
val(1:size(mPCSuti,1),3) = mPCSuti(:,ii);
p = kruskalwallis(val,[],'off');
% Boxplots
figure;
thiscol = [hcol; qcol; mcol];
al_goodplot(val(1:size(hPCSuti,1),1),1,0.5,thiscol(1,:));
al_goodplot(val(1:size(qPCSuti,1),2),2,0.5,thiscol(2,:));
al_goodplot(val(1:size(mPCSuti,1),3),3,0.5,thiscol(3,:));
xticks([1 2 3]), xticklabels({'humans','macaques','mice'}), ylabel('individual PCSs');
title(['Transmodal region pairs (p = ' num2str(p) ')']);
ylim([0.5 2.5]), fnc_setfig;

ii = 3;
val = nan(100,3);
val(1:size(hPCSuti,1),1) = hPCSuti(:,ii);
val(1:size(qPCSuti,1),2) = qPCSuti(:,ii);
val(1:size(mPCSuti,1),3) = mPCSuti(:,ii);
p = kruskalwallis(val,[],'off');
% Boxplots
figure;
thiscol = [hcol; qcol; mcol];
al_goodplot(val(1:size(hPCSuti,1),1),1,0.5,thiscol(1,:));
al_goodplot(val(1:size(qPCSuti,1),2),2,0.5,thiscol(2,:));
al_goodplot(val(1:size(mPCSuti,1),3),3,0.5,thiscol(3,:));
xticks([1 2 3]), xticklabels({'humans','macaques','mice'}), ylabel('individual PCSs');
title(['Corss-modal region pairs (p = ' num2str(p) ')']);
ylim([0.5 2.5]), fnc_setfig; 


