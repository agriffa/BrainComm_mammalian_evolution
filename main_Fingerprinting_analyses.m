%
% This code reproduce the fingerprinting analisys of the article 'The 
% evolution of information transmission in mammalian brain networks' -
% Alessandra Griffa, Mathieu Mach, Julien Dedelley, Daniel Gutierrez-Barragan, 
% Gutierrez-Barragan, Alessandro Gozzi, Gilles Allali, Joanes Grandjean, 
% Dimitri Van De Ville, Enrico Amico.
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
data_dir = '/Users/alli/Desktop/PROJECT_IT/DATASET_BrainComm_mammalian_evolution/DATA';



%% Preliminary settings
% Add directories to Matlab path
gitdir = pwd;
addpath(genpath(fullfile(gitdir,'utility')));

% Structural connectivity weigth 
wei = 'mm';              % Euclidean distance (in millimiters)

% NOTE: the first letter of 'spec_string' and 'spec' entries indicates the
% species (h = human; q = macaque; m = mouse)
spec_string = {'hHCP','qTVB','mAD3'};
spec = {'h','q','m'};  

% Set test-retest datasets dirrectories
htest_dir = fullfile(data_dir, [spec_string{1},'_test']);
hretest_dir = fullfile(data_dir, [spec_string{1},'_retest']);
qtest_dir = fullfile(data_dir, [spec_string{2},'_test']);
qretest_dir = fullfile(data_dir, [spec_string{2},'_retest']);
mtest_dir = fullfile(data_dir, [spec_string{3},'_test']);
mretest_dir = fullfile(data_dir, [spec_string{3},'_retest']);

% Colors
hcol = [25, 118, 210]./255;     % human: blue
qcol = [0, 137, 123]./255;      % macaque: green
mcol = [216, 67, 21]./255;      % mouse: orange
hmap = 'Blues9';
qmap = 'Greens9';
mmap = 'Reds9';



%% LOAD PCS test-retest data, and PCS whole-group data
measure = 'PCS'; 
% Loop over species
for s = 1:length(spec)
    
    disp('');
    disp('');
    disp(['    ' spec_string{s} ' - dataset ' num2str(s) ' of ' num2str(length(spec_string)) ' - LOAD TEST-RETEST DATA']);
    
    % LOAD TEST-RETEST DATA
    eval(['test_dir = ', spec{s}, 'test_dir;']);
    eval(['retest_dir = ', spec{s}, 'retest_dir;']);
    
    test_filename = dir(fullfile(test_dir,[measure,'*.mat']));
    retest_filename = dir(fullfile(retest_dir,[measure,'*']));

    if length(test_filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - PCS TEST FILE']);
        continue
    else
        a = load(fullfile(test_dir,test_filename.name));
        eval(['test_data = a.', measure ,';']);
    end
    
    if length(retest_filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - PCS RETEST FILE']);
        continue
    else
        a = load(fullfile(retest_dir,retest_filename.name));
        eval(['retest_data = a.', measure ,';']);
    end

    % BUILD FEATURE MATRIX FOR FINGERPRINTING
    nn = size(test_data,1);
    ns = size(test_data,3);
    subj_label = [1:ns, 1:ns];
    X = nan(nn * (nn-1) / 2, length(subj_label));
    ii = find(triu(ones(nn),1) == 1);
    for i = 1:ns
        temp = test_data(:,:,i);
        X(:,i) = temp(ii);
        temp = retest_data(:,:,i);
        X(:,ns+i) = temp(ii);
    end
    eval([spec{s},'subj_label = subj_label;']);
    eval([spec{s},'X = X;']);
    figure, subplot(1,2,1), imagesc(X), xlabel('subjects'), ylabel('brain region pairs'), title([measure ' - ' spec_string{s}]), 
    colormap(jet), colorbar, fnc_setfig;
    subplot(1,2,2), imagesc(subj_label'), colorbar, title('subject ID'), fnc_setfig;
    
    % LOAD WHOLE-GROUP DATA
    this_data_dir = fullfile(data_dir, spec_string{s});
    filename = dir(fullfile(this_data_dir,['PCS_*',wei,'*.mat']));
    if length(filename) ~= 1
        warning(['Missing data for ', spec_string{s}, ' - PCS FILE']);
    else
        a = load(fullfile(this_data_dir,filename.name));
        eval([spec{s},'PCS = a.PCS;']);
    end
    
    clear a
    
end



%% Success rate and identifiability matrices
nrand = 500;
close all
disp('PRESS A KEY TO CONTINUE TO NEXT DATASET');
% Loop over species
for s = 1:length(spec)
    
    clc
    
    disp('');
    disp('');
    disp(['    ' spec_string{s} ' - species ' num2str(s) ' of ' num2str(length(spec))]);
    
    eval(['X = ', spec{s}, 'X;']);
    eval(['subj_label = ', spec{s}, 'subj_label;']);
    eval(['cmap = ', spec{s}, 'map;']);
    ns = length(unique(subj_label));
    
    [sr, IM, Idiff, Iothers, Iself] = fnc_identifiability(X, subj_label, spec{s}, cmap, 1, nan, nan);
    
    % Randomized data (random assignment test-retest; nrand = 100)
    sr_rand = zeros(nrand,1);
    for r = 1:nrand
        % Scrumble retest data
        IMrand = IM(:, randperm(ns));
        % Compute Finn Success Rate
        for i = 1:ns
            true_label = subj_label(i);
            test_label = subj_label(find(IMrand(i,:) == max(IMrand(i,:))));
            if true_label == test_label
                sr_rand(r) = sr_rand(r) + 1;
            end
        end
        sr_rand(r) = sr_rand(r) / ns;
    end    
    
    disp([spec_string{s} ' - SR = ' num2str(sr)]);
    disp([spec_string{s} ' - SR NULL = ' num2str(mean(sr_rand)) ' (' num2str(std(sr_rand)) ')']);
    disp(' ');
    disp('  ... press a key to continue to next datase');
    
    pause
    close all
    
end



%% Is subject identifiability driven by selective / parallel communication?
% Compute Jaccard for low and high PCS scores
nrand = 500;
close all
disp('PRESS A KEY TO CONTINUE TO NEXT DATASET');
for s = 1:length(spec)
    
    clc
    
    disp('');
    disp('');
    disp(['    ' spec_string{s} ' - species ' num2str(s) ' of ' num2str(length(spec))]);
    
    eval(['PCS = ', spec{s}, 'PCS;']);
    eval(['X = ', spec{s}, 'X;']);
    eval(['subj_label = ', spec{s}, 'subj_label;']);
    eval(['cmap = ', spec{s}, 'map;']);
    ns = length(unique(subj_label));
    nn = size(PCS,1);
    
    % Median split of group-average PCS values
    PCSavg = mean(PCS,3);
    PCSavg = PCSavg(find(triu(ones(nn),1)));
    %medianval = median(PCSavg);
    medianval = 1.3;
    disp([spec_string{s} ' - median = ' num2str(medianval)]);
    
    Xlow = X(PCSavg <= medianval,:);
    Xhigh = X(PCSavg > medianval,:);
    
    [sr_low, IM_low, Idiff_low, Iothers_low, Iself_low] = fnc_identifiability(Xlow, subj_label, [spec{s}, 'LOW'], cmap, 1, nan, nan);
    close(gcf)
    title('low-PCS');
    [sr_high, IM_high, Idiff_high, Iothers_high, Iself_high] = fnc_identifiability(Xhigh, subj_label, [spec{s}, 'HIGH'], cmap, 1, nan, nan);
    close(gcf)
    title('high-PCS');
    
    % Randomized data (random assignment test-retest; nrand = 100)
    sr_rand_low = zeros(nrand,1);
    sr_rand_high = zeros(nrand,1);
    for r = 1:nrand
        ii = randperm(ns);
        % Scrumble retest data
        IMrand_low = IM_low(:,ii);
        IMrand_high = IM_high(:,ii);
        % Compute Finn Success Rate
        for i = 1:ns
            true_label = subj_label(i);
            test_label = subj_label(find(IMrand_low(i,:) == max(IMrand_low(i,:))));
            if true_label == test_label
                sr_rand_low(r) = sr_rand_low(r) + 1;
            end
            test_label = subj_label(find(IMrand_high(i,:) == max(IMrand_high(i,:))));
            if true_label == test_label
                sr_rand_high(r) = sr_rand_high(r) + 1;
            end
        end
        sr_rand_low(r) = sr_rand_low(r) / ns;
        sr_rand_high(r) = sr_rand_high(r) / ns;
    end
    
    disp([spec_string{s} ' - SR low-PCS = ' num2str(sr_low)]);
    disp([spec_string{s} ' - SR low-PCS NULL = ' num2str(mean(sr_rand_low)) ' (' num2str(std(sr_rand_low)) ')']);
    disp([spec_string{s} ' - SR high-PCS = ' num2str(sr_high)]);
    disp([spec_string{s} ' - SR high-PCS NULL = ' num2str(mean(sr_rand_high)) ' (' num2str(std(sr_rand_high)) ')']);    
    disp('  ... press a key to continue to next datase');
        
    pause
    close all
    
end


