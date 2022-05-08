function [sr, IM, Idiff, Iothers, Iself] = fnc_identifiability(X, subj_label, titlestr, cmap, doplot, valmin, valmax)

%
% Compute identifiability matrix and success rate for input data X,
% subj_label.
%
% NOTE: the script assume a data structure with test-retest data. The X
% matrix is expected to contain first the test data, and then the retest
% data, in the same subject order!
%
% INPUTs:
% X : feature matrix (n_features X total_n_datapoints), with
%     total_n_datapoints = 2 * n_subjects
% subj_label : numberical subject IDs (1 X total_n_datapoints)
% titlestr : string for plots' title
% cmap : colormap
% doplot : 1/0 - generate plots yes/no
%
% OUTPUTs:
% sr : success rate (see Finn et al., 2015)
% IM : identifiability matrix (see Amico et al., 2018)
% Iself : average of IM diagonal (see Amico et al., 2018)
% Iothers : average of IM outside-diagonal (see Amico et al., 2018)
% Idiff : Iself - Iothers (see Amico et al., 2018)
%
% Alessandra Griffa
% University of Geneva
% May 2022
%

% Number of subjects
ns = length(unique(subj_label));

% SIMILARITY MEASURE
CM = squareform(pdist(X'+1,'jaccard')); 
CM = 1 - CM;    % from distance to similarity
CM = CM .* (~eye(size(CM,1)));  % remove correlation of an acquisition with itself
CM(1:ns,1:ns) = nan;

% Identifiability matrix (IM)
ii_test = 1:ns;
ii_retest = ns+1 : 2*ns;
IM = CM(ii_test, ii_retest);

% Compute Finn Success Rate (sr)
sr = 0;
for i = 1:ns
    true_label = subj_label(i);
    test_label = subj_label(find(IM(i,:) == max(IM(i,:))));
    if true_label == test_label
        sr = sr + 1;
    end
end
sr = sr / ns;

% Iself, Iothers, Idiff
ii_other = find(~eye(length(ii_test))); % indexes entries outside diagonal
ii_self = find(eye(length(ii_test)));   % indexes entries on diagonal
Iself = nanmean(IM(ii_self));
Iothers = nanmean(IM(ii_other));
Idiff = Iself - Iothers;

% Plotting
if doplot
    
    temp = othercolor(cmap);
    col1 = temp(256,:);
    col2 = temp(180,:);
    figure;
    al_goodplot(IM(ii_self),1,0.5,col1);
    al_goodplot(IM(ii_other),2,0.5,col2);
    grid off, xticks([1 2]), xticklabels({'self','other'}), ylabel('Jaccard index');
    fnc_setfig;
    if ~isnan(valmin) && ~isnan(valmax)
        ylim([valmin, valmax]);
    end
    
    temp = IM(find(eye(size(IM))==1));
    [a,b] = sort(temp,'descend');
    figure, imagesc(IM(b,b)), colormap(othercolor(cmap)), colorbar, axis equal tight, ylabel('subjects - test'), xlabel('subjects - retest');
    title([{['SR = ' num2str(sr,3)]}]), fnc_setfig;
    if ~isnan(valmin) && ~isnan(valmax)
        caxis([valmin, valmax]);
    end

end
    
