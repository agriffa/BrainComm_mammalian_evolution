function [] = fnc_addYeoLines(yeoinfo, col, lw)

% yeoinfo: yeoinfo.rsnOrder
%          yeoinfo.rsnROIs
% col: line color
% lw: line width


hold on;

% Set line color
%col = [213, 216, 220] ./ 255;

% Nodal RSNs assignment, assuming that nodes are ordered according to
% rsnOrder
a = yeoinfo.rsnROIs(yeoinfo.rsnOrder);


% col_rsn = [130,57,139;...
%             101,154,196;...
%             26,133,53;...
%             159,103,164;...
%             232,237,181;...
%             246,173,70;...
%             216,92,110;...
%             55,85,158;...
%             41,41,41] / 255;
% figure, imagesc(a), colormap(col_rsn), axis off, set(gcf,'color','w')


% Loop over RSNs (do not consider the cerebellum, which has only 2 regions)
nrsn = max(a);
for i = 1:nrsn
    ii = find(a==i);
    % Left line
    plot([ii(1)-0.5 ii(1)-0.5],[ii(1) ii(end)], 'color', col, 'linewidth', lw);
    % Right line
    plot([ii(end)+0.5 ii(end)+0.5],[ii(1) ii(end)], 'color', col, 'linewidth', lw);
    % Top line
    plot([ii(1) ii(end)],[ii(1)-0.5 ii(1)-0.5], 'color', col, 'linewidth', lw);    
    % Bottom line
    plot([ii(1) ii(end)],[ii(end)+0.5 ii(end)+0.5], 'color', col, 'linewidth', lw);
end