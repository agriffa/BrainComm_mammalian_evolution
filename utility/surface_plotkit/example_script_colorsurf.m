
% Example script
% I assume you have loaded structural connectivity brain data, including the 'label' variable


% Add path
plotkit_path = '/Users/alli/Desktop/PROJECT_IT/CODE/Alessandra/surface_plotkit';
addpath(genpath(plotkit_path));


% Set files for cortical suface plot
path_surf_rh = fullfile(plotkit_path, 'data', 'fsaverage6', 'surf', 'rh.pial');
path_surf_lh = fullfile(plotkit_path, 'data', 'fsaverage6', 'surf', 'lh.pial');
path_annot_rh = fullfile(plotkit_path, 'data', 'fsaverage6', 'label', 'rh.Schaefer2018_200Parcels_7Networks_order.annot');
path_annot_lh = fullfile(plotkit_path, 'data', 'fsaverage6', 'label', 'lh.Schaefer2018_200Parcels_7Networks_order.annot');



%%
% Plot CORTICAL REGIONS ONLY (regions 1 to 200)
% EXAMPLE 1: random colormap 
% Call plotting function ~ colorsurf_2hemi_5perspectives_Schaefer
% INPUTS
% path_surf_rh: this is a cortical surface file for the right hemisphere. Options: rh.inflated or rh.pial
% path_annot_rh: this is a file assigning to each vertex of the right cortical surface an anatomical label (e.g., 7Networks_RH_Default_pCunPCC_1, etc.)
% path_surf_lh: this is a cortical surface file for the right hemisphere. Options: lh.inflated or lh.pial
% path_annot_rh: this is a file assigning to each vertex of the left cortical surface an anatomical label (e.g., 7Networks_LH_Default_pCunPCC_1, etc.)
% CM: colormap matrix
% llist: cortical region labels (~ labels(1:200))
% n_view: 4 or 5 brain view (it makes sense to select 5 views only if you sue the xx.pial surface)
% my_lighting: Matlab lighting string - E.g.: 'gouraud'; 'none'
n_corticalnodes = 200;
CM = rand(n_corticalnodes, 3); % ColorMap matrix: one RGB color per node
[k1, k2] = colorsurf_2hemi_5perspectives_Schaefer(path_surf_rh, path_annot_rh, ...
           path_surf_lh, path_annot_lh, CM, labels(1:200), 5, 'gouraud');

% SAVE FIGURE AS HIGH RESOLUTION TIF IMAGE
export_fig(k1, '/Users/alli/Desktop/test_view1.tif', '-transparent', '-m2');



%%
% Plot CORTICAL REGIONS ONLY (regions 1 to 200)
% EXAMPLE 1: user-defined nodal colors [what you need :)]

% Create nodal example data
test = sum(SC.mll(:,:,1));
test = R;
figure, plot(test,'*'), xlabel('node numerical id');

% Assign colors to nodes according to test values and a given colormap
% For the colormap you can use the classic Matlab ones (e.g., this_cm =
% colormap(jet);), or the colormaps degined in the 'othercolor' toolbox
% (https://www.mathworks.com/matlabcentral/fileexchange/30564-othercolor)
this_cm = othercolor('RdBu11'); close all 
% The mapsc2rgb function maps the values in 'test' to RGB colors according to 'this_cm' colormap
max_val = max(test);
min_val = min(test);
CM = squeeze(mapsc2rgb(test, this_cm, min_val, max_val));

% Plot
[k1, k2] = colorsurf_2hemi_5perspectives_Schaefer(path_surf_rh, path_annot_rh, ...
           path_surf_lh, path_annot_lh, CM, labels(1:200), 5, 'gouraud');
