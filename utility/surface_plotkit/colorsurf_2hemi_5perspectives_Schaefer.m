%
%
% 1) read surface as set of nodes and vertices 
% 2) read annotation file, which contains annotation for each vertex of the
%    corresponding surface
% 3) take in input colormap 
%
% INPUT: - path_surf_rh (e.g.: '.../rh.pial')
%        - path_annot_rh (e.g.: '.../lh.pial')
%        - path_surf_lh (e.g.: '.../rh.myaparc_36.annot')
%        - path_annot_lh (e.g.: '.../lh.myaparc_36.annot')
%        - CM: user colormap (it is assumed not to contain subcortical ROIs)
%           3D array (row: ROIs; columns: color (1 or 3 column); 3rd dimension: different time points)
%        - llist: list of regions names in CM order 
%        - nview: can be set equal to 4 (the function outpus a single image
%                 with four, coronal brian viewes), or the 5 (the function
%                 outputs two images, the second one being a transversal
%                 view of the brain)
%        - my_lighting: lighting option. E.g.: 'none'; 'gouraud'
% OUTPUT: save ave file
% REQUIRE: - FreeSurfer .m scripts ()
%          - thight_subplot ()
%
% Alessandra Griffa
% Ecole polytechnique federale de Lausanne EPFL
% Signal Processing Laboratory 5 LTS5  | MIPLab
% Jun 2020
%

function [k1,k2] = colorsurf_2hemi_5perspectives(path_surf_rh, path_annot_rh, path_surf_lh, path_annot_lh, CM, llist, nview, my_lighting)

% Initialize figures handles
k1 = NaN;
k2 = NaN;

% LOAD Right hemisphere
[v_rh,f_rh] = read_surf(path_surf_rh); % v: vertices; f: faces of trangles
nv_rh = size(v_rh,1); % number of vertices
if min(f_rh(:)) == 0  %faces in freesurfer start at 0 (in matlab, they start at 1);
    f_rh = f_rh + 1;
end
[~,label_rh,colortable_rh] = read_annotation(path_annot_rh);

% LOAD Left hemisphere
[v_lh,f_lh] = read_surf(path_surf_lh);
nv_lh = size(v_lh,1); % number of vertices
if min(f_lh(:)) == 0  %faces in freesurfer start at 0 (in matlab, they start at 1);
    f_lh = f_lh + 1;
end
[~,label_lh,colortable_lh] = read_annotation(path_annot_lh);



%% Assign a color to each vertex
% Generate a matrix vColor(#verices x 3), where each line contains the
% color of the corresponding vertex.
vColor_rh = ones(nv_rh,size(CM,2)) .* 211/255;
vColor_lh = ones(nv_lh,size(CM,2)) .* 211/255;
for i = 1:length(llist)
    index_rh = find(strcmp(llist{i},colortable_rh.struct_names));
    if ~isempty(index_rh)
        ii = colortable_rh.table(index_rh,5);
        indexVerticesRoi_rh = find(label_rh == ii);
        thisColor_rh = CM(i,1:end);
        vColor_rh(indexVerticesRoi_rh,1:end) = repmat(thisColor_rh,length(indexVerticesRoi_rh),1);
    end
    index_lh = find(strcmp(llist{i},colortable_lh.struct_names));
    if ~isempty(index_lh)
        ii = colortable_lh.table(index_lh,5);
        indexVerticesRoi_lh = find(label_lh == ii);
        thisColor_lh = CM(i,1:end);
        vColor_lh(indexVerticesRoi_lh,1:end) = repmat(thisColor_lh,length(indexVerticesRoi_lh),1);
    end
end



%% PLOT
% PLOT SINGLE HEMISPHERES
% Surface plot according to vColor colortable
%close all;
k1 = figure;%('Visible','on');

if nview == 4
    ha = tight_subplot(2,2,[0 0.01],[.01 .01],[.01 .01]); 
elseif nview == 5
    ha = tight_subplot(2,2,[0.01 0.2],[.01 .01],[.01 .01]);
else
    ha = tight_subplot(2,2,[0 0.01],[.01 .01],[.01 .01]);
end

% LEFT HEMI SIDE
axes(ha(1)); 
p = patch('faces',f_lh,'vertices',v_lh,'facecolor','flat','edgecolor','none','facealpha',1);
set(p,'FaceVertexCData',vColor_lh);
daspect([1 1 1]);
view([-90 0]);
camlight('headlight');
%lighting gouraud; % specify lighting algorithm        
lighting(my_lighting);
%alpha(0.3)
axis tight off;    

% RIGHT HEMI SIDE
axes(ha(2)); 
p = patch('faces',f_rh,'vertices',v_rh,'facecolor','flat','edgecolor','none','facealpha',1);
set(p,'FaceVertexCData',vColor_rh);
daspect([1 1 1]);
view([90 0]);
camlight('headlight');
lighting(my_lighting);
%alpha(0.3)
axis tight off;

% RIGHT HEMI MEDIAL
axes(ha(4)); 
p = patch('faces',f_rh,'vertices',v_rh,'facecolor','flat','edgecolor','none','facealpha',1);
set(p,'FaceVertexCData',vColor_rh);
daspect([1 1 1]);
view([-90 0]);
camlight('headlight');
lighting(my_lighting); % specify lighting algorithm
%alpha(0.3)
axis tight off;   

% LEFT HEMI MEDIAL
axes(ha(3)); 
p = patch('faces',f_lh,'vertices',v_lh,'facecolor','flat','edgecolor','none','facealpha',1);
set(p,'FaceVertexCData',vColor_lh);
daspect([1 1 1]);
view([90 0]);
camlight('headlight');
lighting(my_lighting);
%alpha(0.3)
axis tight off; 

set(gcf, 'Color', 'None')


% PLOT BOTH HEMISPHERE, TRANSVERSAL VIEW
if nview == 5
    k2 = figure;
    ha = tight_subplot(1,1,[0 0],[.2 .2],[.2 .2]);
    axes(ha(1));
    p = patch('faces',f_rh,'vertices',v_rh,'facecolor','flat','edgecolor','none','facealpha',1);
    set(p,'FaceVertexCData',vColor_rh);
    p = patch('faces',f_lh,'vertices',v_lh,'facecolor','flat','edgecolor','none','facealpha',1);
    set(p,'FaceVertexCData',vColor_lh);
    daspect([1 1 1]);
    view(3); % view(3) sets the default three-dimensional view, az = ???37.5, el = 30
    view([50 -40 70]); % sets the viewpoint to the Cartesian coordinates x, y, and z
    camlight; % creates a light right and up from camera
    lighting(my_lighting);
    set(gcf,'color','w');
    view([0 90]);
    axis off;
    set(gcf, 'Color', 'None');
end
















