function [hFig, hAxes, return_chan_loc] = be_viz_topo(channels, data, varargin)

global GlobalData

colormap(jet)
colorbar
%caxis([50 100]) % change colorbar limits POUR Z SCORE
%caxis([0 0.04]) %meandiffall gamma1-gammaL
%caxis([0.06 0.12]) %meandiffall BETA
%caxis([0 0.45]) %meandiffall basse freq
%caxis([-2 2]) % valeur t après le paired t-test
%caxis([1e-29 3e-29]) % good for one
%caxis([1e-30 1.5e-28])
caxis([0 0.6])


% ===== Procss optional inputs
nD  =   2;
if numel(varargin)>0
    nD  =   varargin{1};
end
Mask    =   [];
if numel(varargin)>1
    Mask    =   varargin{2};
end
hAxes       =   [];
if numel(varargin)>2
    hAxes   =   varargin{3};
end

% ===== SET UP NEW FIGURE =====
FigureId            =   db_template('FigureId');
FigureId.Type       =   '3DViz';
FigureId.SubType    =   '';
FigureId.Modality   =   [];

% Create figure
hFig                =   [];
if isempty(hAxes)
    hFig                = figure_3d('CreateFigure', FigureId);
    hAxes               =   findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hChildren           = get(hAxes, 'Children');
    delete(hChildren(~strcmpi(get(hChildren, 'Type'), 'light')));
    % Set Topography axes as current axes
    set(0,    'CurrentFigure', hFig);
    set(hFig, 'CurrentAxes',   hAxes);
    hold on
    % Reset initial zoom factor and camera position
    figure_3d('ResetView', hFig);
end
    
% ===== MAKE TOPOGRAPHY =====
% commons
chan_loc        =   ( cellfun( @(a) mean(a(:,1:4),2)', {channels.Loc}, 'uni', 0) )';
chan_loc        =   cell2mat( chan_loc );

% Compute best fitting sphere from sensors
[bfs_center, bfs_radius] = bst_bfs(chan_loc);


% ===== CREATE A HIGH-DEF SURFACE =====
% Remove the duplicated positions
precision = 1e5;
Vertices = unique(round(chan_loc * precision)/precision,'rows');

% Tesselate sensor cap
Faces = channel_tesselate(Vertices, 1);
% Clean from some very pathological triangles
Faces = tess_threshold(Vertices, Faces, 3, []);

% Refine mesh
[Vertices, Faces] = tess_refine(Vertices, Faces, [], 1);
if (length(Vertices) < 800)
    [Vertices, Faces] = tess_refine(Vertices, Faces, [], 1);
end
    
% ===== MAKE INTERPOLATION =====
% Get extrapolation options
epsilon                 =   1e-4;
Whitener                =   [];

% Compute normals at each vertex
VerticesNormals         =   tess_normals(Vertices, Faces);

% Create a pseudo-channel file of virtual magnetometers
destChan                =   repmat(db_template('channeldesc'), [1, length(Vertices)]);
for i = 1:length(Vertices)
    destChan(i).Name    =   num2str(i);
    destChan(i).Type    =   'MEG';
    destChan(i).Loc     =   Vertices(i,:)';
    destChan(i).Orient  =   VerticesNormals(i,:)';
    destChan(i).Weight  =   1;
end

% Compute interpolation
% [WExtrap, src_xyz]      = 	be_computeinterp(channels, destChan, bfs_center, bfs_radius, Whitener, epsilon);

% % Apply interpolation matrix sensors => display surface
% DataToPlot              =   WExtrap * data;
% DataToPlot = griddata(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3), data, Vertices(:,1), Vertices(:,2), Vertices(:,3), 'linear'); 
F = scatteredInterpolant(chan_loc, data,'linear', 'linear'); 
DataToPlot = F(Vertices);
% DataToPlot = griddatan(chan_loc, data, Vertices, 'nearest'); 
% if any(DataToPlot)
%     range1                  =   max(data)-min(data);
%     range2                  =   max(DataToPlot)-min(DataToPlot);
%     DataToPlot              =   (DataToPlot - min(DataToPlot)) / range2 * range1 + min(data);
% end

%% === STAT MASK
if ~isempty(Mask)
 
    % % Mask topo
    % Mtp     =   WExtrap * Mask; 
    % Mtp     =   (Mtp - min(Mtp)) / (max(Mtp) - min(Mtp));
    % Srt     =   sort( Mtp, 'descend' );
    % thresh  =   Srt( round(sum(Mask)/numel(Mask)*numel(Mtp)) );
    
    % % Unplottable faces    
    % iMask   =   Mtp>thresh;
    Mask = boolean(Mask);
    chan_loc = chan_loc(Mask,:);
end

if nD==2
    [X,Y] = bst_project_2d(Vertices(:,1), Vertices(:,2), Vertices(:,3));
    % Get 2D vertices coordinates, re-tesselate
    Vertices    = [X, Y, 0*X];
    Faces       = delaunay(X,Y);
    % Clean from some pathological triangles
    %Faces      = tess_threshold(Vertices, Faces, 20, 179.6);
    Faces       = tess_threshold(Vertices, Faces, 20, []);
    % Plot nose / ears
    radii = [Vertices(:,2);Vertices(:,1)];
    figure_topo('PlotNoseEars', hAxes, (max(radii)-min(radii))/4, 1);   
    [chan_loc(:,1),chan_loc(:,2)] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3));
    chan_loc(:,3) = zeros(size(chan_loc(:,1)));
end
return_chan_loc = chan_loc;
% ===== DRAW INTERPOLATION =====
% Get figure colormap
% ColormapInfo            =   getappdata(hFig, 'Colormap');
% sColormap               =   bst_colormaps('GetColormap', ColormapInfo.Type);
    
% ===== Map data on target patch =====
hSurf1 = patch('Faces',     Faces, ...
    'Vertices',  Vertices, ...
    'BackfaceLighting', 'lit', ...
    'Parent',    hAxes, ...
    'FaceVertexCData', DataToPlot, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp');
           
% Surface lighting
material ([.95 0 0])
lighting gouraud
plot3(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3), 'rp', 'MarkerSize', 20)
if exist('hFig', 'var')
    set(hFig, 'Visible', 'on')
end

return