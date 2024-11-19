%SUMMARY
% Author: Connor Gallimore
% 12/26/2023

% This function creates arcs of a circle that represent percentages as the
% proportion of total circumference. I first saw this data visualization
% idea in a Nature paper (below) and wanted to implement it here in MATLAB.

% (Fig 4, Soula et al. 2023, https://www.nature.com/articles/s41593-023-01270-2)

% Required arguments:
    % 'data', a vector or matrix of proportions/percentages. 

    % 'dim',  the dimension containing each sub-cat/ component of the total
    %         e.g. if data series are distributed along the rows, such that
    %              sub-categories / components are arranged along the 
    %              columns, your second argument would be '2'
    %              if data series are distributed along columns and
    %              sub-cats along the rows, your second arg would be '1'

% Optional Name,Value pairs:
    % 'color', an m x 3 vector or matrix specifying RGB triplet(s)
    % 'orientation', 'horizontal' or 'vertical', for more than 1 data 
    %                series, determines whether they are plotted from left-
    %                to-right or top-to-bottom
    % 'precision', specifies the rounding precision for text labels (i.e. 
    %              the max number of decimal places)
    % 'innerRadius', specifies the radius of the circle, indirectly impacting
    %           the width of the line arcs
    % 'scheme', 'category' or 'series', determines the coloring scheme. In
    %           one case, your 'color' matrix may represent the color you 
    %           want each common 'category' to be for all series (default).
    %           In another case, you may be specifying the base color you
    %           want your percentage components to be for each 'series'. In
    %           this option, subsequent percentages will be plotted darker


% Outputs:
% a structure of handles 'H' to the line 'arcs', text labels ('lbls'), and
% 'colors'

% Examples: 
% see 'circPercent_demo.mlx' for usage tips and tricks

% Additional notes:
% The origin and radius are arbitrary, as the width of line arcs will be 
% scaled as a function of radius and number of plot series, with axes 
% removed at the end anyways. I chose (0, 0) and 3, respectively.

% I used the helper function 'getAlignmentFromAngle' from MATLAB's 'pie.m'
% to position text labels, for which I take no credit and all rights go to
% Clay M. Thompson 3-3-94 and The MathWorks, Inc., Copyright 1984-2022 

% This function also uses 'distinguishable_colors.m' by Timothy E. Holy, 
% for which all rights are reserved to them (Copyright (c) 2010-2011, 
% https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
%--------------------------------------------------------------------------


function H = circPercent(data, dim, varargin)

% Initialize input parser here / and or modularize inside a function

% defaults -- center origin (h, k) and radius (r)
inner_r= 0.6;
h= 0;
k= 0;
R= 10;
r= R * inner_r; 

% anonymous functions validating proportions & optional color matrix input
validProps= @(x) isnumeric(x) && ismatrix(x) && all(all(x >= 0) & all(x <= 1));
validColor= @(x) (validProps(x) && size(x, 2) == 3) || all(x >= 1) && size(x, 2) == 1; 

% parse Name,Val pairs
tmp_fco= strcmpi(varargin, 'faceColor'); 
tmp_eco= strcmpi(varargin, 'edgeColor'); 
tmp_txc= strcmpi(varargin, 'textColor'); 
tmp_elw= strcmpi(varargin, 'lineWidth');
tmp_ori= strcmpi(varargin, 'orientation');
tmp_pre= strcmpi(varargin, 'precision'); 
tmp_rad= strcmpi(varargin, 'innerRadius'); 
tmp_sch= strcmpi(varargin, 'scheme'); 
tmp_ang= strcmpi(varargin, 'startAngle');
tmp_prs= strcmpi(varargin, 'patchRes'); 


if any(tmp_ori);  ori= varargin{find(tmp_ori) + 1};
else;             ori= 'horizontal';   % default horizontal
end

if any(tmp_sch);  sc= varargin{find(tmp_sch) + 1};
else;             sc= 'category';      % default color by category
end

if any(tmp_pre);  prec= varargin{find(tmp_pre) + 1};
else;             prec= 2;             % default round to second decimal
end

if any(tmp_rad);  usr_r= varargin{find(tmp_rad) + 1};
    if ~isempty(usr_r) && validProps(usr_r);   r= R * usr_r;
    else;                                      r= R * inner_r; 
    end                                % default 0.6 unless valid usr input
end                          

if any(tmp_ang);  a= varargin{find(tmp_ang) + 1};
else;             a= 0;                % default start ang is 0 (3 o'clock)
end

if any(tmp_prs);  patch_res= varargin{find(tmp_prs) + 1}; 
else;             patch_res= 300;      % default npts for largest arc
end

if any(tmp_eco);  ec= varargin{find(tmp_eco) + 1}; 
else;             ec= 'k'; 
end

if any(tmp_txc);  tc= varargin{find(tmp_txc) + 1}; 
else;             tc= 'k'; 
end

if any(tmp_elw);  lw= varargin{find(tmp_elw) + 1}; 
else;             lw= 1; 
end


% if groups/series distributed along columns, transpose
if dim ~= 2 
    data= data';
    dim= 2; 
end

n_dims= length(size(data)); 
nc_dim= dim; 
np_dim= find(~ismember(1:n_dims, dim));

np= size(data, np_dim);   % num percents to plot
nc= size(data, nc_dim);   % num components to total

% normalize groups that don't sum to 1
if any(sum(data, dim) > 1+sqrt(eps))
    ix2n= sum(data, dim) > 1+sqrt(eps);  % indices to normalize
    data(ix2n, :)= data(ix2n, :) ./ sum(data(ix2n, :), dim);
end

% pre-allocate cells for arcs and cartesian coordinates
[thetas, arcrds, x, y]= deal(cell(np, nc)); 

% pre-allocate txt strings for later based on rounding precision
txt= strcat(string(round(data .* 100, prec)), repmat("%", np, nc)); 

% parse color / color scheme args
if any(tmp_fco)
    col= varargin{find(tmp_fco) + 1};
    if validColor(col)
        fc= col; 
    elseif isempty(col)
        if strcmpi(sc, 'series')
            fc= distinguishable_colors(np*2);
            fc= fc(np+1:end, :);
        else
            fc= distinguishable_colors(nc*2);
            fc= fc(nc+1:end, :);
        end
    else
        error('make sure your color matrix is either mx3 or left empty []')
    end
else
    fc= distinguishable_colors(nc*2);
    fc= fc(nc+1:end, :);
end
dark= linspace(1, 0, nc+1);  % darkness scaling
dark= dark(1:nc)'; 

% anonymous function for color decision rule
colorCat= @(x, y, z) strcmpi(x, 'category') && size(y, 1) == z; 

% resolve potential color / scheme input discrepancies
if ~colorCat(sc, fc, nc)
    if strcmpi(sc, 'series') && size(fc, 1) ~= np
        warning(['series scheme was indicated, but a color matrix was either unspecified or ' ...
                 'dim 1 of color matrix didnt match the number of data series(' num2str(np) ').' ...
                 ' making new colors as default.'])
        fc= distinguishable_colors(np*2);
        fc= fc(np+1:end, :);
        
    elseif strcmpi(sc, 'category') 
        warning(['category scheme was indicated, but dim 1 of color matrix didnt match the number ' ...
                 'of categories(' num2str(nc) '). making new colors as default.'])
        fc= distinguishable_colors(nc*2);
        fc= fc(nc+1:end, :);
    end
end

% compute orientation dependent coordinates, text scaling, and axis limits
[h, k, r, R, txt_r, ax_limits]= getOriDependentCoords(ori, h, k, r, R, np); 


% make arc start and finish points
a0= deg2rad( a * ones(size(data, np_dim), 1) );   % initialize start point
p360= data .* 360;                      % data as proportion of 360 degrees
af= deg2rad( cumsum(p360, nc_dim) + a ); 
a0(:, 2:nc)= af(:, 1:nc-1); 

r0= repmat(r, 1, nc); 
rf= repmat(R, 1, nc); 


% check for zeros
test_mat= [zeros(np, 1) af];
zpos=  diff(test_mat, [], nc_dim) == 0; 
zpres= any(any(zpos));

% compute theta to evaluate arc for max non-zero percentage
inc= zeros(np, 1); 
for n= 1:np
    ix_nonz{n}= find(~zpos(n, :)); 
    [~, loc]= max(data(n, ix_nonz{n})); % index the max non-zero comp
    mnz(n)= ix_nonz{n}(loc);             % to prevent too small increment
    thetas{n, mnz(n)}= linspace(a0(n, mnz(n)), af(n, mnz(n)), patch_res)'; 
    arcrds{n, mnz(n)}= [repmat(r0(n, mnz(n)), patch_res, 1), ...
                        repmat(rf(n, mnz(n)), patch_res, 1)];
    inc(n, 1)=   mean(diff(thetas{n, mnz(n)}));
end

% compute theta for remaining percentages using comparable increment
for n= 1:np
    rest= find(~ismember(1:nc, mnz(n)));
    for p= rest
        n_pts= length( a0(n, p):inc(n):af(n, p) ); 
        thetas{n, p}= linspace(a0(n, p), af(n, p), n_pts)';
        arcrds{n, p}= [repmat(r0(n, p), n_pts, 1), ... 
                       repmat(rf(n, p), n_pts, 1)];
    end
end

%--------------------------------------------------------------------------
% handle unlikely cases in the data (e.g. zeros, one component is 100, etc)
[thetas, perfect_circle]= handleSpecialZeroCases(data, thetas, zpres, zpos, mnz); 


% compute circular arcs
cartConv= @(t, r) pol2cart(t, r); 
for n= 1:np
    [x(n, :), y(n, :)]= cellfun(cartConv, thetas(n, :), arcrds(n, :), 'UniformOutput', false);
    x(n, :)= cellfun(@(t) t + h(n), x(n , :), 'UniformOutput', false);
    y(n, :)= cellfun(@(r) r + k(n), y(n , :), 'UniformOutput', false);
end

% reshape all to single column and "close" the data
x_vtx= cellfun(@(t) [t(:, 1); flipud(t(:, 2)); t(1)], x, 'UniformOutput', false); 
y_vtx= cellfun(@(r) [r(:, 1); flipud(r(:, 2)); r(1)], y, 'UniformOutput', false); 


% pre-allocate structure of graphics objects for each series
if ~perfect_circle
    for s= np:-1:1
        arcs(s).series= gobjects(1, nc);
        lbls(s).series= gobjects(1, nc); 
    end
end

% pre-allocate text coordinates & determine position
[xc, yc]= deal(zeros(np, nc));
centerTheta= cellfun(@median, thetas); 
[xt, yt]= pol2cart(centerTheta, txt_r); 


% plot arcs & percentages
for n= 1:np
    % assign color matrix based on rule
    if strcmpi(sc, 'series')
        colors= repmat(fc(n, :), nc, 1) .* dark; 
    elseif strcmpi(sc, 'category')
        colors= fc;
    end
    % plot elements, skipping zeros (empty cells)
    for j= 1:nc
        if zpres
            if isempty(thetas{n, j})
                continue
            end
        end

        xc(n, j)= xt(n, j) + h(n);
        yc(n, j)= yt(n, j) + k(n); 

        % determine text alignment based on orientation
        if strcmpi(ori, 'concentric')
            halign= 'center'; valign= 'middle'; 
        else
            [halign, valign]= getAlignmentFromAngle(centerTheta(n, j)); % from pie.m
        end

        % a working polyshape function -- NOTE: cannot pre-allocate
        % gobjects for this
        if perfect_circle
            arcs(n).series(1, j)= polyshape([x_vtx{n, j} y_vtx{n, j}]);   
            plot(arcs(n).series(1, j), 'FaceColor', colors(j, :), 'EdgeColor', ec, 'LineWidth', lw); hold on
        else
            arcs(n).series(1, j)= patch('Faces',    1:length(x_vtx{n, j}), ...
                                        'Vertices', [x_vtx{n, j} y_vtx{n, j}], ...
                                        'FaceVertexCData', colors(j, :), ...
                                        'FaceColor', 'flat', 'EdgeColor', ec, 'LineWidth', lw);   hold on
        end

        lbls(n).series(1, j)= text(xc(n, j), yc(n, j), txt(n, j), ...
                                   'color', tc, ...
                                   'HorizontalAlignment', halign, ...
                                   'VerticalAlignment', valign);
    end
end

% position text labels on top of arcs
for n= 1:np
    uistack(lbls(n).series, 'top')
end

set(gcf, 'color', 'w')
axis(ax_limits) 

if np == 1
    axis square
else
    axis image; 
end
axis off; 

H.arcH=  arcs;
H.txtH=  lbls; 
H.color= colors; 

end


%% Helper functions--------------------------------------------------------

function [h, k, r, R, text_radius, ax_limits]= getOriDependentCoords(ori, h, k, r, R, np)

% r = inner radius
% R = outer radius

switch ori
    case 'concentric'    % fix plot at origin & increment radius
        h= zeros(1, np); 
        k= zeros(1, np); 
        r= (2:np+1)';    % override any innerRadius argument
        R= r + 0.975; 
    otherwise            % scale text pos differently for multiple series
        r= repmat(r, np, 1); 
        R= repmat(R, np, 1); 

        % adjust (h, k) for the desired orientation
        if strcmpi(ori, 'horizontal')
            h= h:(3*R+1):(3*R+1)*np-1;       % step right for horz series
            k= zeros(1, np); 
        else
            h= zeros(1, np); 
            k= k:-(3*R+1):-((3*R+1)*np-1);   % step down for vert series
        end
end

% pre-set ax limits
x1= h-R-np;   x2= h+R+np;    y1= k-R-np;    y2= k+R+np;
if strcmpi(ori, 'concentric')
    text_radius= (R + r) / 2;  % average inner & outer radii
    x_lo= x1(np); x_hi= x2(np); y_lo= y1(np); y_hi= y2(np); 
else
    text_radius= R * 1.05;     % position outside the rings
    x_lo= x1(1);  x_hi= x2(np);  y_lo= y1(np);  y_hi= y2(1); 
end

ax_limits= [x_lo x_hi y_lo y_hi]; 

end

%--------------------------------------------------------------------------
function [thetas, fill100]= handleSpecialZeroCases(data, thetas, zpres, zpos, max_nonzero_ix)

fill100= false; 

if zpres    % handle case where zeros are present
    [z_row, z_col]= find(zpos);     % index zeros
    % if data already sums to 1 (i.e., one component accounts for 100%),
    % then complete the circle
    % else, skip the percentage, producing a partial pie
    [zr, ix]= sort(z_row);
    zc= z_col(ix); 
    for i= 1:length(zr)             % for rows with zeros
        if sum(data(zr(i), :)) == 1  
            fill100= true; 

            if zc(i) > max_nonzero_ix(zr(i))
                thetas{zr(i), zc(i)-1}(end)= thetas{zr(i), zc(i)-1}(1);
            else
                thetas{zr(i), zc(i)+1}(end)= thetas{zr(i), zc(i)+1}(1);
            end
        else
            thetas{zr(i), zc(i)}= repmat(thetas{zr(i), zc(i)}, 1, 2); 
        end
    end
end

end

%--------------------------------------------------------------------------
function [halign, valign] = getAlignmentFromAngle(angle)
% Determine the text label alignment based on the angle around the circle.

% Convert the angle to degrees
angle = (180/pi)*angle;

% Round the angles to the nearest 45 degrees
angle = mod(round(angle/45)*45, 360);

% Determine the horizontal alignment
if angle == 90 || angle == 270
    halign = 'center';
elseif angle > 90 && angle < 270
    halign = 'right';
else
    halign = 'left';
end

% Determine the vertical alignment
if angle == 0 || angle == 180
    valign = 'middle';
elseif angle > 180 && angle < 360
    valign = 'top';
else
    valign = 'bottom';
end

end
