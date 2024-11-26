%SUMMARY
% Author: Connor Gallimore
% 12/26/2023

% This function creates arcs of a circle that represent percentages as the
% proportion of total circumference. I first saw this data visualization
% idea in a Nature paper (below) and wanted to implement it here in MATLAB.

% (Fig 4, Soula et al. 2023, https://www.nature.com/articles/s41593-023-01270-2)

% Required arguments:
    % 'data', a vector or matrix of proportions/percentages. 

% Optional positional argument (2nd input after data)
    % 'dim',  the dimension containing each sub-cat/ component of the total
    %         e.g. if data series are distributed along the rows, such that
    %              sub-categories / components are arranged along the 
    %              columns, dim argument would be '2'
    %              if data series are distributed along columns and
    %              sub-cats along the rows, your dim arg would be '1'
    %         if input is scalar or vector, this input is not required
    %         if input is matrix, the function defaults to dim 2 
    %              (i.e. assumes different groups/series are row-wise),
    %              which passing this argument overrides. 

% Optional Name,Value pairs:
    % 'facecolor', an m x 3 vector or matrix specifying RGB triplet(s), or
    %              a 1 x m array setting face color values (e.g. 1:m)
    % 'edgecolor', same as 'facecolor', but only specify one that will be
    %              applied to all patch edges (e.g. 'w', or 'k')
    % 'textcolor', same as 'edgecolor'
    % 'facealpha', scalar in range [0, 1] specifying opacity
    % 'linewidth', a positive scalar value in points (1 pt = 1/72 inches)
    % 'orientation', 'horizontal', 'vertical', or 'concentric' 
    %                this argument only exerts effects for more than 1 data 
    %                series, determining whether they are plotted from 
    %                left-to-right, top-to-bottom, or 'inside-to-outside'
    % 'precision', specifies the rounding precision for text labels (i.e. 
    %              the max number of decimal places)
    % 'innerRadius', scalar in range [0, 1] specifying the inner radius of
    %                patches as a proportion of the outer radius.
    %                a value of '0' creates a pie chart
    %                a value of '1' creates a ring with no visible slices
    %                default = 0.65
    % 'outerRadius', a non-negative scalar specifying outer patch radius
    % 'startAngle',  scalar value in degrees specifying start angle where
    %                patches emanate. 0 degrees corresponds to 3 o'oclock.
    %                positive values rotate counterclockwise, negative
    %                values clockwise. 
    % 'scheme', 'category' or 'series', determines the coloring scheme. In
    %           one case, your 'color' matrix may represent the color you 
    %           want each common 'category' to be for all series (default).
    %           In another case, you may be specifying the base color you
    %           want your percentage components to be for each 'series'. In
    %           this option, subsequent percentages will be plotted darker
    % 'patchRes',  a positive scalar value specifying the number of points
    %              used to represent the largest arc. Subsequent (smaller) 
    %              arcs are represented in proportion to this


% Outputs:
% a structure of handles 'H' to the line 'arcs', text labels ('lbls'), and
% 'colors'

% Examples: 
% see 'donutPlot_demo.mlx' for usage tips and tricks

% Additional notes:
% Used 'getAlignmentFromAngle' from MATLAB's 'pie.m', all rights go to
% Clay M. Thompson 3-3-94 and The MathWorks, Inc., Copyright 1984-2022 

% Also used 'distinguishable_colors.m' by Timothy E. Holy, for which 
% all rights reserved to them (Copyright (c) 2010-2011, 
% https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
%--------------------------------------------------------------------------


function H = donutPlot(data, varargin)

% static defaults -- center origin (h, k)
h= 0;
k= 0;

% input parser defaults
def.R=  10;             % outer radius
def.r=  0.65;           % inner radius as a proportion of outer radius
def.d=  2;              % dim to operate on
def.FC= [];             % face color; resolved in local function
def.FA= 1;              % face alpha
def.EC= 'k';            % edge color
def.TC= 'k';            % text color
def.LW= 1;              % line width
def.SA= 0;              % start angle 3 o'clock ('90' would be midnight)
def.RP= 2;              % rounding precision 
def.PR= 300;            % resolution of 1 arc of largest patch
def.OR= 'horizontal';   % orientation (has no effect if vector input)
def.CS= 'category';     % color scheme
def.ori_types= {'horizontal', 'vertical', 'concentric', 'h', 'v', 'c'};
def.sch_types= {'category', 'series'}; 

% pass defaults to parser obj for validation
[p, f]= validateInputs(data, def, varargin);

% assign parsed user inputs
dim=  p.Results.dim;
R=    p.Results.outerRadius; 
r=    R * p.Results.innerRadius;
a=    p.Results.startAngle; 
ori=  p.Results.orientation; 
sc=   p.Results.scheme; 
col=  p.Results.faceColor; 
fa=   p.Results.faceAlpha;
lw=   p.Results.lineWidth;
ec=   p.Results.edgeColor;
tc=   p.Results.textColor; 
prec= p.Results.precision; 
res=  p.Results.patchRes; 

% force vector inputs to be row-wise
if f.validArray(data)
    data= data(:)'; 
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

% only used if 'series' scheme is requested
dark= linspace(1, 0.1, nc)';  % darkness scaling

% normalize groups that don't sum to 1
if any(sum(data, dim) > 1+sqrt(eps))
    ix2n= sum(data, dim) > 1+sqrt(eps);  % indices to normalize
    data(ix2n, :)= data(ix2n, :) ./ sum(data(ix2n, :), dim);
end

% pre-allocate cells for arcs and cartesian coordinates
[thetas, arcrds, x, y]= deal(cell(np, nc)); 

% pre-allocate txt strings for later based on rounding precision
txt= strcat(string(round(data .* 100, prec)), repmat("%", np, nc)); 

% correct potential mismatches in user facecolor / scheme input
fc= resolveColorSchemeMismatch(col, sc, np, nc, f);

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
    thetas{n, mnz(n)}= linspace(a0(n, mnz(n)), af(n, mnz(n)), res)'; 
    arcrds{n, mnz(n)}= [repmat(r0(n, mnz(n)), res, 1), ...
                        repmat(rf(n, mnz(n)), res, 1)];
    inc(n, 1)=   mean(diff(thetas{n, mnz(n)}));
end

% compute theta for remaining percentages using comparable increment
for n= 1:np
    rest= find(~ismember(1:nc, mnz(n)));
    for c= rest
        n_pts= length( a0(n, c):inc(n):af(n, c) ); 
        thetas{n, c}= linspace(a0(n, c), af(n, c), n_pts)';
        arcrds{n, c}= [repmat(r0(n, c), n_pts, 1), ... 
                       repmat(rf(n, c), n_pts, 1)];
    end
end


% handle unlikely cases in the data (e.g. zeros, one component is 100, etc)
[thetas, perfect_circle]= handleSpecialZeroCases(data, thetas, zpres, zpos, mnz); 


% finally compute the circular arcs
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

% define text coordinates 
centerTheta= cellfun(@median, thetas); 
[xt, yt]= pol2cart(centerTheta, txt_r); 

% determine text alignment based on orientation
switch ori
    case {'concentric', 'c'}
        [halign, valign]= deal( repmat({'center'}, np, nc), repmat({'middle'}, np, nc) ); 
    otherwise    
        [halign, valign]= cellfun(@getAlignmentFromAngle, num2cell(centerTheta), 'UniformOutput', false); 
end

% create facealpha matrix
if isscalar(fa)
    fa= repmat(fa, size(data)); 
end

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

        if perfect_circle   % use polyshape
            arcs(n).series(1, j)= polyshape([x_vtx{n, j} y_vtx{n, j}]);   
            plot(arcs(n).series(1, j), 'FaceColor', colors(j, :), 'FaceAlpha', fa(n, j), ...
                                       'EdgeColor', ec, 'LineWidth', lw); hold on
        else                % use patch
            arcs(n).series(1, j)= patch('Faces',    1:length(x_vtx{n, j}), ...
                                        'Vertices', [x_vtx{n, j} y_vtx{n, j}], ...
                                        'FaceVertexCData', colors(j, :), ...
                                        'FaceVertexAlphaData', fa(n, j), ...
                                        'FaceAlpha', 'flat', 'FaceColor', 'flat', ...
                                        'EdgeColor', ec, 'LineWidth', lw);   hold on
        end

        lbls(n).series(1, j)= text(xt(n, j) + h(n), yt(n, j) + k(n), txt(n, j), ...
                                   'color', tc, ...
                                   'HorizontalAlignment', halign{n, j}, ...
                                   'VerticalAlignment', valign{n, j});
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

function [p, f]= validateInputs(data, defaults, varargs)

% returns structure of parsed inputs and some validation functions

% define anonymous validation functions in a structure
f.scalarNum=  @(x) isnumeric(x) && isscalar(x);
f.isNonNeg=   @(x) all(all(x >= 0)); 
f.validData=  @(x) all((isnumeric(x) | islogical(x))) && isreal(x) && f.isNonNeg(x); 
f.validDim=   @(x) f.scalarNum(x) && f.isNonNeg(x);
f.validArray= @(x) ~isscalar(x) && isvector(x) && isnumeric(x);
f.validAlpha= @(x) (isscalar(x) && (x >= 0) && (x <= 1)) || (ismatrix(x) && ~isscalar(x));
f.validProps= @(x) isnumeric(x) && ismatrix(x) && all(f.isNonNeg(x) & all(x <= 1));
f.validRGB=   @(x) f.validProps(x) && size(x, 2) == 3;
f.validColV=  @(x) f.validArray(x) && all(x >= 1); 
f.validColor= @(x) f.validRGB(x) || ischar(x);
f.validFCol=  @(x) f.validRGB(x) || f.validColV(x) || iscellstr(x) || isstring(x); 
f.colorCat=   @(x, y, z) strcmpi(x, 'category') && size(y, 1) == z; 

p= inputParser; 

addRequired(p, 'data', f.validData);
addOptional(p, 'dim', defaults.d, f.validDim);   
addParameter(p, 'scheme', defaults.CS, @(x) any(validatestring(x, defaults.sch_types)))
addParameter(p, 'faceAlpha', defaults.FA, f.validAlpha);
addParameter(p, 'faceColor', defaults.FC, f.validFCol);
addParameter(p, 'edgeColor', defaults.EC, f.validColor);
addParameter(p, 'textColor', defaults.TC, f.validColor);
addParameter(p, 'lineWidth', defaults.LW, f.validDim);
addParameter(p, 'orientation', defaults.OR, @(x) any(validatestring(x, defaults.ori_types)));
addParameter(p, 'precision', defaults.RP, f.validDim);
addParameter(p, 'innerRadius', defaults.r, f.validAlpha);
addParameter(p, 'outerRadius', defaults.R, f.validDim);
addParameter(p, 'startAngle', defaults.SA, f.scalarNum);
addParameter(p, 'patchRes', defaults.PR, f.validDim);

% parse and assign all inputs
parse(p, data, varargs{:});

end

%--------------------------------------------------------------------------
function facecolor= resolveColorSchemeMismatch(color, sc, n_pcts, n_cats, fxns)

% parse color / color scheme args
if isempty(color)
    if strcmpi(sc, 'series')
        facecolor= distinguishable_colors(n_pcts*2);
        facecolor= facecolor(n_pcts+1:end, :);
    else
        facecolor= distinguishable_colors(n_cats*2);
        facecolor= facecolor(n_cats+1:end, :);
    end
else
    if fxns.validRGB(color)
        facecolor= color; 
    elseif fxns.validColV(color) || isstring(color)
        facecolor= color(:);
    elseif iscellstr(color)
        facecolor= vertcat(color{:});
    end

end

% resolve potential color / scheme input discrepancies
if ~fxns.colorCat(sc, facecolor, n_cats)
    if strcmpi(sc, 'series') && size(facecolor, 1) ~= n_pcts
        warning(['series scheme was indicated, but a color matrix was either unspecified or ' ...
                 'dim 1 of color matrix didnt match the number of data series(' num2str(n_pcts) ').' ...
                 ' making new colors as default.'])
        facecolor= distinguishable_colors(n_pcts*2);
        facecolor= facecolor(n_pcts+1:end, :);
        
    elseif strcmpi(sc, 'category') 
        warning(['category scheme was indicated, but dim 1 of color matrix didnt match the number ' ...
                 'of categories(' num2str(n_cats) '). making new colors as default.'])
        facecolor= distinguishable_colors(n_cats*2);
        facecolor= facecolor(n_cats+1:end, :);
    end
end

end

%--------------------------------------------------------------------------
function [h, k, r, R, text_radius, ax_limits]= getOriDependentCoords(ori, h, k, r, R, np)

% r = inner radius
% R = outer radius

r= repmat(r, np, 1); 
R= repmat(R, np, 1); 

switch ori
    case {'concentric', 'c'}    % fix plot at origin & increment radius
        h= zeros(1, np); 
        k= zeros(1, np); 
        r= (2:np+1)';           % override any innerRadius argument
        R= r + 0.975; 
    case {'horizontal', 'h'}    % step right for horz series (h)
        h= h:(3*R+1):(3*R+1)*np-1;       
        k= zeros(1, np); 
    case {'vertical', 'v'}      % step down for vert series (k)
        h= zeros(1, np); 
        k= k:-(3*R+1):-((3*R+1)*np-1);   
    otherwise
end

% pre-set ax limits
x1= h-R-np;   x2= h+R+np;    y1= k-R-np;    y2= k+R+np;

switch ori
    case {'concentric', 'c'}
        text_radius= (R + r) / 2;  % average inner & outer radii
        x_lo= x1(np); x_hi= x2(np); y_lo= y1(np); y_hi= y2(np); 
    otherwise
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
elseif isscalar(data) 
    fill100= true;
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
