[![DOI](https://zenodo.org/badge/736407888.svg)](https://zenodo.org/doi/10.5281/zenodo.10460350) 
[![View circPercent on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/157171-circpercent)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=cgallimore25/circPercent&file=README.md)

# donutPlot
Use circular arcs to graphically illustrate proportions and percentages in MATLAB.

Catalog of changes in version 2.0.0 (formerly 'circPercent'):
- Converted plotted elements from line to patch objects for both fill and edge customization
- Support for concentric plotting style, embedding data series inside-out
- New name, value pairs for user-defined start angle, donut radius, transparency, linewidth, and text color
- Transitioned to an input parser object for flexible argument handling
- Modularization of coordinate calculations and special cases


## Syntax:
**donutPlot(data)**

**donutPlot(data, Name, Value, ...)**

**h = donutPlot(_)**

## Input Arguments:
*(Required)*

- **data** - A scalar, vector, or matrix of proportions/percentages. The rows are assumed to be groups, or 'series' of data, and the columns are the components to the total.
             [scalar | vector | matrix]

## Output Arguments:
*(Optional)*
- **h**                - Figure handle of spider plot.
                         [structure]

## Name-Value Pair Arguments:
*(Optional)*

-**'facecolor'**   - an m x 3 vector or matrix specifying RGB triplet(s), or a 1 x m array setting face color values (e.g. 1:m), or a cell of character arrays or string array of valid MATLAB color names / short names (e.g. {'r', 'g', 'b'}).
                     [m x 3 vector | matrix | cellstr | string]

-**'edgecolor'**   - same as 'facecolor', but only specify one that will be applied to all patch edges (e.g. 'w', or 'k').
                     [RGB triplet | char]

-**'textcolor'**   - same as 'edgecolor'
                     [RGB triplet | char]

-**'facealpha'**   - scalar in range [0, 1] specifying opacity.

-**'linewidth'**   - a positive scalar value in points (1 pt = 1/72 inches).

-**'orientation'** - 'horizontal', 'vertical', or 'concentric'; this argument only exerts effects for more than 1 data series, determining whether they are plotted from 
                     left-to-right, top-to-bottom, or 'inside-to-outside'. Can also use shorthand forms 'h', 'v', 'c'.
                     ['horizontal' | 'vertical' | 'concentric' | 'h' | 'v' | 'c']

-**'precision'**   - specifies the rounding precision for text labels (i.e. the max number of decimal places to keep)
                     [0 | 1 | ... n | non-negative scalar]

-**'innerRadius'** - scalar in range [0, 1] specifying the inner radius of patches as a proportion of the outer radius. A value of '0' creates a pie chart; a value of '1' creates a ring with no visible slices 
                     [0.65 (default) | scalar]

-**'outerRadius'** - a non-negative scalar specifying outer patch radius. Most applicable for 'horizontal' or 'vertical', whereas 'concentric' defaults to ring widths of radius 1 and this argument is ignored.

-**'startAngle'**  - scalar value in degrees specifying start angle where patches emanate. 0 degrees corresponds to 3 o'oclock. Positive values rotate counterclockwise, negative values clockwise.
                     [0 (default) | scalar]

-**'scheme'**      - 'category' or 'series', determines the coloring scheme. In one case, your 'color' matrix may represent the color you want each common 'category' to be for all series (default). 
                     In another case, you may be specifying the base color you want your percentage components to be for each 'series'. In this option, subsequent percentages will be plotted darker.
                     ['category' | 'series']

-**'patchRes'**    - a positive scalar value specifying the number of points used to represent the largest arc. Subsequent (smaller) arcs are represented in proportion to this.
                     [300 (default) | positive scalar]


## Examples

### Example 1: Minimum working example

```matlab
% Create some data
data= [0.2445	0.2554	0.3237	0.1762;
       0.1924	0.2255	0.1672	0.4147;
       0.1513	0.2769	0.2749	0.2967;
       0.1402	0.3639	0.2439	0.2518]; 

total= sum(data, 2);  % show groups are rows and components are cols

% Create plot
figure;
donutPlot(data);
```