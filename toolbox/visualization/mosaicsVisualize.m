function mosaicsVisualize(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits)
% Co-visualize a treeshrew and a human cone mosaic
%
% Syntax:
%   mosaicsVisualize(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits)
%
% Description:
%    Co-visualize a treeshrew and a human cone mosaic.
%
%
% Inputs:
%    figNo                      Figure number, 1, 2 etc.
%    theTreeShrewConeMosaic     The treeshrew @coneMosaicHex object 
%    theHumanConeMosaic         The human @coneMosaicHex object 
%    displayUnits               The spatial support units, 'degs' or 'microns'
%
% Outputs:
%    None
%

% History:
%    1/24/2019  NPC   Wrote it

    
    showPortionOfMosaic = true;
    if (showPortionOfMosaic)
        if (strcmp(displayUnits, 'degs'))
            spatialRangeDegs = [-1 1];
            ticksDegs = [-2:0.5:2];
            spatialRangeMetersTreeShrew = spatialRangeDegs * theTreeShrewConeMosaic.micronsPerDegree * 1e-6;
            spatialRangeMetersHuman = spatialRangeDegs * theHumanConeMosaic.micronsPerDegree * 1e-6;
            ticksMetersTreeShrew = ticksDegs * theTreeShrewConeMosaic.micronsPerDegree * 1e-6;
            ticksMetersHuman = ticksDegs * theHumanConeMosaic.micronsPerDegree * 1e-6;
        else
            spatialRangeMicrons = theTreeShrewConeMosaic.micronsPerDegree*[-2 2];
            ticksMicrons = [-500:100:500];
            spatialRangeMetersTreeShrew = spatialRangeMicrons * 1e-6;
            spatialRangeMetersHuman = spatialRangeMicrons * 1e-6;
            ticksMetersTreeShrew = ticksMicrons * 1e-6;
            ticksMetersHuman = ticksMicrons * 1e-6;
        end
    end
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 2255 1130], 'Color', [1 1 1]);
    
    subplotRows = 1;
    subplotCols = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
   
    ax = subplot('Position', subplotPosVectors(1,1).v);
    if (strcmp(displayUnits, 'degs'))
        theTreeShrewConeMosaic.visualizeGrid('axesHandle', ax, ...
            'ticksInVisualDegs', true, ...
            'backgroundColor', [1 1 1]); 
         if (showPortionOfMosaic)
            set(gca, 'FontSize', 18, ...
            'XTick', ticksMetersTreeShrew, 'XTickLabel', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksMetersTreeShrew, 'YTickLabel', {}, ...
            'XLim', spatialRangeMetersTreeShrew, 'YLim', spatialRangeMetersTreeShrew ...
            );
         end
    else
        theTreeShrewConeMosaic.visualizeGrid('axesHandle', ax, ...
            'ticksInMicrons', true, ...
            'backgroundColor', [1 1 1]); 
        if (showPortionOfMosaic)
            set(gca, 'FontSize', 18, ...
            'XTick', ticksMetersTreeShrew, 'XTickLabel', sprintf('%2.0f\n', ticksMicrons), ...
            'YTick', ticksMetersTreeShrew, 'YTickLabel', {}, ...
            'XLim', spatialRangeMetersTreeShrew, 'YLim', spatialRangeMetersTreeShrew ...
            );
         end
    end
    
    title(sprintf('treeshrew mosaic %d cones', numel(find(theTreeShrewConeMosaic.pattern>1))));
    ylabel('');
    
    ax = subplot('Position', subplotPosVectors(1,2).v);
    if (strcmp(displayUnits, 'degs'))
        theHumanConeMosaic.visualizeGrid('axesHandle', ax, ...
            'ticksInVisualDegs', true, ...
            'backgroundColor', [1 1 1], 'generateNewFigure', true); 
        if (showPortionOfMosaic)
            set(gca, 'FontSize', 18, ...
            'XTick', ticksMetersHuman, 'XTickLabel', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksMetersHuman, 'YTickLabel', {}, ...
            'XLim', spatialRangeMetersHuman, 'YLim', spatialRangeMetersHuman ...
            );
        end
    else
        theHumanConeMosaic.visualizeGrid('axesHandle', ax, ...
            'ticksInMicrons', true, ...
            'backgroundColor', [1 1 1], 'generateNewFigure', true); 
        if (showPortionOfMosaic)
            set(gca, 'FontSize', 18, ...
            'XTick', ticksMetersHuman, 'XTickLabel', sprintf('%2.0f\n', ticksMicrons), ...
            'YTick', ticksMetersHuman, 'YTickLabel', {}, ...
            'XLim', spatialRangeMetersHuman, 'YLim', spatialRangeMetersHuman ...
            );
         end
    end
    title(sprintf('human mosaic %d cones', numel(find(theHumanConeMosaic.pattern>1))));
    
    ylabel('');
end

