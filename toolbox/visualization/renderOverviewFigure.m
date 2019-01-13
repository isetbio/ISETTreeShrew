function hFig = renderOverviewFigure(figNo, inputImageArray, sceneArray, ...
    oiArray, coneExcitationArray, demosaicedExcitationArray, theConeMosaic, timeAxis)
% Summarize the various ISETBio pipeline steps for a set of lemem images
%
% Syntax:
%   renderOverviewFigure(figNo, inputImageArray, sceneArray, oiArray, coneExcitationArray, demosaicedExcitationArray, theConeMosaic)
%
% Description:
%    Summarize the various ISETBio pipeline steps for a set of lemem images
%
%
% Inputs:
%    figNo                      Figure number, 1, 2 etc.
%    inputImageArray            Array with the input images to be visualized
%    sceneArray                 Array with the scenes to be visualized
%    oiArray                    Array with the ois to be visualized
%    coneExcitationArray        Array with the cone mosaic responses to be visualized
%    demosaicedExcitationArray  Array with the demosaiced cone responses to be visualized
%    theConeMosaic              The employed cone mosaic
%
%
% Outputs:
%    hFig                       Handle to the generated figure
%

% History:
%    1/24/2019  NPC   Wrote it


    hFig = figure(figNo); 
    clf;
    set(hFig, 'Position', [10 10 1450 1500], 'Color', [1 1 1]);
    
    subplotRows = numel(inputImageArray);
    if (numel(demosaicedExcitationArray) == 0)
        subplotCols = 4;
    else
        subplotCols = 6;
    end
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
    
    fontSize = 14;
    ticksDegs = [-10:1:10];
    spatialRangeDegs = max(theConeMosaic.fov)/2*[-1 1];
    
    for imIndex = 1:numel(inputImageArray)
        subplot('Position', subplotPosVectors(imIndex,1).v);
        [mRows, nCols, ~] = size(inputImageArray{imIndex});
        image(1:nCols, 1:mRows, inputImageArray{imIndex}); axis 'image';
        xlabel('\it space (pixels)');
        set(gca, 'FontSize', fontSize);
        if (imIndex==1)
            title('input RGB');
        end
        
        realizedSceneRGB = sceneGet(sceneArray{imIndex}, 'rgb');
        spatialSupportMilliMeters = sceneGet(sceneArray{imIndex}, 'spatial support', 'mm');
        spatialSupportX = squeeze(spatialSupportMilliMeters(1,:,1));
        spatialSupportY = squeeze(spatialSupportMilliMeters(:,1,2));
        
        subplot('Position', subplotPosVectors(imIndex,2).v);
        image(spatialSupportX, spatialSupportY, realizedSceneRGB);
        axis 'image'; axis 'ij'; grid on
        set(gca, 'FontSize', fontSize, 'YTickLabels', {}, 'YTick', [-5:1:5], ...
            'XTick', [-5:1:5], 'XTickLabels', sprintf('%2.1f\n', -5:1:5));
        if (imIndex==1)
            title('ISETBio scene');
        end
        
        if (imIndex == 4)
            xlabel('\it space (mm on presentation display)');
       else
            set(gca, 'XTickLabel', {});
            xlabel('');
        end
    
        oiRGB = oiGet(oiArray{imIndex}, 'rgb');
        spatialSupportMM = oiGet(oiArray{imIndex}, 'spatial support', 'mm');
        optics = oiGet(oiArray{imIndex}, 'optics');
        focalLength = opticsGet(optics, 'focal length');
        mmPerDegree = focalLength*tand(1)*1e3;
        spatialSupportDegs = spatialSupportMM/mmPerDegree;
        spatialSupportXDegs = spatialSupportDegs(1,:,1);
        spatialSupportYDegs = spatialSupportDegs(:,1,2);
        
        
        subplot('Position', subplotPosVectors(imIndex,3).v);
        image(spatialSupportXDegs,spatialSupportYDegs, oiRGB);
        axis 'image'; axis 'ij'; grid on
        set(gca, 'FontSize', fontSize, ...
            'XLim', spatialRangeDegs, 'YLim', spatialRangeDegs, ...
            'XTick', ticksDegs, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksDegs, 'YTickLabels', {} ...
            );
        if (imIndex == 4)
            xlabel('\it space (degs on retina)');
        else
            set(gca, 'XTickLabel', {});
            xlabel('');
        end
        
        if (imIndex==1)
            title('retinal image');
        end
        
        % The cone mosaic excitation
        axHandle = subplot('Position', subplotPosVectors(imIndex,4).v);
        theConeExcitationsAllInstances = coneExcitationArray{imIndex};

        % Mean over all trials
        meanResponse = theConeMosaic.reshapeHex2DmapToHex3Dmap(squeeze(mean(theConeExcitationsAllInstances,1)));
        
        % Determine response range
        responseRange = [min(meanResponse(:)) max(meanResponse(:))];
        timeAxis = timeAxis - mean(timeAxis);
        [~, centerTimeBin] = min(abs(timeAxis));
        meanResponse = squeeze(meanResponse(:,:,centerTimeBin));

        theConeMosaic.renderActivationMap(axHandle, meanResponse, ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 responseRange(2)], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', theConeMosaic.integrationTime*1000));
        
        ticksMeters = ticksDegs * mmPerDegree * 1e-3;
        spatialRangeMeters = spatialRangeDegs * mmPerDegree * 1e-3;
        axis 'image'; axis 'ij'; 
        set(gca, 'FontSize', fontSize, ...
            'XLim', spatialRangeMeters, 'YLim', spatialRangeMeters, ...
            'XTick', ticksMeters, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksMeters, 'YTickLabels', {} ...
            );
        
        if (imIndex == 4)
            xlabel('\it space (degs on retina)');
        else
            set(gca, 'XTickLabel', {});
            xlabel('');
        end
        
        ylabel('');
        if (imIndex==1)
            title('mosaic activation');
        end
        
        if (numel(demosaicedExcitationArray)>0)
            % The demosaiced response maps
            dStruct = demosaicedExcitationArray{imIndex};

            % The demosaiced L-cone excitation
            subplot('Position', subplotPosVectors(imIndex,5).v);
            imagesc(dStruct.demosaicedMapSupportDegs, ...
                    dStruct.demosaicedMapSupportDegs, ...
                    dStruct.theDemosaicedLconeIsomerizationMap);
            set(gca, 'FontSize', fontSize, 'XTick', ticksDegs, 'YTickLabel', {});
            axis 'image'
            colormap(gray);
            if (imIndex == 4)
                xlabel('\it space (degs on retina)');
            else
                set(gca, 'XTickLabel', {});
                xlabel('');
            end

            ylabel('');
            if (imIndex==1)
                title(sprintf('de-mosaiced\nL-cone activation'));
            end
        
            % The demosaiced S-cone excitation
            subplot('Position', subplotPosVectors(imIndex,6).v);
            imagesc(dStruct.demosaicedMapSupportDegs, ...
                    dStruct.demosaicedMapSupportDegs, ...
                    dStruct.theDemosaicedSconeIsomerizationMap);
            set(gca, 'FontSize', fontSize, 'XTick', ticksDegs, 'YTickLabel', {});
            axis 'image'
            colormap(gray);
            if (imIndex == 4)
                xlabel('\it space (degs on retina)');
            else
                set(gca, 'XTickLabel', {});
                xlabel('');
            end

            ylabel('');
            if (imIndex==1)
                title(sprintf('de-mosaiced\nS-cone activation'));
            end
        end
        
        drawnow;
    end % imIndex
end
