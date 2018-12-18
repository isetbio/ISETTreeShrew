function computeMosaicResponses

    % Create presentation display and place it 5 cm in front of the eye
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 5/100);

    imDir = 'imagesLaMem';
    imFileNames = {...
        '00011740' ...
        '00011741' ...
        '00011742' ...
        '00011743' ...
        };
    imageSizePixels = [256 256];
    imageSizeVisualAngle = 3.95;
    meanLuminanceCdM2 = 40;
    coneIntegrationTime = 10/1000;
    
    [inputImageArray, sceneArray] = preProcessImages(imDir, imFileNames,...
        presentationDisplay, ...
        imageSizePixels, ...
        imageSizeVisualAngle, ...
        meanLuminanceCdM2);

    % Generate default treeshrew optics
    theTreeShrewOI = oiTreeShrewCreate();
    
    % Generate default human optics
    theHumanOI = oiCreate('wvf human');
    
    % Load the human cone mosaic
    load('humanConeMosaic4degs.mat', 'theData')
    theHumanConeMosaic = theData;
    theHumanConeMosaic.integrationTime = coneIntegrationTime;
    theHumanConeMosaic.eccBasedMacularPigment = false;
    clear 'theData';
    
   
    % Create same angular size treeshrew cone mosaic
    theTreeShrewConeMosaic = coneMosaicTreeShrewCreate(...
        theTreeShrewOI.optics.micronsPerDegree, ...
        'fovDegs', max(theHumanConeMosaic.fov)*[1 1], ...
        'sConeMinDistanceFactor', 2, ...
        'integrationTime', coneIntegrationTime);

    %visualizeMosaics(theTreeShrewConeMosaic, theHumanConeMosaic);
    
    
    oiArrayTreeShrew = cell(1, numel(sceneArray));
    coneExcitationArrayTreeShrew = cell(1, numel(sceneArray));
    oiArrayHuman = cell(1, numel(sceneArray));
    coneExcitationArrayHuman = cell(1, numel(sceneArray));
    
    demosaicResponses = ~true;
    demosaicingSampleSpacingMicrons = 1;
    
    nRepeats = 1;
    
    doHumanComputation = ~false;
    doTreeShrewComputation = ~true;
    
    if (doHumanComputation)
    tic
    
    for repeat = 1:nRepeats
        fprintf('processing batch %d/%d\n', repeat,nRepeats);
        for imIndex = 1:numel(sceneArray)
            fprintf('\t processing image %d \ %d\n', imIndex, numel(sceneArray));
            % Compute the retinal image
            oiArrayHuman{imIndex} = oiCompute(theHumanOI, sceneArray{imIndex});
            % Compute mosaic excitation
            coneExcitationArrayHuman{imIndex} = theHumanConeMosaic.compute(oiArrayHuman{imIndex});
            if (demosaicResponses)
            % Compute de-mosaiced responses
                coneExcitationsRates = coneExcitationArrayHuman{imIndex}/theHumanConeMosaic.integrationTime;
                [theDemosaicedLconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, lConeXlocsDegs, lConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('L-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                [theDemosaicedMconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theMConeIsomerizations, MConeXlocsDegs, MConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('M-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                [theDemosaicedSconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, sConeXlocsDegs, sConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('S-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);
            end
            
        end % inIndex
    end
    fprintf('Human computation: %f\n', toc/4/nRepeats);
    fprintf('Rendering human figure. This will take a long time\n');
    renderFigure(20, inputImageArray, sceneArray, oiArrayHuman, coneExcitationArrayHuman, theHumanConeMosaic);
    end
    
    
    if (doTreeShrewComputation)
    tic
    for repeat = 1:nRepeats
        fprintf('processing batch %d/%d\n', repeat,nRepeats);
        for imIndex = 1:numel(sceneArray)
            fprintf('\t processing image %d \ %d\n', imIndex, numel(sceneArray));
            % Compute the retinal image
            oiArrayTreeShrew{imIndex} = oiCompute(theTreeShrewOI, sceneArray{imIndex});
            % Compute mosaic excitation
            coneExcitationArrayTreeShrew{imIndex} = theTreeShrewConeMosaic.compute(oiArrayTreeShrew{imIndex});
            
            if (demosaicResponses)
                % Compute de-mosaiced responses
                coneExcitationsRates = coneExcitationArrayHuman{imIndex}/theTreeShrewConeMosaic.integrationTime;
                [theDemosaicedLconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, lConeXlocsDegs, lConeYlocsDegs] = ...
                    theTreeShrewConeMosaic.demosaicConeTypeActivationFromFullActivation('L-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                [theDemosaicedMconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theMConeIsomerizations, MConeXlocsDegs, MConeYlocsDegs] = ...
                    theTreeShrewConeMosaic.demosaicConeTypeActivationFromFullActivation('M-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                [theDemosaicedSconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, sConeXlocsDegs, sConeYlocsDegs] = ...
                    theTreeShrewConeMosaic.demosaicConeTypeActivationFromFullActivation('S-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);
            end
            
        end % inIndex
    end
    fprintf('TreeShrew computation: %f\n', toc/4/nRepeats);
    
    renderFigure(10, inputImageArray, sceneArray, oiArrayTreeShrew, coneExcitationArrayTreeShrew, theTreeShrewConeMosaic);
    end
    
    
end

function visualizeMosaics(theTreeShrewConeMosaic, theHumanConeMosaic)
    spatialRangeDegs = [-2 2];
    ticksDegs = [-2:0.5:2];
    spatialRangeMetersTreeShrew = spatialRangeDegs * theTreeShrewConeMosaic.micronsPerDegree * 1e-6;
    spatialRangeMetersHuman = spatialRangeDegs * theHumanConeMosaic.micronsPerDegree * 1e-6;
    ticksMetersTreeShrew = ticksDegs * theTreeShrewConeMosaic.micronsPerDegree * 1e-6;
    ticksMetersHuman = ticksDegs * theHumanConeMosaic.micronsPerDegree * 1e-6;
    
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1285 745], 'Color', [1 1 1]);
    
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
    theTreeShrewConeMosaic.visualizeGrid('axesHandle', ax, ...
        'displayVisualDegs', true, ...
        'backgroundColor', [1 1 1]); 
    set(gca, 'FontSize', 18, ...
        'XTick', ticksMetersTreeShrew, 'XTickLabel', sprintf('%2.1f\n', ticksDegs), ...
        'YTick', ticksMetersTreeShrew, 'YTickLabel', {}, ...
        'XLim', spatialRangeMetersTreeShrew, 'YLim', spatialRangeMetersTreeShrew ...
        );
    
    title(sprintf('treeshrew mosaic %d cones', numel(find(theTreeShrewConeMosaic.pattern>1))));
    xlabel('space (degs)');
    ylabel('');
    
    ax = subplot('Position', subplotPosVectors(1,2).v);
    theHumanConeMosaic.visualizeGrid('axesHandle', ax, ...
        'displayVisualDegs', true, ...
        'backgroundColor', [1 1 1], 'generateNewFigure', true); 
    set(gca, 'FontSize', 18, ...
        'XTick', ticksMetersHuman, 'XTickLabel', sprintf('%2.1f\n', ticksDegs), ...
        'YTick', ticksMetersHuman, 'YTickLabel', {}, ...
        'XLim', spatialRangeMetersHuman, 'YLim', spatialRangeMetersHuman ...
        );
    title(sprintf('human mosaic %d cones', numel(find(theHumanConeMosaic.pattern>1))));
    xlabel('space (degs)');
    ylabel('');
end


function renderFigure(figNo, inputImageArray, sceneArray, oiArray, coneExcitationArray, theConeMosaic)
    
    subplotRows = numel(inputImageArray);
    subplotCols = 4;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.08, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
   
    idx = find(theConeMosaic.pattern > 1);
    
    ticksDegs = [-10:1:10];
    spatialRangeDegs = max(theConeMosaic.fov)/2*[-1 1];
    
    hFig = figure(figNo); 
    clf;
    set(hFig, 'Position', [10 10 1450 1500], 'Color', [1 1 1]);

    
    for imIndex = 1:numel(inputImageArray)
        subplot('Position', subplotPosVectors(imIndex,1).v);
        [mRows, nCols, ~] = size(inputImageArray{imIndex});
        image(1:nCols, 1:mRows, inputImageArray{imIndex}); axis 'image';
        xlabel('\it space (pixels)');
        set(gca, 'FontSize', 18);
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
        set(gca, 'FontSize', 18, 'YTickLabels', {}, 'YTick', [-5:1:5], ...
            'XTick', [-5:1:5], 'XTickLabels', sprintf('%2.1f\n', -5:1:5));
        if (imIndex==1)
            title('ISETBio scene');
        end
        
        if (imIndex == 4)
            xlabel('\it space (mm on presentation display)');
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
        set(gca, 'FontSize', 18, ...
            'XLim', spatialRangeDegs, 'YLim', spatialRangeDegs, ...
            'XTick', ticksDegs, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksDegs, 'YTickLabels', {} ...
            );
        if (imIndex == 4)
            xlabel('\it space (degs on retina)');
        end
        
        if (imIndex==1)
            title('retinal image');
        end
        
        axHandle = subplot('Position', subplotPosVectors(imIndex,4).v);
        theConeExcitations = coneExcitationArray{imIndex};
        

        responseRange = prctile(theConeExcitations(idx), [5 99]);
        theConeMosaic.renderActivationMap(axHandle, theConeExcitations, ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 responseRange(2)], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', theConeMosaic.integrationTime*1000));
        
        ticksMeters = ticksDegs * mmPerDegree * 1e-3;
        spatialRangeMeters = spatialRangeDegs * mmPerDegree * 1e-3;
        axis 'image'; axis 'ij'; 
        set(gca, 'FontSize', 18, ...
            'XLim', spatialRangeMeters, 'YLim', spatialRangeMeters, ...
            'XTick', ticksMeters, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
            'YTick', ticksMeters, 'YTickLabels', {} ...
            );
        
        if (imIndex == 4)
            xlabel('\it space (degs on retina)');
        else
            xlabel('');
        end
        
        ylabel('');
        if (imIndex==1)
            title('mosaic activation');
        end
        
    end
end


