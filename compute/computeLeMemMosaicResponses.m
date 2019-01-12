function computeLeMemMosaicResponses
% Compute cone mosaic responses to a set of LeMem images
%
% Syntax:
%   computeLeMemTreeShrewMosaicResponses
%
% Description:
%    Compute cone mosaic responses to a set of LeMem images. Also depict
%    the processing pipeline from RGB images to 256x256 ISETBio scenes to
%    retinal images, and finally to cone mosaic excitation images
%
% Inputs:
%    None
%
% Outputs:
%    None
%

% History:
%    12/12/2018  NPC   Wrote it

    % Create presentation display and place it 5 cm in front of the eye
    presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 5/100);

    % Specify location of LaMem input images
    imDir = 'imagesLaMem';
    
    % Specify names of LaMem input images
    imFileNames = {...
        '00011740' ...
        '00011741' ...
        '00011742' ...
        '00011743' ...
        };

    
    % Specify size of processed image in pixels
    imageSizePixels = [256 256];
    
    % Specify size of processed image in degrees
    imageSizeVisualAngle = 3.95;
    
    % Specify mean luminance of generated scene
    meanLuminanceCdM2 = 40;
    
     % 2 mm pupil
    pupilDiameterMM = 2.0;
    
    % Cone integration time in seconds
    coneIntegrationTime = 2/1000;
    
    % How many mosaic response instances (trials) to compute
    nTrialsNum = 1;

    % Pre-process input RGB images to generate ISETBio scenes
    [inputImageArray, sceneArray] = preProcessImages(imDir, imFileNames,...
        presentationDisplay, ...
        imageSizePixels, ...
        imageSizeVisualAngle, ...
        meanLuminanceCdM2);

    % Generate default treeshrew optics
    theTreeShrewOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM);
    
    % Generate default human optics
    theHumanOI = oiCreate('wvf human', pupilDiameterMM);
    
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

    
    % Setup data containers
    oiArrayTreeShrew = cell(1, numel(sceneArray));
    coneExcitationArrayTreeShrew = cell(1, numel(sceneArray));
    oiArrayHuman = cell(1, numel(sceneArray));
    coneExcitationArrayHuman = cell(1, numel(sceneArray));
    
    % Define what is to be computed
    doHumanComputation = ~true;
    doTreeShrewComputation = true;
    visualizeMosaics = ~true;
    visualizeHumanResults = false;
    visualizeTreeShrewResults = true;
    demosaicResponses = ~true;
    demosaicingSampleSpacingMicrons = 1;
    
    if (visualizeMosaics)
        % Visualize the mosaics in degrees
        displayUnits = 'degs'; figNo = 1;
        visualizeMosaics(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits);
        % Visualize the mosaics in microns
        displayUnits = 'microns';  figNo = 2;
        visualizeMosaics(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits);
    end
    
    % Stimulus timing params - 100 msec stimulus
    displayRefreshRate = 66;
    stimulusSamplingIntervalInSeconds = 1.0/displayRefreshRate;
    stimulusDurationInSeconds = round(0.1/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    secondsForResponseStabilization = round(0.05/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    secondsForResponseExtinction = round(0.05/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    
    temporalParams = struct(...
        'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
        'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
        'secondsForResponseStabilization', secondsForResponseStabilization, ...
        'secondsForResponseExtinction', secondsForResponseExtinction, ...
        'addCRTrasterEffect', ~true, ...
        'rasterSamples', 5 ...
    );
    [stimulusTimeAxis, stimulusModulationFunction, rasterModulation] = squareTemporalWindowCreate(temporalParams);
    if (~isempty(rasterModulation))
        stimulusModulationFunction = stimulusModulationFunction .* rasterModulation;
    end
    
    % Select number of repetitions (more than 1 to estimate compute time)
    nRepeats = 1;

    % Action !
    if (doHumanComputation)
    tic
    
    for repeat = 1:nRepeats
        fprintf('processing batch %d/%d\n', repeat,nRepeats);
        parfor imIndex = 1:numel(sceneArray)
            fprintf('\t processing image %d of %d\n', imIndex, numel(sceneArray));
            % Compute the retinal image
            oiArrayHuman{imIndex} = oiCompute(theHumanOI, sceneArray{imIndex});
            % Compute mosaic excitation
            coneExcitationArrayHuman{imIndex} = theHumanConeMosaic.compute(oiArrayHuman{imIndex}, 'emPath', emPath);
          
            if (demosaicResponses)
                if (emPathLength > 1)
                    % Compute mean response over trials
                    meanConeExcitationCount = squeeze(mean(coneExcitationArrayHuman{imIndex},1));
                     % Only use response at first time bin
                    meanConeExcitationCount = squeeze(meanConeExcitationCount(:,:,1));
                else
                    % Compute mean response over trials
                    meanConeExcitationCount = mean(coneExcitationArrayTreeShrew{imIndex},1);
                end
                
                coneExcitationsRates = meanConeExcitationCount/theHumanConeMosaic.integrationTime;
                
                % Compute de-mosaiced L-cone responses
                [theDemosaicedLconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, lConeXlocsDegs, lConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('L-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                % Compute de-mosaiced M-cone responses
                [theDemosaicedMconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theMConeIsomerizations, MConeXlocsDegs, MConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('M-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                % Compute de-mosaiced S-cone responses
                [theDemosaicedSconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, sConeXlocsDegs, sConeYlocsDegs] = ...
                    theHumanConeMosaic.demosaicConeTypeActivationFromFullActivation('S-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);
            end
            
        end % inIndex
    end % for repeat = 1:nRepeats
        
    fprintf('Human computation: %f\n', toc/4/nRepeats);
    if (visualizeHumanResults)
        fprintf('Rendering human figure. This will take a long time\n');
        renderFigure(20, inputImageArray, sceneArray, oiArrayHuman, coneExcitationArrayHuman, theHumanConeMosaic, nTrialsNum, emPathLength);
    end
    end % if (doHumanComputation)
    
    
    % Action !
    if (doTreeShrewComputation)
    tic
    
            
    for repeat = 1:nRepeats
        fprintf('processing batch %d/%d\n', repeat,nRepeats);
        for imIndex = 1:numel(sceneArray)
            fprintf('\t processing image %d of %d\n', imIndex, numel(sceneArray));
            
            % Compute the retinal image to the stimulus
            oiArrayTreeShrew{imIndex} = oiCompute(theTreeShrewOI, sceneArray{imIndex});
            
            
            photons = sceneGet(sceneArray{imIndex}, 'photons');
            figure(22);
            subplot(1,2,1)
            imagesc(sceneGet(sceneArray{imIndex}, 'rgb image'));
            axis image;
            photons = bsxfun(@plus, photons*0, mean(mean(photons,1),2));
            sceneArrayZeroContrast = sceneArray{imIndex};
            sceneArrayZeroContrast = sceneSet(sceneArrayZeroContrast, 'photons', photons);
            photons = sceneGet(sceneArrayZeroContrast, 'photons');
            subplot(1,2,2)
            imagesc(sceneGet(sceneArrayZeroContrast, 'rgb image'));
            axis image;
        
            % Compute the retinal image to the stimulus
            oiArrayTreeShrewNull = oiCompute(theTreeShrewOI, sceneArrayZeroContrast);
    
            % Compute the optical image sequence
            theOIsequence = oiSequence(oiArrayTreeShrewNull, oiArrayTreeShrew{imIndex}, ...
                stimulusTimeAxis, stimulusModulationFunction, 'composition', 'blend');
%             pause
%             theOIsequence.visualize('movie illuminance');
%             pause

            % Compute the emPaths
            emPathLength = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theTreeShrewConeMosaic.integrationTime, ...
                'stimulusSamplingInterval', temporalParams.stimulusDurationInSeconds);
            theEMpaths = zeros(nTrialsNum, emPathLength, 2);
            
            % Compute mosaic excitation - single shot
            %coneExcitationArrayTreeShrew{imIndex} = theTreeShrewConeMosaic.compute(oiArrayTreeShrew{imIndex}, 'emPath', emPath);
            
            % Compute mosaic excitation - oi sequence
            coneExcitationArrayTreeShrew{imIndex} = theTreeShrewConeMosaic.computeForOISequence(theOIsequence, ...
                'stimulusSamplingInterval', temporalParams.stimulusDurationInSeconds, ...
                'emPaths', theEMpaths);

            if (demosaicResponses) 
                
                % This is for single shot compute, not for oiSequence
                % compute - fix it
                if (emPathLength > 1)
                    % Compute mean response over trials
                    meanConeExcitationCount = squeeze(mean(coneExcitationArrayTreeShrew{imIndex},1));
                     % Only use response at first time bin
                    meanConeExcitationCount = squeeze(meanConeExcitationCount(:,:,1));
                else
                    % Compute mean response over trials
                    meanConeExcitationCount = mean(coneExcitationArrayTreeShrew{imIndex},1);
                end
                
                coneExcitationsRates = meanConeExcitationCount/theTreeShrewConeMosaic.integrationTime;
                
                % Compute de-mosaiced L-cone responses
                [theDemosaicedLconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, lConeXlocsDegs, lConeYlocsDegs] = ...
                    theTreeShrewConeMosaic.demosaicConeTypeActivationFromFullActivation('L-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);

                % Compute de-mosaiced S-cone responses
                [theDemosaicedSconeIsomerizationMap, demosaicedMapSupportDegs, ...
                 theLConeIsomerizations, sConeXlocsDegs, sConeYlocsDegs] = ...
                    theTreeShrewConeMosaic.demosaicConeTypeActivationFromFullActivation('S-cones', ...
                        coneExcitationsRates, demosaicingSampleSpacingMicrons);
            end
            
        end % inIndex
    end % for repeat = 1:nRepeats
    fprintf('TreeShrew computation: %f\n', toc/4/nRepeats);
    
    if (visualizeTreeShrewResults)
        renderOverviewFigure(10, inputImageArray, sceneArray, oiArrayTreeShrew, coneExcitationArrayTreeShrew, theTreeShrewConeMosaic, nTrialsNum, emPathLength);
        for imIndex = 1:numel(sceneArray)
            renderInstancesVideo(11, imIndex,oiArrayTreeShrew, coneExcitationArrayTreeShrew, theTreeShrewConeMosaic, nTrialsNum, emPathLength);
        end
    end % visualizeTreeShrewResults
    end % if (doTreeShrewComputation)

end

function visualizeMosaics(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits)

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


% Method to render a video showing the response instances to all optical images
function renderInstancesVideo(figNo, imIndex, oiArray, coneExcitationArray, theConeMosaic, nTrialsNum, emPathLength)
    
    videoFileName = sprintf('LeMem%d_ResponseInstances', imIndex);
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
        
    hFig = figure(figNo); 
    clf;
    set(hFig, 'Position', [10 10 1450 700], 'Color', [1 1 1]);
    
    subplotRows = 1;
    subplotCols = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.08, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
    
    oiRGB = oiGet(oiArray{imIndex}, 'rgb');
    spatialSupportMM = oiGet(oiArray{imIndex}, 'spatial support', 'mm');
    optics = oiGet(oiArray{imIndex}, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportXDegs = spatialSupportDegs(1,:,1);
    spatialSupportYDegs = spatialSupportDegs(:,1,2);
    spatialRangeDegs = max(theConeMosaic.fov)*0.5*[-1 1];
    ticksDegs = [-5:1:5];
    
    subplot('Position', subplotPosVectors(1,1).v);
    image(spatialSupportXDegs,spatialSupportYDegs, oiRGB);
    axis 'image'; axis 'ij'; grid on
    set(gca, 'FontSize', 18, ...
        'XLim', spatialRangeDegs, 'YLim', spatialRangeDegs, ...
        'XTick', ticksDegs, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
        'YTick', ticksDegs, 'YTickLabels', {} ...
        );
    xlabel('\it space (degs on retina)');
    title('retinal image');

        
    axHandle = subplot('Position', subplotPosVectors(1,2).v);
    theConeExcitationsAllInstances = coneExcitationArray{imIndex};
    
    
    responseRange = prctile(theConeExcitationsAllInstances(:), [5 95]);
    
    iTrial = 1;
    tmp = theConeMosaic.reshapeHex2DmapToHex3Dmap(squeeze(theConeExcitationsAllInstances(iTrial,:,:)));
    
    for tBin = 1:emPathLength
        theConeMosaic.renderActivationMap(axHandle, squeeze(tmp(:,:,tBin)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 responseRange(2)], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', theConeMosaic.integrationTime*1000));

        ticksMeters = ticksDegs * mmPerDegree * 1e-3;
        axis 'image'; axis 'ij'; 
        set(gca, 'FontSize', 18, 'YTickLabel', {}, 'XTick', ticksMeters, 'XTickLabels', ticksDegs);

        xlabel('\it space (degs on retina)');
        ylabel('');
        time = theConeMosaic.integrationTime*1000*tBin;
        title(sprintf('mosaic activation (%d ms)', time));
        drawnow;
        
        % Add video frame
        videoOBJ.writeVideo(getframe(hFig));
    end % tBin
    
    % Close video stream
    videoOBJ.close();
    fprintf('File saved in %s\n', videoFileName);
        
end


% Method to render all the processing stages
function renderOverviewFigure(figNo, inputImageArray, sceneArray, oiArray, coneExcitationArray, theConeMosaic, nTrialsNum, emPathLength)
    
    hFig = figure(figNo); 
    clf;
    set(hFig, 'Position', [10 10 1450 1500], 'Color', [1 1 1]);
    
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
   
    nonNullConeIndices = find(theConeMosaic.pattern > 1);
    
    ticksDegs = [-10:1:10];
    spatialRangeDegs = max(theConeMosaic.fov)/2*[-1 1];
    
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
        theConeExcitationsAllInstances = coneExcitationArray{imIndex};

        responseRange = prctile(theConeExcitationsAllInstances(:), [5 95]);
        
        iTrial = 1;
        tmp = theConeMosaic.reshapeHex2DmapToHex3Dmap(squeeze(theConeExcitationsAllInstances(iTrial,:,:)));
        
        timeAxis = theConeMosaic.timeAxis;
        timeAxis = timeAxis - mean(timeAxis);
        [~, centerTimeBin] = min(abs(timeAxis));
        meanResponse = squeeze(tmp(:,:,centerTimeBin));
        

        theConeMosaic.renderActivationMap(axHandle, meanResponse, ...
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




function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreate(temporalParams)

    stimulusSamples = round(temporalParams.stimulusDurationInSeconds/temporalParams.stimulusSamplingIntervalInSeconds);
    if (isfield(temporalParams, 'secondsForResponseStabilization'))
        stabilizingTimeSamples = round(temporalParams.secondsForResponseStabilization/temporalParams.stimulusSamplingIntervalInSeconds);
    else
        stabilizingTimeSamples = 0;
    end

    if (isfield(temporalParams, 'secondsForResponseExtinction'))
        extinctionTimeSamples = round(temporalParams.secondsForResponseExtinction/temporalParams.stimulusSamplingIntervalInSeconds);
    else
        extinctionTimeSamples = 0;
    end

    sampleTimes = 1:(stabilizingTimeSamples + stimulusSamples + extinctionTimeSamples);
    sampleTimes = sampleTimes - stabilizingTimeSamples - round(stimulusSamples/2);
    sampleTimes = sampleTimes * temporalParams.stimulusSamplingIntervalInSeconds;

    squareTemporalWindow = zeros(1,numel(sampleTimes));
    onTime = find(sampleTimes >= -temporalParams.stimulusDurationInSeconds/2);
    squareTemporalWindow(onTime(1)+(0:(stimulusSamples-1))) = 1;

    if (isfield(temporalParams, 'addCRTrasterEffect')) && (temporalParams.addCRTrasterEffect)
        % Add CRT raster effect
        phosphorFunction = crtPhosphorActivationFunction(1/temporalParams.stimulusSamplingIntervalInSeconds, temporalParams.rasterSamples);
        rasterSamples = numel(phosphorFunction.timeInSeconds);
        raster = zeros(1,numel(squareTemporalWindow)*rasterSamples);
        raster(1,1:rasterSamples:end) = squareTemporalWindow*0+1;
        rasterModulation = conv(raster, phosphorFunction.activation);
        rasterModulation = rasterModulation(1:numel(squareTemporalWindow)*rasterSamples);

        tmp = zeros(1,numel(squareTemporalWindow)*rasterSamples);
        for i = 1:numel(squareTemporalWindow)*rasterSamples-1
            tmp(i) = squareTemporalWindow(floor((i-1)/rasterSamples)+1);
        end
        squareTemporalWindow = tmp;
        sampleTimes = linspace(sampleTimes(1), sampleTimes(end), numel(squareTemporalWindow));
    else
        rasterModulation = [];
    end

    figure()
    plot(sampleTimes, squareTemporalWindow, 'ks-'); hold on;
    if (~isempty(rasterModulation ))
        plot(sampleTimes, rasterModulation, 'rs-');
    end

end

function phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
% phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
%
% Create a phoshor activarion function with sharp rise, shower decline
%
%  7/7/16  npc Wrote it.
%

    alpha = 1.9; t_50 = 0.02/1000; n = 2;
    phosphorFunction.timeInSeconds = linspace(0,1.0/refreshRate, samplesPerRefreshCycle);
    phosphorFunction.activation = (phosphorFunction.timeInSeconds.^n)./(phosphorFunction.timeInSeconds.^(alpha*n) + t_50^(alpha*n));
    phosphorFunction.activation = phosphorFunction.activation - phosphorFunction.activation(end);
    phosphorFunction.activation(phosphorFunction.activation<0) = 0;
    phosphorFunction.activation = phosphorFunction.activation / max(phosphorFunction.activation);
    figure(3);
    plot(phosphorFunction.timeInSeconds, phosphorFunction.activation, 'ko-');
    
end