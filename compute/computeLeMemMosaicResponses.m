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
    displayRefreshHz = 66;
    
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
    
    % Generate stimulus modulation
    stimulusDurationSeconds = 100/1000;

    % 2 mm pupil
    pupilDiameterMM = 2.0;
    
    % Cone integration time in seconds
    coneIntegrationTime = 2/1000;
    
    % How many mosaic response instances (trials) to compute
    nTrialsNum = 1;

    % Select number of repetitions (more than 1 to estimate compute time)
    nRepeats = 1;
    
    % Define what is to be computed
    doHumanComputation = ~true;
    doTreeShrewComputation = true;
    
    % Define what is to be visualized
    visualizeMosaics = ~true;
    generateVideo = true;
    generateOverviewFigure = true;
    
    
    
    if (doTreeShrewComputation)
        species = 'treeshrew';
        
        % Generate default treeshrew optics
        theOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM);
        
        % Create the  treeshrew cone mosaic
        theConeMosaic = coneMosaicTreeShrewCreate(...
            theOI.optics.micronsPerDegree, ...
            'fovDegs', imageSizeVisualAngle*[1 1], ...
            'sConeMinDistanceFactor', 2, ...
            'integrationTime', coneIntegrationTime);
        
    elseif (doHumanComputation)
        species = 'human';
        
        % Generate default human optics
        theOI = oiCreate('wvf human', pupilDiameterMM);
    
        % Load the human cone mosaic
        load('humanConeMosaic4degs.mat', 'theData')
        theConeMosaic = theData;
        theConeMosaic.integrationTime = coneIntegrationTime;
        theConeMosaic.eccBasedMacularPigment = false;
        clear 'theData';
        
    else
        return;
    end
    
    % Generate stimulus modulation
    [stimulusTimeAxis, stimulusModulationFunction, temporalParams] = ...
        stimulusModulationGenerate(stimulusDurationSeconds,displayRefreshHz);
    
    % Pre-process input RGB images to generate ISETBio scenes
    [inputImageArray, sceneArray] = preProcessImages(imDir, imFileNames,...
        presentationDisplay, ...
        imageSizePixels, ...
        imageSizeVisualAngle, ...
        meanLuminanceCdM2);
    
    % Compute
    [oiArray, coneExcitationArray, demosaicedResponseStructArray, timeAxis] = ...
        runPipeline(species, nRepeats, nTrialsNum, sceneArray, theOI, theConeMosaic, ...
        stimulusTimeAxis, stimulusModulationFunction, temporalParams, imageSizePixels);
    
    % Visualize
    if (visualizeMosaics)
        % Visualize the mosaics in degrees
        displayUnits = 'degs'; figNo = 1;
        mosaicsVisualize(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits);
        % Visualize the mosaics in microns
        displayUnits = 'microns';  figNo = 2;
        mosaicsVisualize(figNo, theTreeShrewConeMosaic, theHumanConeMosaic, displayUnits);
    end
    
    if (generateOverviewFigure)
        figNo = 10;
        hFig = renderOverviewFigure(figNo, inputImageArray, sceneArray, ...
            oiArray, coneExcitationArray, demosaicedResponseStructArray, theConeMosaic, timeAxis);
        NicePlot.exportFigToPDF('TreeShrewMosaicAnalysis.pdf', hFig, 300);
    end

    if (generateVideo)
        figNo = 11;
        for imIndex = 1:numel(sceneArray)
            renderInstancesVideo(figNo, imIndex, oiArray, coneExcitationArray, theConeMosaic, timeAxis);
        end
    end
end


function [oiArray, coneExcitationArray, demosaicedResponseStructArray, timeAxis] = runPipeline(...
    species, nRepeats, nTrialsNum, sceneArray, theOI, theConeMosaic, ...
    stimulusTimeAxis, stimulusModulationFunction, temporalParams, imageSizePixels)
    
    % Setup data containers
    oiArray = cell(1, numel(sceneArray));
    coneExcitationArray = cell(1, numel(sceneArray));
    demosaicedResponseStructArray = cell(1, numel(sceneArray));
    
    tic
        
    for repeat = 1:nRepeats
        fprintf('processing batch %d/%d\n', repeat,nRepeats);
        parfor imIndex = 1:numel(sceneArray)
            fprintf('\t ISETBio processing image %d of %d\n', imIndex, numel(sceneArray));

            % Generate a spatially-uniform scene whose radiance (photon rate)
            % matches that of the stimulus scene
            visualizeScenes = ~true;
            nullScene = zeroContrastSceneMatchingSceneRadianceGenerate(sceneArray{imIndex}, visualizeScenes);

            % Compute the retinal image to the stimulus
            oiArray{imIndex} = oiCompute(theOI, sceneArray{imIndex});

            % Compute the retinal image to the zero contrast stimulus
            oiNull = oiCompute(theOI, nullScene);

            % Compute the optical image sequence
            theOIsequence = oiSequence(oiNull, oiArray{imIndex}, ...
                stimulusTimeAxis, stimulusModulationFunction, 'composition', 'blend');

            % Generate the emPaths, here all zero movements
            emPathLength = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
                theConeMosaic.integrationTime, ...
                'stimulusSamplingInterval', temporalParams.stimulusDurationInSeconds);
            theEMpaths = zeros(nTrialsNum, emPathLength, 2);

            % Compute mosaic cone excitation response
            coneExcitationArray{imIndex} = ...
                theConeMosaic.computeForOISequence(theOIsequence, ...
                'stimulusSamplingInterval', temporalParams.stimulusDurationInSeconds, ...
                'emPaths', theEMpaths);
            timeAxis{imIndex} = theConeMosaic.timeAxis;
            
            % Demosaic response
            % Demosaicing takes a long time, so just do a single frame at t = 0 msec
            demosaicLatencySeconds = 0/1000;
            intepolationMethod = 'linear';
            interpolationMethod = 'nearest';
            interpolationMethod = struct(...
                'name', 'kernel smoothing', ...
                'kernelSigmaMicrons', [3.0 0 7.0] ...
            );
            
            demosaicedResponseStructArray{imIndex} = demosaicResponse(...
                    coneExcitationArray{imIndex}, theConeMosaic, ...
                    demosaicLatencySeconds, imageSizePixels, ...
                    interpolationMethod);

        end % imIndex
    end % for repeat = 1:nRepeats
    
    timeAxis = timeAxis{1};
    
    fprintf('''%s'' pipeline took: %f seconds/image to complete\n', species, toc/numel(sceneArray)/nRepeats);
end