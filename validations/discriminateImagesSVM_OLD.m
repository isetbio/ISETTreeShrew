% In this script, we use the most up-to-date tree shrew model to get the
% optical image and generate the cone responses to Kell images for our
% match-to-sample task. We then use an SVM to test whether it can identify
% targets vs distractors. This script is outdated because it uses the old
% coneMosaicHex rather than the updated, faster cMosaic.

% History:
%   08/05/2022  eem  Created based on treeShrewSceneDisplay and
%                    ls_inferenceBinarySVM live scripts

% Configuration.  Execute the tbUseProjec command
% below to configure, if you use ToolboxToolbox.
%{
    tbUseProject('ISETTreeShrew');
%}

%% Initialize
clear; close all;

%% Get optical image from ISETTreeShrew
oi = oiTreeShrewCreate('opticsType', 'wvf', ...
    'name', 'wvf-based optics');

%% Create presentation display
% Set number of trials and images to use
nTrialsNum = 100;
imgidx = 1:2; %1:length(targetFile)

% Define any set of targets and distractors with imgidx to display  
% generated scene, optical image, and cone responses. Current path is for
% match-to-sample camels and wrenches
folderName = 'C:/Users/eemeyer/Documents/GitHub/arcaro_rotation/imagematch_test/';
imgFolder{1} = dir([folderName,'camels/target*']);
imgFolder{2} = dir([folderName,'wrenches/distractor*']);

% Set mean luminance and FOV
meanLuminanceCdPerM2 = 20;
sceneFOVdegs = 5.0;

% Generate default treeshrew cone mosaic with a 5 msec integration time.
% Set microns/deg - might need to change this?
umPerDegree = 76;
imageProps = [1 1];
integrationTime = 1/1000;
sizeConeExcit = 564;

% Set which figures should be plotted for the scene
displayContrastProfiles = false;
displayRadianceMaps = false;
displayRetinalContrastProfiles = false;
displayConeMosaic = false;
displayImages = false;

% Generate object using display measurements from match-to-sample
% experiment.
presentationDisplay = generateTreeShrewDisplayObject_eemParams();

% Initialize figures for generated scene and optical images
if displayImages
    figure;
    genScene = tiledlayout('flow','TileSpacing','Compact');
    figure;
    opticalImage = tiledlayout('flow','TileSpacing','Compact');
end

% Initialize coneExcitations. Not sure where the size of coneExcitations is
% determined to initilize the size, but sizeConeExcit can be replaced.
coneExcitations = {};
coneExcitations{1} = zeros(length(imgidx),nTrialsNum,sizeConeExcit,sizeConeExcit);
coneExcitations{2} = zeros(length(imgidx),nTrialsNum,sizeConeExcit,sizeConeExcit);

% create cone mosaic structure
theConeMosaic = coneMosaicTreeShrewCreate(...
    umPerDegree, ...
    'fovDegs', sceneFOVdegs*imageProps, ... % double check these values are correct
    'integrationTime', integrationTime);

% Loop through the number of images and loop through the targets and
% distractors (itype). Resulting cone excitations for all images are stored
% in the cell array coneExcitations with cell 1: targets, cell 2:
% distractors
for img = 1:length(imgidx)
    for itype = 1:2
        % Load the camels (targets) and realized it on presentation display. 
        imgFile = fullfile(imgFolder{itype}(img).folder,imgFolder{itype}(img).name);
        realizedStimulusScene = sceneFromFile(imgFile, 'monochrome', ...
            meanLuminanceCdPerM2, presentationDisplay);

        % Set FOV in scene structure
        realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', sceneFOVdegs);

        % Compute the retinal image
        oi = oiCompute(oi, realizedStimulusScene);
        
        if displayImages
        % Visualize generated scene and computed optical image
        % Visualize different aspects of the 
        ax_temp = nexttile(genScene);
        visualizeScene(realizedStimulusScene, 'displayContrastProfiles', displayContrastProfiles,...
            'displayRadianceMaps',displayRadianceMaps,'axesHandle',ax_temp);

        % Visualize different aspects of the computed optical image
        % nexttile(opticalImage);
        ax_temp = nexttile(opticalImage);
        visualizeOpticalImage(oi, 'displayRetinalContrastProfiles', displayRetinalContrastProfiles,...
            'displayRadianceMaps',displayRadianceMaps,'axesHandle',ax_temp);
        end

        % Generate cone responses and visualize
        emPath = zeros(nTrialsNum, 1, 2);
        % Compute mosaic excitation responses
        coneExcitations{itype}(img,:,:,:) = theConeMosaic.compute(oi, 'emPath', emPath);

        % Visualize cone mosaic and its cone excitation responses
        if displayConeMosaic
            theConeMosaic.displayInfo();
            theConeMosaic.visualizeGrid('backgroundColor', [1 1 1], 'generateNewFigure', true); 

            visualizeConeMosaicResponses(theConeMosaic, squeeze(coneExcitations{itype}(img,:,:,:)), 'R*/cone/tau');
        end
    end
end

%% Attempt at SVM, code from ls_inferenceBinarySVM in ISETLiveScripts
classLoss = zeros(length(imgidx),length(imgidx));
for tt = 1:length(imgidx)
    for dd = 1:length(imgidx)
        % Simulate a 2AFC task 
        taskIntervals = 1;
        [classificationMatrix, classLabels] = generateSetUpForClassifier(theConeMosaic, ...
            squeeze(coneExcitations{1}(tt,:,:,:)), squeeze(coneExcitations{2}(dd,:,:,:)), taskIntervals);

        % Find principal components of the responses
        [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

        % Project the responses onto the space formed by the first 4 PC vectors
        pcComponentsNumForClassification = 2;
        classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

        % Visualize the first 4 principal components.
        % visualizePrincipalComponents(pcVectors, varianceExplained, theConeMosaic);

        % Train a binary SVM classifier and visualize the support vectors in 2
        % dimensions
        svm = fitcsvm(classificationMatrixProjection,classLabels);
        
        % Visualize the data along with the hyperplane computed by the SVM 
        visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);

        % Cross validation, default is 10-fold
        CVSVMModel = crossval(svm);

        % Estimate out-of-sample misclassification rate
        classLoss(tt,dd) = kfoldLoss(CVSVMModel);
    end
end

figure;
imagesc(classLoss)

%% Necessary functions for SVM classification and visualization
function [classificationMatrix, classLabels] = generateSetUpForClassifier(theMosaic, ...
    coneExcitationsTest, coneExcitationsNull, taskIntervals)

% Obtain the indices of the grid nodes that contain cones
[~,~,~, nonNullConeIndices] = theMosaic.indicesForCones;

% Extract the response vectors for nodes containing cones
[nTrials, nRows, mCols, nTimeBins] = size(coneExcitationsTest);
coneExcitationsTestReshaped = reshape(coneExcitationsTest, [nTrials nRows*mCols nTimeBins]);
coneExcitationsNullReshaped = reshape(coneExcitationsNull, [nTrials nRows*mCols nTimeBins]);
testResponses = coneExcitationsTestReshaped(:, nonNullConeIndices, :);
nullResponses = coneExcitationsNullReshaped(:, nonNullConeIndices, :);

% Collapse response vectors across space and time
responseSize = numel(nonNullConeIndices)*nTimeBins;
testResponses = reshape(testResponses, [nTrials responseSize]);
nullResponses = reshape(nullResponses, [nTrials responseSize]);
    
% Assemble the response vectors into a classification matrix simulating either 
% a one interval task or a two-interval task.
if (taskIntervals == 1)
    % In the one interval task, the null and test response instances are labelled as the 2 classes.
    % Allocate matrices
    classificationMatrix = nan(2*nTrials, responseSize);
    classLabels = nan(2*nTrials, 1);
    % Class 1
    classificationMatrix(1:nTrials,:) = nullResponses;
    classLabels((1:nTrials)) = 0;
    % Class 2
    classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
    classLabels(nTrials+(1:nTrials)) = 1;
elseif (taskIntervals == 2)
    % In the two inteval task, we concatenate [null test] as one class and [test null] as the other. 
    % Allocate matrices
    classificationMatrix = nan(nTrials, 2*responseSize);
    classLabels = nan(nTrials, 1);
    halfTrials = floor(nTrials/2);
    % Class 1
    classificationMatrix(1:halfTrials,:) = [...
        nullResponses(1:halfTrials,:) ...
        testResponses(1:halfTrials,:)];
    classLabels((1:halfTrials)) = 0;
    % Class 2
    idx = halfTrials+(1:halfTrials);
    classificationMatrix(idx,:) = [...
        testResponses(idx,:) ...
        nullResponses(idx,:)];
    classLabels(idx) = 1;
else
    error('Task can have 1 or 2 intervals only.')
end
end

function visualizePrincipalComponents(pcVectors, varianceExplained, theMosaic)
    figure(); clf;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 2, ...
       'heightMargin',  0.1, ...
       'widthMargin',    0.01, ...
       'leftMargin',     0.01, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.1, ...
       'topMargin',      0.1);
    for pcaComponentIndex = 1:4
        title = sprintf('PCA %d, variance\nexplained: %2.2f%%', ...
            pcaComponentIndex, varianceExplained(pcaComponentIndex));
        r = floor((pcaComponentIndex-1)/2)+1;
        c = mod(pcaComponentIndex-1,2)+1;
        ax = subplot('Position', subplotPosVectors(r,c).v);
        theMosaic.renderActivationMap(ax, pcVectors(:,pcaComponentIndex), ...
             'fontSize', 14, ...
              'titleForMap', title);
        ylabel('')
        set(ax, 'YTick', [])
        if (pcaComponentIndex < 3)
            xlabel('')
            set(ax, 'XTick', [])
        end
    end
end

function visualizeSVMmodel(svmModel, data, classes)
    sv = svmModel.SupportVectors;
    h = max(abs(data(:)))/100; % Mesh grid step size
    r = -h*100:h:h*100;
    [X1,X2] = ndgrid(r, r);
    [~,score] = predict(svmModel,[X1(:),X2(:)]);
    scoreGrid = reshape(score(:,1),numel(r), numel(r));

    figure(); clf;
    class0Indices = find(classes == 0);
    class1Indices = find(classes == 1);
    
    plot(data(class0Indices,1),data(class0Indices,2),'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c'); hold on
    plot(data(class1Indices,1), data(class1Indices,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    contourf(X1,X2,scoreGrid, 50); 
    plot(data(class0Indices,1),data(class0Indices,2),'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'c')
    plot(data(class1Indices,1), data(class1Indices,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(sv(:,1),sv(:,2),'ks','MarkerSize',12, 'LineWidth', 1.5);
    colormap(brewermap(1024, 'RdBu'))
    hold off
    xlabel('PC component #1 activation')
    ylabel('PC component #2 activation')
    legend('null stimulus', 'test stimulus')
    set(gca, 'FontSize',14)
    axis 'square'
end