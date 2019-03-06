%% t_TestSVM
% Illustrate basic principles of how utilizing SVM for a 2AFC utilizing tree shrew optics...
%
% Description:
%    ...
%
% See also:
% ls_inferenceBinarySVM.mlx

% History:
%   02/17/19 jsc Wrote initial version.

%% Define Parameters
nTrialsNum = 50;

emPathLength = 1;

taskIntervals = 2;

theMosaic = coneMosaicTreeShrewCreate(75, ...%theOI.optics.micronsPerDegree, ...
    'fovDegs', 14, ...        % match mosaic width to stimulus size
    'integrationTimeSeconds', 10/1000);

emPath = zeros(nTrialsNum, emPathLength, 2);

%% Create presentation display and place it 5 cm in front of the eye
presentationDisplay = displayCreate('LCD-Apple', 'viewing distance', 5/100);

%% begin loop
v = 1:1:3;
acc = zeros(1,length(v));
time = zeros(1,length(v));
contrast = .05;
c = 1;
for i = v
    
    tic
    
    %% parameter struct for a low spatial frequency Gabor stimulus
    stimParams = struct(...
        'spatialFrequencyCyclesPerDeg', 3, ... % changing cycles/deg
        'orientationDegs', 0, ...               % 0 degrees
        'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
        'sizeDegs', 14, ...                     % 14 x 14 size
        'sigmaDegs', 100, ...                   % sigma of Gaussian envelope
        'contrast', contrast,...                   % 0.05 Michelson contrast
        'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
        'pixelsAlongWidthDim', [], ...          % pixels- width dimension
        'pixelsAlongHeightDim', [] ...          % pixel- height dimension
        );
    
    %%
    
    % Generate a scene representing the Gabor stimulus with the above params as
    % realized on the presentationDisplay
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);

    %%
    % zero contrast for the null stimulus
    stimParams.contrast = 0.0;
    
    % Generate a scene representing the 10% Gabor stimulus as realized on the presentationDisplay

    nullScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay);

    %% Generate wavefront-aberration derived ts optics
    theOI = oiTreeShrewCreate();
    
    %% Compute the retinal image of the test stimulus
    theTestOI = oiCompute(theOI, testScene);
    
    %% Compute the retinal image of the null stimulus
    theNullOI = oiCompute(theOI, nullScene);
    
    %%
    %
    % Compute mosaic excitation responses to the test stimulus
    coneExcitationsTest = theMosaic.compute(theTestOI, 'emPath', emPath);
    %%
    % Compute mosaic excitation responses to the null stimulus
    coneExcitationsNull = theMosaic.compute(theNullOI, 'emPath', emPath);
    
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
    
    %% Form the classification matrix.
    
    %
    if (taskIntervals == 1)
        % In the one interval task, the null and test response instances are labelled as the 2 classes.
        % Allocate matrices
        classificationMatrix = nan(2*nTrials, responseSize);
        classes = nan(2*nTrials, 1);
        % Class 1
        classificationMatrix(1:nTrials,:) = nullResponses;
        classes((1:nTrials)) = 0;
        % Class 2
        classificationMatrix(nTrials+(1:nTrials),:) = testResponses;
        classes(nTrials+(1:nTrials)) = 1;
    elseif (taskIntervals == 2)
        % In the two inteval task, we concatenate [null test] as one class and [test null] as the other.
        % Allocate matrices
        classificationMatrix = nan(nTrials, 2*responseSize);
        classes = nan(nTrials, 1);
        halfTrials = floor(nTrials/2);
        % Class 1
        classificationMatrix(1:halfTrials,:) = [...
            nullResponses(1:halfTrials,:) ...
            testResponses(1:halfTrials,:)];
        classes((1:halfTrials)) = 0;
        % Class 2
        idx = halfTrials+(1:halfTrials);
        classificationMatrix(idx,:) = [...
            testResponses(idx,:) ...
            nullResponses(idx,:)];
        classes(idx) = 1;
    else
        error('Task can have 1 or 2 intervals only.')
    end
    
    % Find principal components of the responses
    [pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);
    
    % Project the responses onto the space formed by the first 2 PC vectors
    pcComponentsNum = 2;
    classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNum);
    
    %% Train a binary SVM classifier
    
    svm = fitcsvm(classificationMatrixProjection,classes);
    
    % Perform a 10-fold cross-validation on the trained SVM model
    kFold = 10;
    CVSVM = crossval(svm,'KFold',kFold);
    
    % Compute classification loss for the in-sample responses using a model
    % trained on out-of-sample responses
    %%
    fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
    %%
    % Average percent correct across all folds
    
    percentCorrect = mean(fractionCorrect) * 100;
    acc(c) = percentCorrect;
    time(c) = toc;
    
    responses = coneExcitationsTest;
    responseRange = [min(responses(:)) max(responses(:))];
    
    c=c+1
end
%%
figure()
semilogx(v,acc,'o')
%xlim([0,max(v)])
ylim([50,100])
xlabel('Spatial Frequency (Cyc/Deg)')
ylabel('SVM Accuracy (%Correct)')
title({'SVM Tree Shrew CSF',sprintf('%.2f Michelson Contrast, N=%.0f',contrast,nTrialsNum)})

sum(time)