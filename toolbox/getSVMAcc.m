function [percentCorrect,trialsPerFold,time] = getSVMAcc(theMosaic, theOI, testScene, nullScene, varargin)
% Calculate the accuracy of an SVM trained to discriminate between a test
% stimulus and a null stimulus using a given cone mosaic and optical
% structure.
%
% Syntax:
%   [percentCorrect,time] = getSVMAcc(theMosaic, testScene, nullScene)
%
% Input:
%    theMosaic              An object of class 'coneMosaic'
%    theOI                  An optical image structure
%    testScene              An scene structure
%    nullScene              Scene to be compared against testScene
%
% Output:
%     percentCorrect        Vector of length 'nFolds' containing accuracy of
%                           each validation, in units of percent, total 
%                           accuracy computed by averaging this vector   
%     trialsPerFold         Number of trials per fold for the accuracy
%     time                  Amount of time taken to run the SVM, in seconds
%
% Optional key/value pairs:
%     'nTrialsNum'          Number of trial iterations to be trained upon 
%                           (default 250)
%      'nFolds'             Number of folds for validation (default 10)
%      'emPathN'             Number of steps of eye movement (default 1)
%      'taskIntervals'      Number of task intervals:
%                           2 (default) - The 'subject' is shown both scenes 
%                           for each trial, and is asked to determine which 
%                           scene is the test scene
%                           1 - The 'subject' is shown both scenes at the 
%                           beginning of the task, and a trial consists of
%                           the subject being shown one scene and asked to
%                           determine which scene it is
%
% See Also:
%   coneMosaic
%   oiCreate
%   sceneCreate
%
% History:
%   02/08/19  jsc drafted
%   04/26/19  jsc formatted 

p = inputParser;
p.addParameter('nTrialsNum', 250, @isnumeric);
p.addParameter('nFolds', 10, @isnumeric);
p.addParameter('taskIntervals', 2, @isnumeric);
p.addParameter('emPathN', 1, @isnumeric);


% Parse input
p.parse(varargin{:});
nFolds = p.Results.nFolds;
nTrialsNum = p.Results.nTrialsNum;
taskIntervals = p.Results.taskIntervals;
emPathN = p.Results.emPathN;

tic;

% Generate various instances of eye movement paths (zero movement of emPathN = 1)
emPath = zeros(nTrialsNum, emPathN, 2);

% Compute the retinal image of the test stimulus
testOI = oiCompute(theOI, testScene);

% Compute the retinal image of the null stimulus
nullOI = oiCompute(theOI, nullScene);

% Compute mosaic excitation responses to the test stimulus
coneExcitationsTest = theMosaic.compute(testOI, 'emPath', emPath);

% Compute mosaic excitation responses to the null stimulus
coneExcitationsNull = theMosaic.compute(nullOI, 'emPath', emPath);

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

% Form the classification matrix.

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
    % In the two interval task, we concatenate [null test] as one class and [test null] as the other.
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
pcVectors = pca(classificationMatrix);

% Project the responses onto the space formed by the first 2 PC vectors
pcComponentsNum = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNum);

% Train a binary SVM classifier
svm = fitcsvm(classificationMatrixProjection,classes);

% Perform a k-fold cross-validation on the trained SVM model
CVSVM = crossval(svm,'KFold',nFolds);

% Compute classification loss for the in-sample responses using a model
% trained on out-of-sample responses
fractionCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');

% Calculate percent correct for each fold
percentCorrect = fractionCorrect.*100;

trialsPerFold = nTrialsNum/nFolds;

time = toc;

end
