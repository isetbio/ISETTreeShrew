function searchResults = binarySearchOverFeature(feature,featureRange,theMosaic,theOI,varargin)
% Find an estimate for a performance threshold using a binary search
% approach
%
% Syntax:
%   searchResults = binarySearchOverFeature('contrast',[0.001,0.03],theMosaic,theOI)
% Inputs:
%   feature - character vector, the feature of the stimulus that is being varied
%   featureRange - the range of this feature to sample (assume threshold is within)
%   theMosaic - cone mosaic object
%   theOI - optical image structure
%
% Outputs:
%   searchResults: structure with fields:
%       feature - same feature as passed in
%       samples - the levels of the feature that were analyzed
%       performances - the SVM accuracy for each sample
%       finalAccuracy - the final performance (hopefully near the
%       threshold)
%       finalSE - the standard error of the final SVM result
%       thresholdFeature - the level of the feature that led to the
%       finalAccuracy

%
% Optional key/value pairs:
%   stimParams - uses csfStimParamsDefault as default. What stimulus are you
%   showing the SVM? Only one feature will change.
%   presentationDisplay - uses csfStimParamsDefault as default. 
%   desiredAccRange - Default [74.5,75.5] (%). What is the performance threshold?
%   maxCycles - Maximum number of steps before taking best available
%   approximation.
%   nTrialsNum - SVM information: How many trials do you want to train the
%   SVM on? Default is 250.
%   nFolds - SVM information: How many folds do you want to use for k-fold
%   validation? Default is 10.

% See also: t_BinarySearchCSF, plotbinarySearch

% History:
%   04/25/19 jsc  Wrote initial version.

if ~ischar(feature)
    error('The choice of feature should be a character vector.')
end

[defaultStimParams,defaultPresentationDisplay] = csfStimParamsDefault();
nullStimParams = defaultStimParams;
nullStimParams.contrast = 0;

%Parse inputs

p = inputParser;
p.StructExpand = false;
p.addParameter('stimParams', defaultStimParams, @isstruct)
p.addParameter('presentationDisplay', defaultPresentationDisplay, @isstruct)
p.addParameter('desiredAccRange', [74.5 , 75.5],@isnumeric)
p.addParameter('maxCycles', 8, @isnumeric)
p.addParameter('nTrialsNum', 250, @isnumeric)
p.addParameter('nFolds', 10, @isnumeric)

p.parse(varargin{:});
stimParams = p.Results.stimParams;
presentationDisplay = p.Results.presentationDisplay;
maxCycles = p.Results.maxCycles;
desiredAccRange = p.Results.desiredAccRange;
nTrialsNum = p.Results.nTrialsNum;
nFolds = p.Results.nFolds;

if isfield(stimParams,feature)
    stimParams = rmfield(stimParams,feature);
end

% Generate null scene (contrast = 0)

nullScene = generateGaborScene(...
    'stimParams', nullStimParams,...
    'presentationDisplay', presentationDisplay);

maxFeature = max(featureRange);
minFeature = min(featureRange);
samples = NaN(1,maxCycles);
performances = NaN(1,maxCycles);
standardError = NaN(1,maxCycles);
cycle = 0;

while 1
    
    % Generate test stimulus for this sample
    currentFeature = mean([minFeature,maxFeature]);
    
    stimParams.(feature) = currentFeature;
    
    testScene = generateGaborScene(...
        'stimParams', stimParams,...
        'presentationDisplay', presentationDisplay, ...
        'minimumPixelsPerHalfPeriod', 5);
    
    % Use a binary SVM to get a measure of performance (ability to
    % discriminate between null and test stimulus)
    svm_results = getSVMAcc(theMosaic, theOI, testScene, nullScene, 'nTrialsNum', nTrialsNum,'nFolds',nFolds);

    acc = mean(svm_results);
    acc_SE = std(svm_results)/sqrt(length(svm_results));

    % Save variables
    samples(cycle+1) = currentFeature;
    performances(cycle+1) = acc;
    standardError(cycle+1) = acc_SE;
    
    % Determine parameters for next cycle of search
    if acc < min(desiredAccRange)
        minFeature = currentFeature;
    elseif acc > max(desiredAccRange)
        maxFeature = currentFeature;
    else
        % if reached desired range
        thresholdFeature = currentFeature;
        finalAcc = acc;
        finalSEtemp = acc_SE;
        break;
    end
    
    % If you don't make it, find the closest sample- if 75%
    % is within a 99% CI of that SVM, accept that sample
    if cycle > maxCycles
        temp = [abs(75-rmmissing(performances)) ; rmmissing(performances) ; rmmissing(samples) ; rmmissing(standardError)]';
        temp = sortrows(temp,1);
        
        tempAcc = temp(1,2);
        tempSE = temp(1,4);
        
        if (tempAcc-1.6*tempSE) < 75 && (tempAcc+1.6*tempSE) > 75
            finalAcc = tempAcc;
            finalSEtemp = tempSE;
            thresholdFeature = temp(1,3);
            break
        else
            samples = NaN(1,maxCycles);
            performances = NaN(1,maxCycles);
            standardError = NaN(1,maxCycles);
            cycle = -1;
        end
        
    end
    cycle = cycle + 1;
end

% Create structure for results
searchResults.feature = feature;
searchResults.samples = rmmissing(samples);
searchResults.performances = rmmissing(performances);
searchResults.finalAccuracy = finalAcc;
searchResults.finalSE = finalSEtemp;
searchResults.thresholdFeature = thresholdFeature;
end