function [sampledData, success] = resampleThresholdBinarySearch(expInfo,feature,featureSample,featureRange,finalAccuracy,varargin)
% Sample the psychometric threshold space indicated by a binary search
%
% Syntax:
%     thresholdSample = sampledThresholdBinarySearch(expInfo, 
%         'contrast', ...
%         [.02, .003, .01], ...
%         [.001, .03], ...
%         75.4, ...
%         'stepsBeforePlotting',3, ...
%         'acceptedAccRange',[74.5,75.5], ...
%         'minInputPoints',5, ...
%         'numOtherPoints',5, ...
%         'numOutputPoints',5, ...
%         'featureRange',[0,1] ...
%         );
%
% Inputs:
%   expInfo - struct with information about simulated experiment. Contains
%   theMosiac, theOI, stimParams, nTrialNum and nFolds
%   feature - character vector, the feature of the stimulus that is being varied
%   featureSample - numeric array, levels of the feature that were sampled
%   in the search
%   featureRange - the range of this feature to sample (assume threshold is within)
%   finalAccuracy - what was the eventual performance the search reported?
 
%
% Optional key/value pairs:
%   stepsBeforePlotting - steps to ignore in search in order to sample
%   threshold
%   acceptedAccRange - Not all searches reached the same accuracy. What
%   would you consider acceptable?
%   minInputPoints - Number of unique points after the stepsBeforePlotting
%   numOtherPoints - Number of points across the full feature range to
%   also plot
%   numOutputPoints - Number of points to plot in the threshold
 
% See also: t_BinarySearchCSF, t_psychometricFromBinarySearch, plotPsychometricFunction
 
% History:
%   05/15/19 jsc  Wrote initial version.
 
% Parse inputs
p = inputParser;
p.addParameter('stepsBeforePlotting', 3, @isnumeric)
p.addParameter('minInputPoints', 5, @isnumeric)
p.addParameter('numOutputPoints', 5, @isnumeric)
p.addParameter('numOtherPoints', 5, @isnumeric)
p.addParameter('acceptedAccRange', [74.5 , 75.5],@isnumeric)
p.parse(varargin{:});
stepsBeforePlotting = p.Results.stepsBeforePlotting;
minInputPoints = p.Results.minInputPoints;
numOutputPoints = p.Results.numOutputPoints;
numOtherPoints = p.Results.numOtherPoints;
acceptedAccRange = p.Results.acceptedAccRange;
 
 
% Pull information about the simulated experiment to ensure the same
% parameters are used here
 
theMosaic = expInfo.theMosaic;
theOI = expInfo.theOI;
presentationDisplay = expInfo.presentationDisplay;
 
nTrialsNum = expInfo.nTrialsNum;
nFolds = expInfo.nFolds;
 
stimParams = expInfo.stimParams;
if isfield(stimParams,feature)
    stimParams = rmfield(stimParams,feature);
end
 
% Now, we need the results of the binary search
 
success = true;
 
if acceptedAccRange(1) < finalAccuracy && acceptedAccRange(2) > finalAccuracy
    % Great, move on
else
    success = false;
    disp("This search didn't have a final performance within your range.")
end

thresholdFeaturesTemp = unique(round(featureSample(stepsBeforePlotting : length(featureSample)),4));
 
if length(thresholdFeaturesTemp) > minInputPoints
    thresholdFeatures = thresholdFeaturesTemp;
else
    success = false;
    warning("This search didn't have enough unique samples given your constraints.")
end
 
if success
    % Determine the features to analyze
    minC = min(thresholdFeatures)-range(thresholdFeatures)/2;
    maxC = max(thresholdFeatures)+range(thresholdFeatures)/2;
    
    thresholdSamples = minC:(maxC-minC)/numOutputPoints:maxC;
    sampleFullRange = min(featureRange):range(featureRange)/numOtherPoints:max(featureRange);
    
    sampleList = sort([sampleFullRange,thresholdSamples]);
    
    %% Using SVM to calculate accuracy for each feature level
     
    % Generate null scene (contrast = 0)
    nullStimParams = stimParams;
    nullStimParams.contrast = 0;
    
    nullScene = generateGaborScene(...
        'stimParams', nullStimParams,...
        'presentationDisplay', presentationDisplay);
    
    performance = zeros(1,length(sampleList));
    performanceSE = zeros(1,length(sampleList));
    %% Cycle through sampleList, using the SVM to calculate accuracy for each.
    
    T = tic;
    
    count = 1;
    N = length(sampleList);
    
    D = parallel.pool.DataQueue;
    
    w = waitbar(0, [feature,' Space Resampling Progress']);
    afterEach(D, @nUpdateWaitbar);
    
    parfor i = 1:length(sampleList)
        
        % Set the feature for this iteration
        stimParamsTemp = stimParams;  
        stimParamsTemp.(feature) = sampleList(i);
        
        % Generate a Gabor patch with that feature
        testScene = generateGaborScene(...
            'stimParams', stimParamsTemp,...
            'presentationDisplay', presentationDisplay,...
            'minimumPixelsPerHalfPeriod', 5);
        
        % Determine how well the specified optics can discriminate between that
        % patch and a null stimulus
        acc = getSVMAcc(theMosaic, theOI, testScene, nullScene, 'nTrialsNum',nTrialsNum,'nFolds',nFolds);
        performance(i) = mean(acc);
        performanceSE(i) = std(acc)/sqrt(length(acc));
        
        send(D,i)
    end
    
    delete(gcp)
    close(w)
    
    sampledData.time = toc(T)/60;
    sampledData.feature = feature;
    sampledData.stimParams = stimParams;
    sampledData.samples = sampleList;
    sampledData.performance = performance;
    sampledData.performanceSE = performanceSE;
    

else
    % If the search didn't match the requirements, return an empty
    % structure and success = false;
    sampledData = struct;
end
 
%% Functions

function nUpdateWaitbar(~)
waitbar(count/N, w);
count = count + 1;
end

end
