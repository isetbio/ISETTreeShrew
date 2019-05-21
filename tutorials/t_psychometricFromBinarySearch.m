function t_psychometricFromBinarySearch(varargin)
%
%% t_psychometricFromBinarySearch (Human/Treeshrew)
% Use resampleThresholdBinarySearch, getPsychometricFit and
% plotPsychometricFunction to confirm that binary searches over certain
% stimulus features actually estimate the threshold

% Syntax:
%   t_psychometricFromBinarySearch('csfData','load('bestCSFData.mat'))
%   t_psychometricFromBinarySearch('psychDataToPlot',load('psychFn_250_trials.mat'))
%

% Description:
% In Casagrande's 1984 treeshrew CSF paper, the researchers fit their data
% of contrast vs accuracy with a psychometric function, and used it to find
% their desired contrast threshold.
%
% When using an SVM in ISETBio, we don't have a perfect a priori idea of
% threshold location, so we set a wide feature range to search. If we
% evenly sampled the full range, we would need a very large number of data
% points in order to sample the psychometric function. Instead, we use a
% binary search algorithm to find the desired contrast threshold, as can be
% seen in t_BinarySearchCSF and binarySearchOverFeature.
%
% This process assumes that there is an underlying psychometric
% function. Once we've found the stimulus threshold, we can sample
% around that threshold and fit the local data to a psychometric function,
% proving that our binary search does, in fact, approximate the correct
% threshold values.

% Optional key/value pairs:
%   csfData - Binary search data (structure)
%   psychDataToPlot - Previously gathered psychometric function data, if you just
%       want to plot (structure)
%   dataName - If you have a specific choice for the name of the generated
%       data, you can input dataName as a character array.
%   overwrite - If there is already a dataset with the same name, do you
%       want to overwrite it? true/false.
%   saveToWS - If you're just exploring the code, you may want to save the
%       generated data to the workspace after it is created. Logical.

% See also: t_BinarySearchCSF
%

% Tools used: resampleThresholdBinarySearch, getPsychometricFit,
% plotPsychometricFunction
%

% History:
%   04/02/19 jsc  Wrote initial version.

% Required csfData dataset can be created using t_BinarySearchCSF. Data must
% have the variables:
%       expInfo: struct w/
%           theMosaic
%           theOI
%           nTrialsNum
%           nFolds
%           psfSigma
%       binaryResults: struct w/
%           contrasts
%           contrastRange
%           accuracies
%           frequencyRange

% What data do you want to use? This could be overwritten by function

data = load('testCSFData.mat');

% Parse optional input from function
p = inputParser;
p.StructExpand = false;
p.addParameter('csfData', data, @isstruct)
p.addParameter('psychDataToPlot', struct, @isstruct)
p.addParameter('dataName', char.empty, @ischar)
p.addParameter('overwrite', false, @islogical)
p.addParameter('saveToWS', false, @islogical)

p.parse(varargin{:});
data = p.Results.csfData;
psychDataToPlot = p.Results.psychDataToPlot;
dataName = p.Results.dataName;
overwrite = p.Results.overwrite;
saveToWS = p.Results.saveToWS;

% Assume that we are creating a new dataset, unless indicated otherwise by
% passing a previously created psych dataset into the function.
compute = true;

if ~ismember('psychDataToPlot',p.UsingDefaults)
    compute = false;
    if ~isempty(dataName)
        warning('You are loading data instead of generating new data, so no new data will be saved under your chosen name.')
    end
end

if ~ismember('csfData',p.UsingDefaults) && compute == false
    error('Either choose a psychometric function data structure to plot, or create one from binary search data.')
end

expInfo = data.expInfo;
binaryResults = data.binaryResults;

%% Parameters to determine which search to use
%
% In the binary search dataset, there will be multiple searches, for
% different
% spatial frequencies. We want to find a search that fits a few criteria in
% order to ensure we have enough different data points, and that the search
% matches set standards.

%First, we don't want to sample the entire field, and want to take advantage of
%the binary search we've already done. How many steps do we want to take
%before deciding we're close enougn to the threshold?

% 1 : 1/2 of the original contrast range
% 2 : 1/4 of the original contrast range
% ...
stepsBeforePlotting = 3;

% What's the minimum number of unique points you want to sample? So, the
% points after the previous parameter.
minInputPoints = 3;

% What's the range of satisfactory accuracies you want to have reached in the original
% search?
acceptedAccRange = [74.5,75.5];

% How many points within the threshold range do you want to
% compute?
numOutputPoints = 5;

% How many points outside the threshold range would you like to compute?
% These will be evenly spaced, spanning the range of contrasts searched.
numOtherPoints = 5;

if compute
    
    %% Loop through the binary searches in the CSF dataset 
    %
    % Loop until you find one that matches the requirements, then resample
    % using that search
    i = 1;
    while 1
        expInfo.stimParams.spatialFrequencyCyclesPerDeg = binaryResults.frequencyRange(i);
        
        % resampleThresholdBinarySearch returns success = false if the
        % search passed in does not match the requirements passed in,
        % success = true if the requirements were matched and the search
        % area was resampled
        [sampledDataTemp,success] = resampleThresholdBinarySearch(...
            expInfo, ...
            'contrast', ...
            binaryResults.contrasts{1,i}, ...
            binaryResults.contrastRange, ...
            binaryResults.finalAccuracy(i), ...
            'stepsBeforePlotting',stepsBeforePlotting, ...
            'acceptedAccRange',acceptedAccRange, ...
            'minInputPoints',minInputPoints, ...
            'numOtherPoints',numOtherPoints, ...
            'numOutputPoints',numOutputPoints ...
            );
        
        if success
            thresholdSample = sampledDataTemp;
            break
        end
        
        if (i+1) > length(binaryResults.frequencyRange)
            error("None of the searches in the CSF data matched the requirements.")
        end  
        
        i = i + 1;
    end
    
    % Will data be saved to the workspace?
    if saveToWS
        assignin('base','thresholdSample',thresholdSample)
    end
    
    %% Fit the data with a psychometric function
    %
    % Determine number of data points used to compute each accuracy
    % measure
    expInfo = data.expInfo;
    nTrialsNum = expInfo.nTrialsNum;
    nTrialsPerMeasure = nTrialsNum/expInfo.nFolds;
    
    % Fit the data
    psychFitResults = getPsychometricFit( ...
        thresholdSample.samples, ...
        thresholdSample.performance, ...
        nTrialsPerMeasure, ...
        'performanceThreshold', mean(acceptedAccRange));
    
    %% Save the data, including the results of the fit
    %
    % Will data be saved to the workspace?
    if saveToWS
        assignin('base','psychFitResults',psychFitResults)
    end
    
    % Will the data be saved under a generic or chosen name?
    if isempty(dataName)
        psychDataToSave = sprintf('psychFn_%.0f_trials.mat',nTrialsNum);
    elseif endsWith(dataName,'.mat')
        psychDataToSave = dataName;
    else
        psychDataToSave = strcat(dataName,'.mat');
    end
    
    % Will data of the same name be overwritten?
    if ~overwrite
        k = 1;
        psychDataToSaveTemp = psychDataToSave;
        while exist(psychDataToSaveTemp, 'file')
            psychDataToSaveTemp = strcat('v',num2str(k),'_',psychDataToSave);
            k = k + 1;
        end
        psychDataToSave = psychDataToSaveTemp;
    end
    
    save(psychDataToSave,'expInfo','thresholdSample','psychFitResults')
else
    % If you passed in a previously gathered dataset with the proper fields,
    % we'll expand that here.
    thresholdSample = psychDataToPlot.thresholdSample;
    psychFitResults = psychDataToPlot.psychFitResults;
end

%% Plot the psychometric function and the contrast threshold
%
% Either way, plot the data
plotPsychometricFunction(thresholdSample.samples, ...
    thresholdSample.performance, ...
    'Contrasts (Michelson)', ...
    'svmError',thresholdSample.performanceSE, ...
    'fitResults', psychFitResults)
end