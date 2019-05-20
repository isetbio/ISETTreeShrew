function fitResults = getPsychometricFit(feature,performance,nTrials,varargin)
% Plot the psychometric function for a given binary search and weibull fit
%
% Syntax:
%   getPsychometricFit(feature,performance,nTrials,'performanceUnits', 'percent')

% Inputs:
%   feature - x values, the feature of the stimulus that is being varied
%   performance - y values, the performance (accuracy) as a function of the
%   feature
%   nTrials - number of total trials for each feature level
  
%
% Optional key/value pairs:
%   performanceThreshold - Default 0.75. What performance threshold are you
%   interested in reporting the feature for?

% See also: sampledThresholdBinarySearch, plotPsychometricFunction

% History:
%   unknown nc  Wrote initial version.
%   05/14/19 jsc  Formatted for csf.

p = inputParser;

p.addParameter('performanceThreshold', .75, @isnumeric)
p.addParameter('performanceUnits', 'percent', @ischar)
p.parse(varargin{:});

performanceThreshold = p.Results.performanceThreshold;
performanceUnits = p.Results.performanceUnits;

switch performanceUnits
    case 'percent'
        multiplier = 1/100;
    case 'fraction'
        multiplier = 1;
    otherwise
        error("performanceUnits can be 'percent' or 'fraction' of trials.")
end

if performanceThreshold > 1
    performanceThreshold = performanceThreshold/100;
end

% Set up psychometric function model. Here we use a cumulative Weibull function
psychometricFunctionModel = @PAL_Weibull;

fractionCorrect = performance * multiplier;

% Set up search grid
gridLevels = 100;
searchGridParams.alpha = logspace(log10(min(feature)),log10(max(feature)),gridLevels);
searchGridParams.beta = 10.^linspace(-4,4,gridLevels);
searchGridParams.gamma = 0.5;
searchGridParams.lambda = 0.0;

% Optimization settings for the fit
optionsParams             = optimset('fminsearch');
optionsParams.TolFun      = 1e-09;
optionsParams.MaxFunEvals = 1000;
optionsParams.MaxIter     = 1000;
optionsParams.Display     = 'off';

% Parameters for the curve fitting
% Parameters that are allowed to vary
% The parameters are: threshold, slope, guess-rate, lapse-rate
paramsFree = [1 1 0 0];
trialsNumCorrectPerFeatureLevel = round(nTrials*fractionCorrect);
trialsNumPerFeatureLevel = repmat(nTrials,1,length(fractionCorrect));

%% Fit the data and get the best fit params
paramsValues = PAL_PFML_Fit(feature(:), trialsNumCorrectPerFeatureLevel(:), trialsNumPerFeatureLevel(:), ...
    searchGridParams, paramsFree, psychometricFunctionModel, 'SearchOptions', optionsParams);

% Obtain the threshold at which performance cross a threshold performance, here 75%
fitResults.featureThreshold = psychometricFunctionModel(paramsValues, performanceThreshold, 'inverse');

%
% Obtain a high resolution version of the fitted function
fitResults.performanceThreshold = performanceThreshold*100;
fitResults.hiResFeatures = searchGridParams.alpha;
fitResults.hiResPerformance = (1/multiplier)*PAL_Weibull(paramsValues, fitResults.hiResFeatures);
end