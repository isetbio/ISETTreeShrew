function plotPsychometricFunction(featureSamples,performance,xLabel,varargin)
% Plot the psychometric function for a given binary search and weibull fit
%
% Syntax:
%   plotPsychometricFunction(featureSamples,performance,'Contrasts (Michelson)','fitResults',psychFitResults)
%
% Inputs:
%   featureSamples - numeric array, levels of the feature that were sampled
%   in the search
%   performance - numeric array, SVM results associated with the feature
%   xLabel

%
% Optional key/value pairs:
%   fitResults - Matlab structure detailing fit of this data with fields
%   'hiResFeature', 'hiResPerformance', 'performanceThreshold' and 'featureThreshold'
%   svmError - error in SVM accuracy- if included, will be displayed as
%   error bars

% See also:

% History:
%   05/14/19 jsc  Wrote initial version.

p = inputParser;
p.addParameter('svmError', 0 ,@isnumeric);
p.addParameter('fitResults',struct,@isstruct);
p.addParameter('ylim',[40,105],@isnumeric);
p.parse(varargin{:});
fitResults = p.Results.fitResults;
svmError = p.Results.svmError;

fit = false;

if ~ismember('fitResults',p.UsingDefaults)
    try
        hiResFeatures = fitResults.hiResFeatures;
        hiResPerformance = fitResults.hiResPerformance;
        featureThreshold = fitResults.featureThreshold;
        performanceThreshold = fitResults.performanceThreshold;
        fit = true;
    catch
        disp('fitResults did not have the proper fields.')
    end
end

figure()

if fit    
    plot(hiResFeatures, hiResPerformance, 'r-', 'LineWidth', 1.5);
    hold on
    
    line([featureThreshold,featureThreshold],[40,performanceThreshold],'Color','black','LineStyle','-')
    line([min(featureSamples),featureThreshold],[performanceThreshold,performanceThreshold],'Color','black','LineStyle','-')
    
    featureSamplesLog = log10(featureSamples);
    featureThresholdLog = log10(featureThreshold);
    % Line seems to plot slightly to the left, so the 1.1 is a temporary fix
    featureThresholdNorm = 1.1*(featureThresholdLog-min(featureSamplesLog))/range(featureSamplesLog);
    performanceThresholdNorm = (performanceThreshold-min(ylim))/range(ylim);
    x = [0.6 featureThresholdNorm];
    y = [0.4 performanceThresholdNorm];
    annotation('textarrow',x,y,'String',sprintf('(%.5f,%.0f)',featureThreshold,performanceThreshold),'FontSize',14,'FontWeight','bold')
    
    set(gca,'xlim',[min(featureSamples),max(featureSamples)], 'XScale', 'log')
    set(gca, 'XTick', unique([featureThreshold, min(featureSamples),max(featureSamples),get(gca, 'XTick')]));
    
    plot(featureSamples,performance, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
else    
    plot(featureSamples,performance, 'ko', 'MarkerSize', 8, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0]);
    hold on
    
    set(gca,'xlim',[min(featureSamples),max(featureSamples)], 'XScale', 'log')
    set(gca, 'XTick', unique([min(featureSamples),max(featureSamples),get(gca, 'XTick')]));
end

if ~ismember('svmError',p.UsingDefaults)
    errorbar(featureSamples,performance,svmError)
end

xt = get(gca, 'XTick');
xtl = cell([1,length(xt)]);
for i = 1:length(xt)
    xtl{i} = num2str(round(xt(i),5));
end
set(gca, 'XTickLabels', xtl);
set(gca, 'YLim', ylim)
set(gca, 'FontSize', 16)
xlabel(xLabel);
ylabel('Performance (% Detected)');
title('Psychometric Function')

hold off
end