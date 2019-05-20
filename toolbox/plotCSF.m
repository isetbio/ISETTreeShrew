function plotCSF(binaryResults,varargin)
% Plot the contrast sensitivity function for a given simulated visual
% system front-end.
%
% Syntax:
%   plotCSF(data,'toDivide','20','plotCasagrandeData','TRUE')
%
% Inputs:
%   binaryResults - Matlab structure with fields 'thresholdContrasts' and
%   'frequencyRange'
%
% Optional key/value pairs:
%   expInfo - Matlab structure with fields 'nTrialsNum', 'psfSigma' and
%   'theMosaic', to include in title 
%   toDivide - Factor to shift sensitivity function downwards, if desired
%   plotCasagrandeData - Logical. Plot Casagrande's experimental data?

% See also:
%   t_BinarySearchCSF

% History:
%   04/08/19 jsc  Wrote initial version.
%   04/28/19 jsc  Formatted.

p = inputParser;
p.addParameter('expInfo',struct,@isstruct);
p.addParameter('toDivide', 1 , @isnumeric);
p.addParameter('plotCasagrandeData', false , @islogical);

% Parse input
p.parse(varargin{:});
expInfo = p.Results.expInfo;
toDivide = p.Results.toDivide;
plotCasagrandeData = p.Results.plotCasagrandeData;

% Shifting the sensitivity downwards by a factor of 'toDivide', if desired
sensitivity = (1./binaryResults.thresholdContrasts)./toDivide;
frequencyRange = binaryResults.frequencyRange;

% If details about the simulated experiment are passed in, use them to
% inform the title
if ismember('expInfo',p.UsingDefaults)
    plotTitle = 'Contrast Sensitivity Function';  
else
    plotTitle = sprintf('Contrast Sensitivity Function \n %.0f Trials, Sigma_{PSF} = %.0f, Mosaic Size = %.0f',expInfo.nTrialsNum,expInfo.psfSigma,round(expInfo.theMosaic.fov(1)));
end

figure()

plot(frequencyRange, sensitivity,'b.-','MarkerSize',20)

% Plotting data from Casagrande (1984), if indicated
if plotCasagrandeData
    hold on
    load csfData_Casagrande1984.mat ts_CSF_M
    
    ts1 = ts_CSF_M{1};
    ts2 = ts_CSF_M{2};
    ts3 = ts_CSF_M{3};
    
    plot(ts1(:,1),ts1(:,2),'k.-')
    plot(ts2(:,1),ts2(:,2),'kx-')
    plot(ts3(:,1),ts3(:,2),'ko-')
    
    legend('ISETTreeShrew CSF','Casagrande TS0','Casagrande TS1','Casagrande TS3','Location','northwest')
end

if max(frequencyRange) > 2
    set(gca,'xlim',[0.1,max(frequencyRange)])
    contrastTicks = [0.1 0.2 0.5 1.0 2.0, max(frequencyRange)];
    contrastTickLabels = {'.1', '.2', '.5', '1', '2', sprintf('%.0f',max(frequencyRange))};
else
    set(gca,'xlim',[0.1,2])
    contrastTicks = [0.1 0.2 0.5 1.0 2.0 ];
    contrastTickLabels = {'.1', '.2', '.5', '1', '2'};
end

set(gca, ...
    'ylim',[0,max(sensitivity)], ...
    'XTick', contrastTicks, ...
    'XTickLabel', contrastTickLabels, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'FontSize', 16)

xlabel('\it Spatial Frequency (Cycles/Degree)');
ylabel('\it Contrast Sensitivity');
title(plotTitle);

end
