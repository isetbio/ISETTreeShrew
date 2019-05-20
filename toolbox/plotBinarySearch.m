function plotBinarySearch(binaryResults)
% Plot the contrast sensitivity function for a given simulated visual
% system front-end.
%
% Syntax:
%   plotBinarySearch(binaryResults)
%
% Inputs:
%   binaryResults - Matlab structure with fields 'frequencyRange',
%   'thresholdContrasts', 'finalAccuracy', 'contrasts',
%   'accuracies', and 'finalSE'

% See also:
%   t_BinarySearchCSF

% History:
%   04/13/19 jsc  Wrote initial version.

% Determine the total contrast range
maxCont = max(cellfun(@(x) max(x),binaryResults.contrasts));
minCont = min(cellfun(@(x) min(x),binaryResults.contrasts));

figure()
hold on

color = zeros(length(binaryResults.frequencyRange),3);
for i = 1:length(binaryResults.frequencyRange)
    color(i,:) = rand(1,3);
end

dataLen = length(binaryResults.frequencyRange);

C = {'r',[1 0.5 0],'y','g','b','c','m'};
if dataLen > 6
    Ctemp = cell([1 dataLen - 6]);
    for i = 1 : dataLen - 6
        Ctemp{i} = rand(1,3);
    end
    C = [C Ctemp];
end


% Plot the individual search points and their corresponding SVM accuracy
for j = 1:length(binaryResults.frequencyRange)
    tempMat = [binaryResults.contrasts{1,j};binaryResults.accuracies{1,j}];
    mat = sortrows(tempMat');
    plot(mat(:,1),mat(:,2),'LineWidth',2,'color',C{j})
    % Use the SE from the SVM results as an error bar for the final
    % accuracy
    errorbar(binaryResults.thresholdContrasts(j),binaryResults.finalAccuracy(j), binaryResults.finalSE(j),'k.','MarkerSize',10)
    text(binaryResults.thresholdContrasts(j)*1.05,binaryResults.finalAccuracy(j),['f = ' num2str(binaryResults.frequencyRange(j)) newline ' cpd'])
end

% Expand the x-axis if the max contrast is higher than 0.03
if maxCont > .031
    set(gca,'xlim',[minCont,maxCont])
    contrastTicks = [0.005,0.01 0.02 0.03, maxCont ];
    contrastTickLabels = {'0.005','.01', '.02', '.03', sprintf('%.0f',maxCont)};
else
    set(gca,'xlim',[minCont,.03])
    contrastTicks = [0.005 0.01 0.02 0.03];
    contrastTickLabels = {'0.005', '.01',  '.02',  '.03'};
end

set(gca, ...
    'XTick', contrastTicks, ...
    'XTickLabel', contrastTickLabels, ...
    'XScale', 'log', ...
    'FontSize', 16 ...
    );

xlabel('Contrast (Michelson)');
ylabel('Performance (% Accuracy)');
title(['Psychometric Function Approximations' newline 'Using a Binary Search']);

hold off
end