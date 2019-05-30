function t_BinarySearchCSF(varargin)
%
%% t_BinarySearchCSF (Human/Treeshrew)
% Use binarySearchOverFeature, plotBinarySearch and plotCSF to generate a
% contrast sensitivity function.

%
% Syntax:
%   t_BinarySearchCSF()
%   t_BinarySearchCSF('dataName','bestDataEver')
%   t_BinarySearchCSF('csfDataToPlot',load('treeshrew_csf_250_trials.mat'))
%

% Description:
%   Simulating Casagrande's 1984 CSF method (used with
%   treeshrews) in ISETBio. However, the steps can be used for simulating
%   the task on humans optics as well. The key steps are as follows:
%
%   1) Choose the specific spatial frequencies for which you will
% determine the sensitivity to contrast. Currently, ISETBio can only model
% the high-frequency CSF dropoff, not the low-frequency dropoff. Casagrande
% found that the high-frequency dropoff for tree shrews was between 0.75
% and 2 cycles/degree.
%
%   2) For each spatial frequency, determine the contrast level (AKA
% the contrast threshold) that causes the tree shrews to accurately discern
% the contrast 75% of the time.
%
%     * In ISETBio, this step is accomplished by using a support vector
%       machine (SVM) approach to determine the fraction of the time that
%       the algorithm can detect a difference in the mosaic excitation
%       between a null scene and a scene with a contrast. A binary search
%       algorithm is used to identify the specific contrast that leads to
%       75% of the SVM predictions being accurate.
%
%       The binary search works as follows: The spatial frequency of the
%       gabor stimulus is set at the beginning of the search, and is held
%       constant throughout. A minimum and maximum contrast is chosen, with
%       the assumption that the 75% contrast is within the chosen range.
%       The steps involved in the binary search process are well-described
%       in https://en.wikipedia.org/wiki/Binary_search_algorithm.
%
%   3) Calculate the contrast sensitivity at each spatial frequency by
% taking the reciprocal of the threshold contrast.
%
%   4) Plot the sensitivity as a function of the spatial frequency.
%
% Optional key/value pairs:
%   csfData - Character array or structure input referencing previously
%       gathered data, if you are specifically interested in plotting that.
%   dataName - If you have a specific choice for the name of the generated
%       data, you can input dataName as a character array.
%   overwrite - If there is already a dataset with the same name, do you
%       want to overwrite it? true/false.
%   saveToWS - If you're just exploring the code, you may want to save the
%       generated data to the workspace after it is created. Logical.

% Tools used:
%   binarySearchOverContrast
%   plotBinarySearch
%   plotCSF
%   getSVMAcc

% See also:
%   t_psychometricFromBinarySearch

% History:
%   03/14/19 jsc  Wrote initial version.

%% Initialize workspace and close old figures
%clear; close all;
%ieInit;

%% What do you want to do?
%
% If you're not interested in computing a new CSF, but instead want to
% simply plot prevously gathered data, then input either the name of the
% CSF dataset or the data structure itself into the function. Otherwise,
% leave the input field empty.

%% Parameters

% What species are you interested in simulating? Choose 'treeshrew' or
% 'human'.
species = 'human';

% If you're simulating tree shrews, what do you want the cone density
% of your mosaic to be? This is controlled by changing the minimum cone
% separation.
%   - min cone separation of 6 um     ->      ~ 32,000 cones/mm^2
%   - min cone separation of 7.5 um   ->      ~ 22,000 cones/mm^2
%   - min cone separation of 8.5 um   ->      ~ 16,000 cones/mm^2
cone_spacing = 7.5; %um

% What do you want the standard deviation of the point-spread function to
% be? For tree shrews, our preliminary analyses have indicated the value is
% around 12 um.
psfSigma = 12; %um


% How large do you want the stimulus, and therefore the activated cone
% mosaic, to be? The Casagrande paper used a stimulus of ~14 x 14
% degrees. However, this would take a very long time to run the SVM, so
% we recommend using a size of between 3x3 and 6x6 degrees. For testing,
% however, you can use an even smaller mosaic.
sizeDegs = 3; % degrees per side

% What discrete spatial frequencies do you want to find the
% sensitivities for? The Casagrande paper showed a treeshrew CSF
% dropoff in the range of 0.75 - 2 cycles/degree, while for humans it
% is in the range of 1 - 20 cycles/degree.
frequencyRange = 0.75:0.25:2; %cycles per degree

% What range of contrasts do you want to explore for the spatial
% frequencies? Make it wide enough to include accuracies other than
% 100% and 50%!
contrastRange = [.001,.03];

% For the binary search, what range do you want it to reach before
% stopping?
desiredAccRange = [74.5,75.5];

% For the binary search, how many iterations do you want before
% stopping? As appealing as it may be to set a high number here, it is
% important to remember that a mere 10 iterations will lead to 1/1024
% of the original search range- likely past the precision of the SVM.
maxCycles = 8;

% When plotting, it is important to remember that the sensitivity of the
% ideal observer will be much higher than the behaviorally-observed
% sensitivity. This has been shown to be at least a factor of 20. The
% parameter toDivide simply shifts the sensitivity function downward before
% plotting.
toDivide = 20;

%% SVM Specifications:

% How much data do you want to use to train the SVM?
nTrialsNum = 250;

% How many folds do you want to use for cross-validation?
nFolds = 10;

%% Parse optional input from function

p = inputParser;
p.StructExpand = false;
p.addParameter('csfDataToPlot', struct, @isstruct)
p.addParameter('dataName', char.empty, @ischar)
p.addParameter('overwrite', false, @islogical)
p.addParameter('saveToWS', false, @islogical)

p.parse(varargin{:});
csfDataToPlot = p.Results.csfDataToPlot;
dataName = p.Results.dataName;
overwrite = p.Results.overwrite;
saveToWS = p.Results.saveToWS;

compute = true;

% Load a previously generated set of data, if indicated.
if ~ismember('csfDataToPlot',p.UsingDefaults)
    compute = false;
    if ~isempty(dataName)
        warning('You are loading data instead of generating new data, so no new data will be saved under your chosen name.')
    end
end

if compute
    
    %% Initialize Variables
    %
    % Initialize cell array to store approximated psychometric function data
    contrastsTotal = cell(1,length(frequencyRange));
    accuraciesTotal = cell(1,length(frequencyRange));
    
    % Initialize vector to store threshold contrasts for each spatial
    % frequency
    thresholdContrasts = zeros(1,length(frequencyRange));
    finalAccuracy = zeros(1,length(frequencyRange));
    finalSE = zeros(1,length(frequencyRange));
    
    %% Create a cone mosaic
    %
    % Create a cone mosaic of the appropriate size and cone density. If
    % you're interested in modeling the optics of another species, this is
    % the first place you would make changes.
    
    % For either species, you're going to want a) to set the psfSigma in the
    % optical image structure by passing 'inFocusPSFsigmaMicrons' through
    % speciesOI and b) set the cone mosaic size by passing 'fovDegs'
    % through speciesMosaic.
    
    % For treeshrews specifically, you might want to change the cone
    % spacing by passing 'customLambda' through speciesMosaic as well.
    
    switch species
        
        case 'treeshrew'
            theMosaic = coneMosaicTreeShrewCreate(75, ...   % microns/degree for scaling
                'fovDegs', sizeDegs, ...                    % match mosaic width to stimulus size
                'customLambda', cone_spacing);              % what's the minimum distance between the cones?
            theOI = oiTreeShrewCreate(...
                'inFocusPSFsigmaMicrons', psfSigma ...
                );
            
        case 'human'
            theMosaic = coneMosaicHex(7, ...               % hex lattice sampling factor
                'fovDegs', sizeDegs, ...                   % match mosaic width to stimulus size
                'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
                'eccBasedConeQuantalEfficiency', true, ... % cone quantal efficiency varies with eccentricity
                'integrationTime', 10/1000, ...            % 10 msec integration time
                'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations
            theOI = oiCreate('wvf human');
           
        otherwise
            error('species should be treeshrew or human')
            
    end
    
    %% Set stimulus parameters and display type to default for the csf experiment
    [stimParams,presentationDisplay] = csfStimParamsDefault();
    
    % We don't need the default contrast or frequency, because
    % we will set the frequency and find a contrast based on it
    stimParams = rmfield(stimParams,{'contrast','spatialFrequencyCyclesPerDeg'});
    stimParams.sizeDegs = sizeDegs;
    
    %% Main loop
    %
    % Begin looping over the discrete spatial frequencies chosen

    % Matlab generates warnings when you create variables during parallel
    % for loops, so we suppress this specific warning during the loop.
    %warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')
    
    count = 1;
    N = length(frequencyRange);
    
    D = parallel.pool.DataQueue;
    
    w = waitbar(0, 'CSF Progress');
    afterEach(D, @nUpdateWaitbar);
    
    T = tic;
    
    parfor i = 1:length(frequencyRange)
        
        stimParamsTemp = stimParams;
        stimParamsTemp.spatialFrequencyCyclesPerDeg = frequencyRange(i);
        
        %% Binary Search
        
        % Perform the search, using the binarySearchOverFeature function
        searchResults = binarySearchOverFeature(...
            'contrast', ...
            contrastRange, ...
            theMosaic, ...
            theOI, ...
            'presentationDisplay',presentationDisplay, ...
            'nTrialsNum',nTrialsNum, ...
            'nFolds',nFolds, ...
            'stimParams',stimParamsTemp, ...
            'maxCycles',maxCycles, ...
            'desiredAccRange',desiredAccRange ...
            );
        
        % Save variables
        contrastsTotal{1,i} = searchResults.samples;
        accuraciesTotal{1,i} = searchResults.performances;
        thresholdContrasts(1,i) = searchResults.thresholdFeature;
        finalAccuracy(1,i) = searchResults.finalAccuracy;
        finalSE(1,i) = searchResults.finalSE;
        
        send(D,i)
        
    end
    
    delete(gcp)
    close(w)
    
    time = toc(T)/60;
    disp(['The binary search took ', num2str(time), ' minutes.'])
    
    %% Save Data
    
    if isempty(dataName)
        csfDataToSave = [species,sprintf('_csf_%.0f_trials',nTrialsNum)];
    elseif endsWith(dataName,'.mat')
        csfDataToSave = erase(dataName,'.mat');
    else
        csfDataToSave = dataName;
    end
    
    % Will the data be overwritten?
    if ~overwrite
        k = 2;
        csfDataToSaveTemp = [csfDataToSave,'.mat'];
        while exist(csfDataToSaveTemp, 'file')
            csfDataToSaveTemp = [csfDataToSave,'_v',num2str(k),'.mat',];
            k = k + 1;
        end
        csfDataToSave = csfDataToSaveTemp;
    end
    
    %create binary search data structures
    
    expInfo.theMosaic = theMosaic;
    expInfo.stimParams = stimParams;
    expInfo.theOI = theOI;
    expInfo.presentationDisplay = presentationDisplay;
    expInfo.psfSigma = psfSigma;
    
    expInfo.nTrialsNum = nTrialsNum;
    expInfo.nFolds = nFolds;
    
    binaryResults.frequencyRange = frequencyRange;
    binaryResults.contrastRange = contrastRange;
    
    binaryResults.contrasts = contrastsTotal;
    binaryResults.accuracies = accuraciesTotal;
    
    binaryResults.thresholdContrasts = thresholdContrasts;
    binaryResults.finalAccuracy = finalAccuracy;
    binaryResults.finalSE = finalSE;
    
    binaryResults.time = time;
    
    if saveToWS
        assignin('base','binaryResults',binaryResults)
        assignin('base','expInfo',expInfo)
    end
    
    save(csfDataToSave,'species','expInfo','binaryResults')
else
    % If you chose to input 'csfDataToPlot' instead of generating new data,
    % extract the simulated experiment details and search results.
    expInfo = csfDataToPlot.expInfo;
    binaryResults = csfDataToPlot.binaryResults;
end

%% Plot the contrast sensitivity data
%
% Visualize relationship between contrast and accuracy.

% First, plot the iterations of each binary search in order to make sure
% the SVM results are monotonic. The threshold contrasts are marked by
% asterisks.

plotBinarySearch(binaryResults)

% If the SVM results are monotonic, it means that the ~75% contrast found
% in the binary search is the only contrast with that accuracy. Now, we can
% plot each threshold contrast as a function of the spatial frequency.

% Do you want to plot Casagrande's data as well?
%plotCasagrandeData = true;

plotCasagrandeData = true;

plotCSF(binaryResults, 'expInfo', expInfo, 'toDivide', toDivide, 'plotCasagrandeData', plotCasagrandeData)

%% Functions

    function nUpdateWaitbar(~)
        waitbar(count/N, w);
        count = count + 1;
    end

end