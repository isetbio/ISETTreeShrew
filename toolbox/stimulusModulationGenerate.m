function [stimulusTimeAxis, stimulusModulationFunction, temporalParams] = stimulusModulationGenerate(stimulusDurationSeconds,displayRefreshRate)
% Generate the timeAxis and temporal modulation function for the oiSequence
%
% Syntax:
%   [stimulusTimeAxis, stimulusModulationFunction, temporalParams] = ...
%     stimulusModulationGenerate(stimulusDurationSeconds,displayRefreshRate)
%
% Description:
%    Generate the timeAxis and temporal modulation function. These are used
%    by the oiSequence instantiation function which generates the spatio-
%    temporal sequence of optical images evoked by the stimulus
%
%
% Inputs:
%    stimulusDurationSeconds        Seconds the stimulus is presented for
%    displayRefreshRate             The refresh rate of the presentation display (Hz)
% 
% Outputs:
%    temporalParams                 Struct with all temporal parameters
%    stimulusTimeAxis               Temporal support of the stimulus modulation function
%    stimulusModulationFunction     The stimulus modulation function
%

% History:
%    1/24/2019  NPC   Wrote it

    % Seconds before stimulus presentation
    secondsForResponseStabilization = 50/1000;
    % Seconds after stimulus presentation
    secondsForResponseExtinction = 50/1000;
    
    stimulusSamplingIntervalInSeconds = 1.0/displayRefreshRate;
    stimulusDurationInSeconds = round(stimulusDurationSeconds/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    secondsForResponseStabilization = round(secondsForResponseStabilization/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    secondsForResponseExtinction = round(secondsForResponseExtinction/stimulusSamplingIntervalInSeconds)*stimulusSamplingIntervalInSeconds;
    
    temporalParams = struct(...
        'stimulusDurationInSeconds', stimulusDurationInSeconds, ...
        'stimulusSamplingIntervalInSeconds', stimulusSamplingIntervalInSeconds, ...
        'secondsForResponseStabilization', secondsForResponseStabilization, ...
        'secondsForResponseExtinction', secondsForResponseExtinction ...
    );
    [stimulusTimeAxis, stimulusModulationFunction] = squareTemporalWindowCreate(temporalParams);
end

function [sampleTimes, squareTemporalWindow, rasterModulation] = squareTemporalWindowCreate(temporalParams)

    stimulusSamples = round(temporalParams.stimulusDurationInSeconds/temporalParams.stimulusSamplingIntervalInSeconds);
    if (isfield(temporalParams, 'secondsForResponseStabilization'))
        stabilizingTimeSamples = round(temporalParams.secondsForResponseStabilization/temporalParams.stimulusSamplingIntervalInSeconds);
    else
        stabilizingTimeSamples = 0;
    end

    if (isfield(temporalParams, 'secondsForResponseExtinction'))
        extinctionTimeSamples = round(temporalParams.secondsForResponseExtinction/temporalParams.stimulusSamplingIntervalInSeconds);
    else
        extinctionTimeSamples = 0;
    end

    sampleTimes = 1:(stabilizingTimeSamples + stimulusSamples + extinctionTimeSamples);
    sampleTimes = sampleTimes - stabilizingTimeSamples - round(stimulusSamples/2);
    sampleTimes = sampleTimes * temporalParams.stimulusSamplingIntervalInSeconds;

    squareTemporalWindow = zeros(1,numel(sampleTimes));
    onTime = find(sampleTimes >= -temporalParams.stimulusDurationInSeconds/2);
    squareTemporalWindow(onTime(1)+(0:(stimulusSamples-1))) = 1;

    if (isfield(temporalParams, 'addCRTrasterEffect')) && (temporalParams.addCRTrasterEffect)
        % Add CRT raster effect
        phosphorFunction = crtPhosphorActivationFunction(1/temporalParams.stimulusSamplingIntervalInSeconds, temporalParams.rasterSamples);
        rasterSamples = numel(phosphorFunction.timeInSeconds);
        raster = zeros(1,numel(squareTemporalWindow)*rasterSamples);
        raster(1,1:rasterSamples:end) = squareTemporalWindow*0+1;
        rasterModulation = conv(raster, phosphorFunction.activation);
        rasterModulation = rasterModulation(1:numel(squareTemporalWindow)*rasterSamples);

        tmp = zeros(1,numel(squareTemporalWindow)*rasterSamples);
        for i = 1:numel(squareTemporalWindow)*rasterSamples-1
            tmp(i) = squareTemporalWindow(floor((i-1)/rasterSamples)+1);
        end
        squareTemporalWindow = tmp;
        sampleTimes = linspace(sampleTimes(1), sampleTimes(end), numel(squareTemporalWindow));
    else
        rasterModulation = [];
    end

    figure()
    plot(sampleTimes, squareTemporalWindow, 'ks-'); hold on;
    if (~isempty(rasterModulation ))
        plot(sampleTimes, rasterModulation, 'rs-');
    end

end

function phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
% phosphorFunction = crtPhosphorActivationFunction(refreshRate, samplesPerRefreshCycle) 
%
% Create a phoshor activarion function with sharp rise, shower decline
%
%  7/7/16  npc Wrote it.
%

    alpha = 1.9; t_50 = 0.02/1000; n = 2;
    phosphorFunction.timeInSeconds = linspace(0,1.0/refreshRate, samplesPerRefreshCycle);
    phosphorFunction.activation = (phosphorFunction.timeInSeconds.^n)./(phosphorFunction.timeInSeconds.^(alpha*n) + t_50^(alpha*n));
    phosphorFunction.activation = phosphorFunction.activation - phosphorFunction.activation(end);
    phosphorFunction.activation(phosphorFunction.activation<0) = 0;
    phosphorFunction.activation = phosphorFunction.activation / max(phosphorFunction.activation);
    figure(3);
    plot(phosphorFunction.timeInSeconds, phosphorFunction.activation, 'ko-');
    
end
