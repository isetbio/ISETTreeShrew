function [stimParamsDefault, displayDefault] = csfStimParamsDefault(~)
% Creates a default stimulus parameters for the gabor patches used in
% Casagrande's 1984 tree shrew paper (except smaller, Casagrande used 16
% degree patches). Vertical stripes, no overlaid gaussian (or: very large
% gaussian width), and a mean luminance of 35 cd/m^2. 

stimParamsDefault = struct(...
    'spatialFrequencyCyclesPerDeg', 1.5, ...% Spatial frequency (cpd)
    'contrast', .01, ...                    % Contrast (Michelson)
    'orientationDegs', 0, ...               % Vertical Orientation
    'sizeDegs', 5, ...                      % Mosaic size (degrees)
    'phaseDegs', 90, ...                    % spatial phase degrees, 0 = cos, 90 = sin
    'sigmaDegs', 100, ...                   % wide Gaussian envelope
    'meanLuminanceCdPerM2', 35, ...         % 35 cd/m2 mean luminance
    'pixelsAlongWidthDim', [], ...          % pixels- width dimension
    'pixelsAlongHeightDim', [] ...          % pixel- height dimension
    );

displayDefault = displayCreate('LCD-Apple', 'viewing distance', 20/100);

end
