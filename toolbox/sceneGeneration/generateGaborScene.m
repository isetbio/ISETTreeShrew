function scene = generateGaborScene(stimParams, varargin)
% Method to generate an ISETBio scene representing a Gabor stimulus
%
% Syntax:
%   scene = generateGaborScene(stimParams, [varargin])
%
% Description:
%    This function generates an ISETBio scene of a Gabor stimulus based on 
%    the passes stimulus parameters. We use a built-in scene generation 
%    method in ISETBio which generates sinusoidal images. This function 
%    requires a params struct which we generate based on the passed 
%    stimParams info
%
% Inputs:
%    stimParams                - Stimulus parameters for the Gabor
%
% Outputs:
%    scene                     - The ISETBio scene representing the Gabor
%
% Optional key/value pairs:
%   renderingDisplay           - Display object, If passed, the returned
%                                scene will be the scene as realized on
%                                that display
% History:
%    11/23/18  NPC  ISETBIO TEAM, 2018

%% parse input
p = inputParser;
p.addParameter('presentationDisplay', [], @validateDisplayArgument);
p.parse(varargin{:});
presentationDisplay = p.Results.presentationDisplay;

if (~isempty(presentationDisplay))
    stimParams = updateStimParamsForDisplay(stimParams, presentationDisplay);
end

% Transform stimParams into imageHarmonicParams
imageHarmonicParams = struct(...
        'freq', stimParams.sizeDegs * stimParams.spatialFrequencyCyclesPerDeg, ...
        'ang', stimParams.orientationDegs/180*pi, ...
        'ph', stimParams.phaseDegs/180*pi, ...
        'contrast', stimParams.contrast, ...
        'GaborFlag', stimParams.sigmaDegs/stimParams.sizeDegs, ...
        'row', stimParams.pixelsAlongHeightDim , ...
        'col', stimParams.pixelsAlongWidthDim);
    
% Generate a scene using these params
scene = sceneCreate('harmonic', imageHarmonicParams);

% Match viewing distance and mean luminance params
scene = sceneSet(scene, 'wangular', stimParams.sizeDegs);
scene = sceneAdjustLuminance(scene, stimParams.meanLuminanceCdPerM2);

if (~isempty(presentationDisplay))
    % Place the scene at the same distance as the display
    viewingDistanceMeters = displayGet(presentationDisplay, 'viewing distance');
    scene = sceneSet(scene, 'distance', viewingDistanceMeters);
    % Realize the scene into the presentation display
    scene = realizeSceneOnDisplay(scene, presentationDisplay);
end
end

% Method to compute the number of pixels for the stimulus given the
% stimulus size, and the viewing distance & pixel size of the display
function stimParams = updateStimParamsForDisplay(stimParams, presentationDisplay)
    % retrieve the display's pixel size 
    displayPixelSizeMeters = displayGet(presentationDisplay, 'sample spacing');
    % retrieve the display's viewing distance
    viewingDistanceMeters = displayGet(presentationDisplay, 'distance');
    % compute pixel size in visual degrees
    displayPixelSizeDegrees = ...
        2 * atand(0.5*displayPixelSizeMeters/viewingDistanceMeters);
    % divide by the stimulus size in degrees to get the pixels along the width
    stimParams.pixelsAlongWidthDim = ...
        round(stimParams.sizeDegs/displayPixelSizeDegrees(1));
    stimParams.pixelsAlongHeightDim = ...
        round(stimParams.sizeDegs/displayPixelSizeDegrees(2));
end


% Method to validation the passed display argument
function vStatus = validateDisplayArgument(x)

   vStatus = true;
   if isempty(x)
       % OK, we can have an empty argument and we do nothing
   elseif ~isstruct(x)
       % If non-empty it must be a struct object
       error('Input is not an ISETBio display object struct');
   else
       if ~isfield(x, 'type')
           error('Input is not an ISETBio display object struct');
       end
       if ~strcmp(x.type, 'display')
            error('Input is not an ISETBio display object struct');
       end
   end
end