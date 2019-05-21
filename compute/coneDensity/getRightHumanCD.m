function [rightAngHuman,rightConeDensityHuman] = getRightHumanCD(varargin)
% 
% Use the coneDensityReadData function to get cone Density as a function of
% visual angle for the right human eye

% Optional Parameters:
%   'rightEyeFOV'       FOV in degrees, numeric
%   'focalLength'       focal length in mm, numeric
%   'coneDensitySource' Source for cone density estimate
%       'Curcio1990'         From Figure 6 of Ref 1 below (default).
%       'Song2011Old'        From Table 1 of Ref 2 below, old subjects data.
%       'Song2011Young'      From Table 1 of Ref 2 below, young subjects data.


% Parse inputs
p = inputParser;
p.addParameter('rightEyeFOV', 160, @isnumeric)
p.addParameter('focalLength', 17, @isnumeric)
p.addParameter('coneDensitySource','Curcio1990',@ischar)

p.parse(varargin{:});
rightEyeFOV = p.Results.rightEyeFOV;
focalLength = p.Results.focalLength;
coneDensitySource = p.Results.coneDensitySource;

% Let's convert the cone density to a unit that can be compared across eye
% sizes: cones per visual angle. Basically, we need a measure of how many
% external degrees of visual angle is captured by 1 mm inside the eye. We
% need an approximation of the total degreees of the FOV, and the total mm
% of photoreceptors (approximated by pi*focal length)

conePerRetinalDistanceToConePerVisAngle = (focalLength*3.14)/rightEyeFOV;
 
% Calculate Cone Density
ang = (-rightEyeFOV/2):(rightEyeFOV/2);
rightConeDensityHumanTemp = zeros(1,rightEyeFOV);
for i = 1:length(ang)
        rightConeDensityHumanTemp(i) = coneDensityReadData('eccentricity', ...
            abs(ang(i)), 'eccentricityUnits', 'deg', 'whichEye', 'right','coneDensitySource',coneDensitySource);
end

rightConeDensityHuman = conePerRetinalDistanceToConePerVisAngle.*rightConeDensityHumanTemp;

rightAngHuman = ang;
% most datasets for human CD only measure for within ~10 degrees of the
% fovea, so we assume constant CD outside that range
rightConeDensityHuman(isnan(rightConeDensityHuman)) = min(rightConeDensityHuman);
end
