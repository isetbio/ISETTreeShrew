function [rightAngTS,rightConeDensityTS] = getRightTreeshrewCD(varargin)
% 
% Use data from Muller, 1989 to get cone Density as a function of
% visual angle for the right tree shrew eye

% Optional Parameters:
%   retinaChoice - char - 'A','B','C'. which flattened retina to use
%   rightEyeFOV - FOV in degrees, numeric
%   focalLength - focal length in mm, numeric
%   axialLength - axial length in mm, numeric
%   angleShift - angle that eye is shifted away from forward, numeric
%   interpolMethod - recommend 'spline' for smooth curve, 'pchip' if you
%   don't want to extrapolate past gradient lines from flattened retina

% Parse inputs
p = inputParser;
p.addParameter('retinaChoice', 'A', @ischar)
p.addParameter('rightEyeFOV', 175, @isnumeric)
p.addParameter('focalLength', 5.81, @isnumeric)
p.addParameter('axialLength', 7.8, @isnumeric)
p.addParameter('angleShift',58,@isnumeric)
p.addParameter('interpolMethod','spline',@ischar)

p.parse(varargin{:});
rightEyeFOV = p.Results.rightEyeFOV;
focalLength = p.Results.focalLength;
retinaChoice = p.Results.retinaChoice;
axialLength = p.Results.axialLength;
angleShift = p.Results.angleShift;
interpolMethod = p.Results.interpolMethod;

%% Data is all right here

% Let's convert the cone density to a unit that can be compared across eye
% sizes: cones per visual angle. Basically, we need a measure of how many
% external degrees of visual angle is captured by 1 mm inside the eye. We
% need an approximation of the total degreees of the FOV, and the total mm
% of photoreceptors (approximated by pi*focal length)

conePerRetinalDistanceToConePerVisAngle = (focalLength*3.14)/rightEyeFOV;

switch retinaChoice
    case 'A'
        exp_eye = 'right';
        ts_fixed_loc = [0, 1, 3, 6, 9, 10, 12, 14, 15]; %mm
        coarseTreeShrewCD = 1000.*[16, 20, 24, 28, 32, 32, 28, 24, 20]; %cones/mm^2
    case 'B'
        exp_eye = 'left';
        ts_fixed_loc = [0, 0.5, 2.5, 3, 3.5, 5,7, 8, 10, 11.5, 13, 14, 14.25, 14.5, 15]; %mm 
        coarseTreeShrewCD = 1000.*[16, 20, 24, 28, 32, 32, 32, 28, 28, 32, 32, 28, 24, 20, 16]; %cones/mm^2
    case 'C'
        exp_eye = 'left';
        ts_fixed_loc = [0,0.5,1,2,3,4.5,8,9,12,13,14,15]; %mm
        coarseTreeShrewCD = 1000.*[16, 20, 24, 28, 32, 36, 36, 32, 28, 24, 20, 16]; %cones/mm^2
end

coarseTreeShrewCD = conePerRetinalDistanceToConePerVisAngle.*coarseTreeShrewCD;

%%
% Looking at the figures: where is the center of the pupil? Not explicitly
% stated in paper, but seems to be in center of slice.
central_Loc = (max(ts_fixed_loc)-min(ts_fixed_loc))/2;

% Now, we want the locations to be with respect to the central location. We
% want the left side of the visual field to be negative, and the right side
% to be positive. We therefore need to discriminate between the left and
% right eye.

switch exp_eye
    case 'right'
        mmFromCenter = ts_fixed_loc - central_Loc;
    case 'left'
        mmFromCenter = central_Loc - ts_fixed_loc;
end

%% Back of envelope determination of visual angle
%
%Tranforming mm from flattened dimensions to angle wrt focal point)

% Determine angle from center in radians, treating mmFromCenter as distance
% along perimeter of circle in mm
angleFromCenter = (mmFromCenter ./ (axialLength/2));

% Project this angle and point on retina into coordinates in mm from center of
% eye
eccProj = (axialLength/2) .* sin(angleFromCenter);
axialProj = (axialLength/2) .* cos(angleFromCenter);

% We're interested in the geometry of the point on the retina with respect
% to the focal point of the eye, so we calculate the coordinates with
% respect to the focal point. Then, we find the angle from the focal point
% to the point on the retina.

focal_dist = axialProj - (axialLength/2) + focalLength;
vis_degrees = atand(eccProj ./ focal_dist) + angleShift;

%% Interpolation
%
% Now, use the interp1 function to create a discrete vector of cone
% density per visual angle as a function of degree eccentricity.

% Visual degrees:
rightAngTS = min(round(vis_degrees)):(10^(-1)):max(round(vis_degrees));
% Cone density:
rightConeDensityTS = interp1(vis_degrees,coarseTreeShrewCD,rightAngTS,interpolMethod);

end