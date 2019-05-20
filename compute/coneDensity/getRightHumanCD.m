function [rightAngHuman,rightConeDensityHuman] = getRightHumanCD(~)
% 
% Use the coneDensityReadData function to get cone Density as a function of
% visual angle for the right human eye

% Parameters:

% In order to get units of cones/mm to units of cones per visual angle, we
% need to multiply by some internal distance (in this case, an
% approximation for the size of the retina) and divide by an external
% correspondent (here, the total field of view of the eye).

singleEyeFOV = 160;
focalLength = 17;
 
% Calculate Cone Density
ang = (-singleEyeFOV/2):(singleEyeFOV/2);
rightConeDensityHuman = zeros(1,singleEyeFOV);
for i = 1:length(ang)
        rightConeDensityHuman(i) = ((focalLength*3.14)/singleEyeFOV).*coneDensityReadData('eccentricity', ...
            abs(ang(i)), 'eccentricityUnits', 'deg', 'whichEye', 'right');
end

rightAngHuman = ang;
rightConeDensityHuman(isnan(rightConeDensityHuman)) = min(rightConeDensityHuman);
end
