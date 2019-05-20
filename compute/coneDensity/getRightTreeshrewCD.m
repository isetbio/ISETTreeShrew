function [rightAngTS,rightConeDensityTS] = getRightTreeshrewCD(ts_choose)
% 
% Use the coneDensityReadData function to get cone Density as a function of
% visual angle for the right human eye

%%
% Parameters
% Axial length
ts_AL = 7.8; %mm 

% Posterior focal length
ts_FL = 5.81;

% What is the angular displacement for the treeshrew eyes? Barnes Jannuzi
% measured the angle of the skull as 58 degrees.
angleShift = 58; %degrees

ts_FOV = 175; %degrees

% 0/1/2: step size of 1.0/0.1/0.01... 1 -> step of 0.1 is usally enough.
fineness = 1;

% How do you want to interpolate the experimental data? Muller has gradient
% lines, but we want to fill in the rest. I recommend 'pchip' if you don't
% want to extrapolate past the experimental data, 'spline' if you want a
% smoother cone density curve over eccentricity.
interpolationMethod = 'spline';

% These vectors are directly adapted from Muller et al 1989. Locations in
% mm are with respect to left side of the slice.

% In order to get units of cones/mm to units of cones per visual angle, we
% need to multiply by some internal distance (in this case, an
% approximation for the size of the retina) and divide by an external
% correspondent (here, the total field of view of the eye).

switch ts_choose
    case 'A'
        exp_eye = 'right';
        ts_fixed_loc = [0, 1, 3, 6, 9, 10, 12, 14, 15]; %mm 
        ts_CD_coarse = (((3.14 * ts_FL)/ts_FOV)*1000).*[16, 20, 24, 28, 32, 32, 28, 24, 20]; %cones/degree
    case 'B'
        exp_eye = 'left';
        ts_fixed_loc = [0, 0.5, 2.5, 3, 3.5, 5,7, 8, 10, 11.5, 13, 14, 14.25, 14.5, 15]; %mm 
        ts_CD_coarse = (((3.14 * ts_FL)/ts_FOV)*1000).*[16, 20, 24, 28, 32, 32, 32, 28, 28, 32, 32, 28, 24, 20, 16]; %cones/degree
    case 'C'
        exp_eye = 'left';
        ts_fixed_loc = [0,0.5,1,2,3,4.5,8,9,12,13,14,15]; %mm 
        ts_CD_coarse = (((3.14 * ts_FL)/ts_FOV)*1000).*[16, 20, 24, 28, 32, 36, 36, 32, 28, 24, 20, 16]; %cones/degree
end

% (range(ts_fixed_loc)/2*1000/(angleShift+binocularField/2))
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
        rel_R = ts_fixed_loc - central_Loc;
    case 'left'
        rel_R = central_Loc - ts_fixed_loc;
end

% Back of envelope determination of visual angle (tranforming mm from
% flattened dimensions to wrapping around sphere, then finding angle)
a = (rel_R ./ (ts_AL/2));
proj_ecc = (ts_AL/2) .* sin(a);
dist = (ts_AL/2) .* cos(a);
focal_dist = dist - (ts_AL/2) + ts_FL;
vis_degrees = atand(proj_ecc ./ focal_dist) + angleShift;

% Now, use the interp1 function to create a discrete vector of cone density
% as a function of visual angle.

% Visual degrees:
rightAngTS = min(round(vis_degrees)):(10^(-1*fineness)):max(round(vis_degrees));
% Cone density:
rightConeDensityTS = interp1(vis_degrees,ts_CD_coarse,rightAngTS, interpolationMethod);
end