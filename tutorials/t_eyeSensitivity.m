%% t_eyeSensitivity
% Illustrate basic principles of how eye parameters affect isomerizations.
%
% Description:
%    Demonstrate how changing the focal length, the pupil diameter 
%    and the inner segment aperture affect photoreceptor isomerization rate
%    computed by ISETBio, and compare this with the analytical analysis presented
%    Animal Eyes, pages 65-66 ff.
%
%   (Land, M. F., & Nilsson, D. E. (2012). Animal Eyes. OUP Oxford.)
%
% See also:
%

% History:
%   02/01/19 jsc  Wrote initial version.
%   02/08/19 jsc, dhb Documentation formatting etc.
%   04/09/19 jsc, Further work on intercept issues

%% Initialize workspace and close old figures
clear; close all;
ieInit;

%% Parameters
%
% How many data points to compute and plot
nPointsToCompute = 5;

% Base eye parameters.  These are reasonable starting points for a tree
% shrew eye.
baseFocalLengthMM = 4.35;
basePupilDiameterMM = 2.0;
baseInnerSegmentDiameterUM = 7.0;

% Level of change of varied parameter.  This also gets scaled
% by sqrt of integer change before being added to the base. 
% Numerical choices are not fundamental, just to give us a reasonable
% range.
deltaFocalLengthMM = 2.0;
deltaPupilDiameterMM = 1.0;
deltaInnerSegmentDiameterUM = 3.0;

% Specify which eye parameter we'll study in this run. The way this works
% is that the vector specifise which of the parameters will be incremented
% in the main loop below, where we recompute isomerizations across eye
% parameter variation.
%   [1,0,0] - Vary focal length
%   [0,1,0] - Vary pupil diameter
%   [0,0,1] - Vary inner segment diameter
v = [0,1,0];

% Specify cone densities similar to Peichl 1989 and set up a mosaic. Tree
% shrew mosaic dominated by L cones, so  we'll approximate with only L cones.
% For historical reasons, ISETBio parameters cone types in a vector
% "black", L, M and S. Tree shrews have no M cones. 
%
% Variable whichConeType determines which cone type we'll use
% to estimate isomerizations.  2 -> L, 3 -> M 4 -> S.  It would
% be a bad idea to use M, since we specify a mosaic with no M cones.
spatialLMSdensities = [0 1 0 0];
whichConeType = 2;

% Size of mosaic in degrees.
fovDegs = 0.4 * [1 1];

%% Initialize variables
meanRetinalIlluminance = zeros(1,nPointsToCompute + 1);
tMosaicExcitationMean = zeros(1,nPointsToCompute + 1);
eyeParameterValue = zeros(1,nPointsToCompute + 1);

% Analytic relative sensitivity as specified from the treatment in Animal
% Eyes. Basically, isomerizaiton rate should go up with the square of pupil
% diameter, up with the square of inner segment diameter, and down with the
% square of the focal length. These are the relations we are going to
% verify here.
s_AnimalEyes = zeros(1,nPointsToCompute + 1);
s_Iset = zeros(1,nPointsToCompute + 1);

%% Scene creation
%
% This one emits equal photon rates at all wavelengths
testScene = sceneCreate('uniformEqualPhoton');

%% Main loop.
%
% Compute isomerizations across parameter variation and save up the results
% for plotting.

%
% First, we will calculate the sensitivities when the parameter of interest
% is zero.

vT = ~ v;
focalLengthMM = vT(1)*baseFocalLengthMM;
pupilDiameterMM = vT(2)*basePupilDiameterMM;
innerSegmentDiameterUM = vT(3)*baseInnerSegmentDiameterUM;

[s_IsetZero,s_AnimalEyesZero] = getSensitivity(pupilDiameterMM,focalLengthMM,innerSegmentDiameterUM,testScene,spatialLMSdensities,fovDegs,whichConeType);    

s_Iset(1) = s_IsetZero;
s_AnimalEyes(1) = s_AnimalEyesZero;

for n = 1:nPointsToCompute
    % Get a vector that lets us decide the size of each of the parameters.  This
    % starts with v and increments on each loop iteration.
    vector = [n-1,n-1,n-1].*v;
        
    % Set up parameters for this itereation. We vary delta according to
    % sqrt of iteration for each parameter, just because we like the way
    % plots look when we do that.
    focalLengthMM = baseFocalLengthMM + sqrt(vector(1)) * deltaFocalLengthMM;
    pupilDiameterMM = basePupilDiameterMM + sqrt(vector(2)) * deltaPupilDiameterMM;
    innerSegmentDiameterUM = baseInnerSegmentDiameterUM + sqrt(vector(3)) * deltaInnerSegmentDiameterUM;
    
    % This switch statement sets up information for plotting, which depends
    % on which eye parameter we are varying. We also store the value of the
    % parameter for each interation, to be used in labeling individual
    % points in the plot.
    switch find(v)
        case 1
            eyeParameterValue(n + 1) = focalLengthMM;
            parameterName = 'Focal Length'; 
            shortParameterName = 'F'; 
            parameterUnits = 'mm';
        case 2
            eyeParameterValue(n + 1) = pupilDiameterMM;
            parameterName = 'Pupil Diameter';
            shortParameterName = 'P_D';
            parameterUnits = 'mm';
        case 3
            eyeParameterValue(n + 1) = innerSegmentDiameterUM;
            parameterName = 'Inner Segment Aperture Diameter';
            shortParameterName = 'IS_D';
            parameterUnits = 'um';
    end

    % Now, calculate the two measures of sensitivity for the parameters
    % chosen. The specific steps can be seen in the Functions section.
    [s_Iset(n + 1),s_AnimalEyes(n + 1)] = getSensitivity(pupilDiameterMM,focalLengthMM,innerSegmentDiameterUM,testScene,spatialLMSdensities,fovDegs,whichConeType);    
        
end

% Now, use pupilDiameter



%ISETBio's average cone excitation for a given image is conceptually the
%same as the "eye sensitivity" given in Animal Eyes

%% Plotting

x = s_AnimalEyes;
y = s_Iset;

ft = fitlm(x,y);
interceptCI = [ ft.Coefficients.Estimate(1) - 2 * ft.Coefficients.SE(1), ...
    ft.Coefficients.Estimate(1) + 2 * ft.Coefficients.SE(1)];

plot(x,y,'o')
hold on
line([0,max(x)],[ft.Coefficients.Estimate(1), ft.Coefficients.Estimate(1) + max(x) * ft.Coefficients.Estimate(2)])

xlabel('Sensitivity (Animal Eyes)')
ylabel('Sensitivty (ISETBio)')
xlim([0,inf])
ylim([0,inf])

title([{'Relationship Between ISETBIO and Animal Eyes Sensitivity'}, ...
    {sprintf('As %s Changes',parameterName)}])

labels = cell(1,nPointsToCompute + 1);
for i=1:(length(eyeParameterValue))
labels(i) = cellstr(sprintf('%s= %.2f %s', shortParameterName, eyeParameterValue(i), parameterUnits));
end

text(x,y,labels,'VerticalAlignment','top','HorizontalAlignment','right')
text(max(x)/2,max(y)/6,sprintf('Intercept CI = [%.2f,%.2f]', interceptCI(1),interceptCI(2)),'VerticalAlignment','top','HorizontalAlignment','left')

%% Functions

function [isetSensitivity, geomSensitivity] = getSensitivity(pupilDiameterMM, ...
    focalLengthMM,innerSegmentDiameterUM,testScene,spatialLMSdensities,fovDegs, ...
    whichConeType)        

% Create optical image object for using the pupil diameter and focal length.
    tOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM, 'focalLengthMM', ...
        focalLengthMM);
    
    % Create the mosaic using the inner segment diameter
    tMosaic = coneMosaicTreeShrewCreate(tOI.optics.micronsPerDegree, ...
        'spatialDensity', spatialLMSdensities, ...
        'customInnerSegmentDiameter', innerSegmentDiameterUM, ...
        'integrationTime', 5/1000, ...
        'fovDegs', fovDegs);
    
    % Compute the retinal image
    tOI = oiCompute(tOI, testScene);
    
    % Compute the mosaic responses (for more precise responses, use more trials)
    nTrialsNum = 1;
    emPath = zeros(nTrialsNum, 1, 2);
    
    % Compute *treeshrew* mosaic excitation responses to treeshrew optical image
    tMosaicExcitation = tMosaic.compute(tOI, 'emPath', emPath);
    
    % This function reshapes the mosaic excitation data, which is only necessary 
    % if nTrialsNum > 1. Either way this is just finding the mean photoreceptor 
    % excitations to the optical image, which correponds to the ISETBio version 
    % of eye sensitivity
    isetSensitivity = ...
        meanResponseToOpticalImage(tMosaic, tMosaicExcitation, whichConeType);
    
    % Calculate estimated sensitivity according to Animal Eyes
    geomSensitivity = 0.62 * (pupilDiameterMM^2 * innerSegmentDiameterUM^2)/ ...
        (focalLengthMM^2);
    
end

function meanResponse = meanResponseToOpticalImage(coneMosaic, coneMosaicResponse, ...
    targetConeType)
nTrialsNum = size(coneMosaicResponse,1);
coneMosaicResponse  = reshape(coneMosaicResponse, [nTrialsNum numel(coneMosaic.pattern)]);
idx = find(coneMosaic.pattern == targetConeType);
meanResponse = mean(mean(coneMosaicResponse(:,idx)));
end

