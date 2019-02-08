%% t_eyeGeometry
% Illustrate basic principles of how eye parameters affect isomerizations.
%
% Description:
%    Demonstrate how changing the focal length, the pupil diameter 
%    and the inner segment aperture affect photoreceptor isomerization rate
%    computed by ISETBio, and compare this with the analytical analysis presented
%    Animal Eyes (ADD BIBLIGRAPHIC INFO), pages XX ff.
%
% See also:
%

% History:
%   02/01/19 jsc  Wrote initial version.
%   02/08/19 jsc, dhb Documentation formatting etc.

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
v = [1,0,0];
v = [0,0,1];

% Specify cone densities similar to Peichl 1989 and set up a mosaic.
% For historical reasons, ISETBio parameters cone types in a vector
% "black", L, M and S. Tree shrews have no M cones. 
%
% Variable whichConeType determines which cone type we'll use
% to estimate isomerizations.  2 -> L, 3 -> M 4 -> S.  It would
% be a bad idea to use M, since we specify a mosaic with no M cones.
spatialLMSdensities = [0 .9 0 .1];
whichConeType = 2;

% Size of mosaic in degrees.
fovDegs = 0.4*[1 1];

%% Initialize variables
meanRetinalIlluminance = zeros(1,nPointsToCompute);
tMosaicExcitationMean = zeros(1,nPointsToCompute);
eyeParameterValue = zeros(1,nPointsToCompute);

% Analytic relative sensitivity as specified from the treatment in Animal
% Eyes. Basically, isomerizaiton rate should go up with the square of pupil
% diameter, up with the square of inner segment diameter, and down with the
% square of the focal length. These are the relations we are going to
% verify here.
s_AnimalEyes = zeros(1,nPointsToCompute);

%% Create a scene
%
% This one emits equal photon rates at all wavelengths
testScene = sceneCreate('uniformEqualPhoton');

%% Main loop.
%
% Compute isomerizations across parameter variation and save up the results
% for plotting.
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
            eyeParameterValue(n) = focalLengthMM;
            parameterName = 'Focal Length'; 
            shortParameterName = 'F'; 
            parameterUnits = 'mm';
        case 2
            eyeParameterValue(n) = pupilDiameterMM;
            parameterName = 'Pupil Diameter';
            shortParameterName = 'P_D';
            parameterUnits = 'mm';
        case 3
            eyeParameterValue(n) = innerSegmentDiameterUM;
            parameterName = 'Inner Segment Aperture Diameter';
            shortParameterName = 'IS_D';
            parameterUnits = 'um';
    end
    
    % Create optical image object for current parameters.
    tOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM, 'focalLengthMM', ...
        focalLengthMM);
    
    % Create the mosaic
    tMosaic = coneMosaicTreeShrewCreate(tOI.optics.micronsPerDegree, ...
        'spatialDensity', spatialLMSdensities, ...
        'customInnerSegmentDiameter', innerSegmentDiameterUM, ...
        'integrationTime', 5/1000, ...
        'fovDegs', fovDegs);
    
    % Compute the retinal image
    tOI = oiCompute(tOI, testScene);
    
    % Get and store the retinal irradiance
    meanRetinalIlluminance(n) = oiGet(tOI, 'mean illuminance');
    
    % Compute the mosaic responses
    nTrialsNum = 3;
    emPath = zeros(nTrialsNum, 1, 2);
    
    % Compute *treeshrew* mosaic excitation responses to treeshrew optical image
    tMosaicExcitation = tMosaic.compute(tOI, 'emPath', emPath);
    
    % Find mean excitations
    tMosaicExcitationMean(n) = ...
        meanResponseToOpticalImage(tMosaic, tMosaicExcitation, whichConeType);
%     
%     wT = tMosaic.wave;
%     tPigment = tMosaic.pigment;
%     treeShrewAbsorbance = tPigment.absorbance;
%     
%     %.2f mmigure();
%     %plotActionSpectra(wT, treeShrewAbsorbance, 'absorbance', 'treeshrew')
%     
%     axialDensities = tPigment.opticalDensity;
%     treeShrewAxialAbsorbance = treeShrewAbsorbance * diag(axialDensities);
%     
%     treeShrewAbsorptance = 1 - 10 .^ (-treeShrewAxialAbsorbance);
%     
%     peakEfficiencies = tPigment.peakEfficiency;
%     treeShrewQuantalEfficiency = treeShrewAbsorptance * diag(peakEfficiencies);
%     
%     treeShrewInnerSegmentArea = tPigment.pdArea*1e12;
%     treeShrewIntegratedQuantalEfficiency = treeShrewQuantalEfficiency * ...
%         treeShrewInnerSegmentArea;
%     dWT = wT(2)-wT(1);
%     treeshrewSpectrallyIntegratedQuantalEfficiencies = ...
%         sum(treeShrewIntegratedQuantalEfficiency,1) * dWT/sum(wT);
%      
%     % Get number of L/S cones
%     pattern = tMosaic.pattern;
%     edges = unique(pattern);
%     counts = histc(pattern(:), edges);
%     cone_count = [counts(2),0,counts(3)];
    
    %E(n,1:3) = treeshrewSpectrallyIntegratedQuantalEfficiencies;
    
    %find average QE over mosaic (weight by cone occurance then average)
%     E(n) = sum(cone_count.*treeshrewSpectrallyIntegratedQuantalEfficiencies) ...
%         /sum(cone_count);
    
    % Calculate estimated sensitivity according to Animal Eyes
    s_AnimalEyes(n) = 0.62 * (pupilDiameterMM^2 * innerSegmentDiameterUM^2)/ ...
        (focalLengthMM^2);
    
end

s_Iset = tMosaicExcitationMean;
x = s_AnimalEyes;
y = s_Iset;
plot(x,y,'o')
xlabel('Sensitivity (Animal Eyes)')
ylabel('Sensitivty (ISETBio)')

title([{'Relationship Between ISETBIO and Animal Eyes Sensitivity'}, ...
    {sprintf('As %s Changes',parameterName)}])

labels = cell(1,nPointsToCompute);
for i=1:nPointsToCompute
labels(i) = cellstr(sprintf('%s= %.2f %s', shortParameterName, eyeParameterValue(i), parameterUnits));
end

text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

I1 = meanRetinalIlluminance./meanRetinalIlluminance(1); %illuminance for each level

E1 = tMosaicExcitationMean./tMosaicExcitationMean(1); %efficiency for each cone (LMS) for each level

S1 = s_Iset./s_Iset(1);


%% Functions

function meanResponse = meanResponseToOpticalImage(coneMosaic, coneMosaicResponse, ...
    targetConeType)
nTrialsNum = size(coneMosaicResponse,1);
coneMosaicResponse  = reshape(coneMosaicResponse, [nTrialsNum numel(coneMosaic.pattern)]);
idx = find(coneMosaic.pattern == targetConeType);
meanResponse = mean(mean(coneMosaicResponse(:,idx)));
end
