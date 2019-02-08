%   t_eyeGeometry:
% In this code, I determine how changing the focal length, the pupil diameter 
%   and the inner segment aperture affect both the sensitivity as defined by 
%   Animal Eyes and the sensitivity as determined by multiplying Illuminance 
%   and Efficiency as determined through ISETBio

max_multiplier = 5; %added in code that made it slightly harder to change this, 
% as I currently define five labels specifically

I = zeros(1,max_multiplier);
E = zeros(1,max_multiplier);
change_Vector = zeros(1,max_multiplier);
s_AnimalEyes = zeros(1,max_multiplier);

v = [0,0,1]; %%%% focal length^2 (mm^2), pupil area (mm^2), aperture area (um^2)
% for example, [1,0,0] will change focal length^2 and keep other variables
% constant

for n = 1:max_multiplier
    vector = [n,n,n].*v;
    vector(vector==0) = 1;
    
    % Default treeshrew optics, mult * 2.0 mm pupil
    %default FL = 4.35
    focalLengthMM = sqrt(vector(1) * 4.35);
    
    pupilDiameterMM = sqrt(vector(2) * 2.0);
    
    innerSegmentDiameter = sqrt(vector(3) * 7.0);
    
    switch find(v,n)
        case 1
            change_Vector(n) = focalLengthMM;
            var = 'Focal Length';
            var_1 = 'F';
            units = 'mm';
        case 2
            change_Vector(n) = pupilDiameterMM;
            var = 'Pupil Diameter';
            var_1 = 'P_D';
            units = 'mm';
        case 3
            change_Vector(n) = innerSegmentDiameter;
            var = 'Inner Segment Aperture Diameter';
            var_1 = 'IS_D';
            units = 'um';
    end
    
    
    
    
    tOI = oiTreeShrewCreate('pupilDiameterMM', pupilDiameterMM, 'focalLengthMM', ...
        focalLengthMM);
    % Specify cone densities similar to Peichl 1989
    spatialLMSdensities = [0 .9 0 .1];
    % 0.4 x 0.4 deg mosaic for the treeshrew
    fovDegs = 0.4*[1 1];
    
    tMosaic = coneMosaicTreeShrewCreate(tOI.optics.micronsPerDegree, ...
        'spatialDensity', spatialLMSdensities, ...
        'customInnerSegmentDiameter', innerSegmentDiameter, ...
        'integrationTime', 5/1000, ...
        'fovDegs', fovDegs);
    % Step 2. Create a scene that emits equal photon rates at all wavelengths
    % and compute the optical image for the human and treeshrew optics
    
    % Create the equal photon rate test scene
    testScene = sceneCreate('uniformEqualPhoton');
    
    % Compute the retinal images
    tOI = oiCompute(tOI, testScene);
    
    % Retrieve the retinal irradiances
    meanIlluminanceTreeShrewRetina = oiGet(tOI, 'mean illuminance');
    
    %output retinal illuminance
    I(n) = meanIlluminanceTreeShrewRetina;
    
    % Compute the mosaic responses
    nTrialsNum = 3;
    emPath = zeros(nTrialsNum, 1, 2);
    
    % Compute *treeshrew* mosaic excitation responses to treeshrew optical image
    tMosaicExcitation = tMosaic.compute(tOI, 'emPath', emPath);
    
    % Find mean excitations
    coneType = 2;
    tMosaicExcitationMean = ...
        meanResponseToOpticalImage(tMosaic, tMosaicExcitation, coneType);
    
    wT = tMosaic.wave;
    tPigment = tMosaic.pigment;
    treeShrewAbsorbance = tPigment.absorbance;
    
    %.2f mmigure();
    %plotActionSpectra(wT, treeShrewAbsorbance, 'absorbance', 'treeshrew')
    
    axialDensities = tPigment.opticalDensity;
    treeShrewAxialAbsorbance = treeShrewAbsorbance * diag(axialDensities);
    
    treeShrewAbsorptance = 1 - 10 .^ (-treeShrewAxialAbsorbance);
    
    peakEfficiencies = tPigment.peakEfficiency;
    treeShrewQuantalEfficiency = treeShrewAbsorptance * diag(peakEfficiencies);
    
    treeShrewInnerSegmentArea = tPigment.pdArea*1e12;
    treeShrewIntegratedQuantalEfficiency = treeShrewQuantalEfficiency * ...
        treeShrewInnerSegmentArea;
    dWT = wT(2)-wT(1);
    treeshrewSpectrallyIntegratedQuantalEfficiencies = ...
        sum(treeShrewIntegratedQuantalEfficiency,1) * dWT/sum(wT);
    
    %Get number of L/S cones
    pattern = tMosaic.pattern;
    edges = unique(pattern);
    counts = histc(pattern(:), edges);
    cone_count = [counts(2),0,counts(3)];
    
    %E(n,1:3) = treeshrewSpectrallyIntegratedQuantalEfficiencies;
    
    %integrate QE over mosaic (so, weight by cone occurance then sum)
    E(n) = sum(cone_count.*treeshrewSpectrallyIntegratedQuantalEfficiencies) ...
        /sum(cone_count);
    
    s_AnimalEyes(n) = 0.62 * (pupilDiameterMM^2 * innerSegmentDiameter^2)/ ...
        (focalLengthMM^2);
    
end



s_Iset = E.*I;
x = s_AnimalEyes;
y = s_Iset;
plot(x,y,'o')
xlabel('Sensitivity (Animal Eyes)')
ylabel('Sensitivty (ISETBio)')

title([{'Relationship Between ISETBIO and Animal Eyes Sensitivity'}, ...
    {sprintf('As %s Changes',var)}])

labels = {sprintf('%s= %.2f %s', var_1, change_Vector(1), units), ...
    sprintf('%s= %.2f %s', var_1, change_Vector(2), units), ...
    sprintf('%s= %.2f %s', var_1 ,change_Vector(3), units), ...
    sprintf('%s= %.2f %s', var_1, change_Vector(4), units) ...
    sprintf('%s= %.2f %s', var_1, change_Vector(5), units)};

text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')

I1 = I./I(1); %illuminance for each level

E1 = E./E(1); %efficiency for each cone (LMS) for each level

S1 = s_Iset./s_Iset(1);

clear('max_multiplier','v','vector','focalLengthMM','pupilDiameterMM', ...
    'innerSegmentDiameter','tOI','spatialLMSdensities','fovDegs','tMosaic', ...
    'testScene','tOI','meanIlluminanceTreeShrewRetina','nTrialsNum','emPath',...
    'tMosaicExcitation','tMosaicExcitationMean','wT','tPigment','treeShrewAbsorbance', ...
    'axialDensities','treeShrewAbsorptance','peakEfficiencies','treeShrewQuantalEfficiency', ...
    'treeShrewInnerSegmentArea','treeShrewIntegratedQuantalEfficiency','dWT', ...
    'treeshrewSpectrallyIntegratedQuantalEfficiencies','pattern','edges', ...
    'counts','cone_count','coneType')

%% Functions

function meanResponse = meanResponseToOpticalImage(coneMosaic, coneMosaicResponse, ...
    targetConeType)
nTrialsNum = size(coneMosaicResponse,1);
coneMosaicResponse  = reshape(coneMosaicResponse, [nTrialsNum numel(coneMosaic.pattern)]);
idx = find(coneMosaic.pattern == targetConeType);
meanResponse = mean(mean(coneMosaicResponse(:,idx)));
end
