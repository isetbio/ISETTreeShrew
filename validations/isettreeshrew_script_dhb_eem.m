% Replicate as best we can the results in
% Sajdak et al., 2019, Experimental Eye Research,
% 185, 107683.
%
% Thanks to Austin Roorda for clarifications and for
% providing tabulated data and example calculations.
%
% In this version, we unpack some of the encapsulated
% tree shrew routines to expose more parameters.

% History:
%   07/22/2022  dhb  Move from original Nicolas/Emily live script versions.
%   07/26/2022  eem  Move from David script to simplify/add cone mosaic.

% Configuration.  Execute the tbUseProjec command
% below to configure, if you use ToolboxToolbox.
%{
    tbUseProject('ISETTreeShrew');
%}

%% Initialize
clear; close all;

%% Unpacked oiTreeShrewCreate

% Specify wavelengths to compute on, and which individual
% tree shrews to analyze here.
targetWavelength = 550;

% Normally this is a range of wavelengths, but we are only
% trying to make sense of results at one wavelength here,
% so just do one for speed.
wavelengthSupport = targetWavelength;

% Can specify any contiguous range between 1 and 11.
%
% There were 11 shrews measured, and we have from Roorda the
% tabulated Zernike coefficients in an Excel spreadsheet.
TSindex = 1;

% Can turn on/off which figures to plot
plotFigureOfMerit = false;
plotPsfOtf = false;
plotMTFCompare = true;

% Figure of merit to optimize.
% Choices are:
%   'None'
%   'PSFPeak';
%   'MTFArea';
whichFigureOfMerit = 'PSFPeak';
if (~strcmp(whichFigureOfMerit,'None')) && plotFigureOfMerit
    PSFPeakFig = figure; clf;
    set(gcf,'Position',[1224 460 1126 1129]);
    MTFAreaFig = figure; clf;
    set(gcf,'Position',[1224 460 1126 1129]);
end

% Pupil diameter
%
% We're trying to reproduce Figure 2,
% which is for a 4 mm pupil. So we calculate for that.
% Measurements were also for a 4 mm pupil.
measuredPupilDiameterMM = 4.0;
calcPupilDiameterMM = 4.0;

% Set up focal length and conversion factors
focalLengthMM = 4.35;
focalLengthMeters = focalLengthMM / 1000;
posteriorNodalDistanceMM = focalLengthMM;
micronsPerDegree = posteriorNodalDistanceMM * 1000 * tand(1);

% Some spatial parameters for visualization and calculation.
spatialsamples = 801;
visualizedSpatialSfrequencyCPD = 10.0;

% Pixels per minute for the PSF. Need to make this
% large enough to capture PSF support, but small enough
% to model variations in the PSF.  A bit of plotting and
% hand fussing to chose.
psfSamplesPerMinute = 0.05;

% 840 nm light used in measurements of Sajdak et al 2019 (section 2.2)
% but we will put in best focus coeffs and are just computing
% at some wavelength without adjustment.
%
% The way we are doing this, we want to use measured wavelength to
% be the target wavelength
measuredWavelength = targetWavelength;

% Read in Zernicke coefficients as matrix. The original table in
% in the spreadsheet (and plot in the paper) has an uncorrected
% number for defocus. We separately read in the defocus that
% maximizes contrast from a separate row of the spreadsheet.  This is the
% updated spreadsheet provided by Austin Roorda with the corrected defocus
% coefficient. The hope is that this will allow us to match MTFs.
zCoeffs = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verB.xlsx','Sheet','Aberration Summaries','Range','B7:L71'));
zCoeff4 = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verB.xlsx','Sheet','Aberration Summaries','Range','B97:L97'));

% Optional zeroing of astigmatism coefficients
NOASTIG = false;
if (NOASTIG)
    zCoeffs(:,3) = 0;
    zCoeffs(:,5) = 0;
end

% Optional diffraction limited PSF
DIFFRACTIONLIMITED = false;
if (DIFFRACTIONLIMITED)
    zCoeffs = zeros(size(zCoeffs));
    zCoeff4 = zeros(size(zCoeff4));
end

%% Set up tree shrew lens absorption
lensAbsorbanceFile = 'treeshrewLensAbsorbance.mat';
targetWavelenth = wavelengthSupport;
% TreeShrew lens absorption. Start with the human lens.
theLens = Lens();

% Load TreeShrew lens unit-density
load(lensAbsorbanceFile, 'wavelength', 'data');

if (isempty(targetWavelenth))
    targetWavelenth = wavelength;
    unitDensity = data;
else
    % Interpolate to optics wavelength support
    unitDensity = interp1(wavelength,data,targetWavelenth, 'pchip');
end

%% Allocate space for storing an MTF slice for each shrew
mtfSlice = zeros(1,spatialsamples);

% Define figures to cycle through if plotPsfOtf is turned on
if plotPsfOtf
psfMeshFig = figure;
otfMeshFig = figure;
end

%% Loop over shrews

% Set up Zernikes. These get offset by 1 from the array, as we
% understand the reporting and storage conventions. The first OSA
% coefficient is piston and that's not included. Then we have zero
% for tip and tilt, also effectively not included. First reported
% coefficient is astigmatism, leading to defocus dropping into
% element 5 of our array.
zCoeffs_TreeShrew = zeros(1,size(zCoeffs,1)+1);
zCoeffs_TreeShrew(2:end) = zCoeffs(:,TSindex);
zCoeffs_TreeShrew(5) = zCoeff4(TSindex);

% Create the wavefront object with desired zcoeefs
wvfP = wvfCreate(...
    'spatialsamples', spatialsamples, ...
    'measured wl', measuredWavelength, ...
    'calc wavelengths', wavelengthSupport, ...
    'zcoeffs', zCoeffs_TreeShrew, ...
    'name', sprintf('treeshrew-%d', calcPupilDiameterMM), ...
    'umPerDegree', micronsPerDegree,...
    'customLCA', @treeShrewLCA);
defocus0 = wvfGet(wvfP,'zcoeffs','defocus');
if (defocus0 ~= zCoeff4(TSindex))
    error('Stored defocus not as intended')
end

% Other wvf properties
wvfP = wvfSet(wvfP, 'measured pupil size', measuredPupilDiameterMM);
wvfP = wvfSet(wvfP, 'calc pupil size', calcPupilDiameterMM);
wvfP = wvfSet(wvfP, 'ref psf sample interval',psfSamplesPerMinute);

% Get PSF angular samples in minutes
wvfPsfAngularSamplesMinutes = wvfGet(wvfP,'psf angular samples');
[wvfXGridMinutes,wvfYGridMinutes] = meshgrid(wvfPsfAngularSamplesMinutes,wvfPsfAngularSamplesMinutes);

% Create optics structure with desired wvf object, and set other optics
% properties.
wvfP = wvfComputePSF(wvfP);
optics = oiGet(wvf2oi(wvfP), 'optics');
optics.name = 'wvf-based tree-shrew optics';
optics = opticsSet(optics, 'focalLength', focalLengthMeters);
optics = opticsSet(optics, 'fnumber', focalLengthMeters*1000/calcPupilDiameterMM);

% Off axis method. We may want to defeat this eventually,
% but it should not affect the PSF itself, nor the computation
% of MTF from PSF.
optics = opticsSet(optics, 'offAxisMethod', 'cos4th');

% Turn pixel vignetting off
optics.vignetting =  0;

% Set up oi
oi.type = 'opticalimage';
oi = oiSet(oi, 'optics', optics);
oi = oiSet(oi, 'name', 'treeshrew');

oi = oiSet(oi, 'bit depth', 32);
oi = oiSet(oi, 'diffuser method', 'skip');
oi = oiSet(oi, 'consistency', 1);

% Get PSF slice at target wavelength
wavePSF = opticsGet(optics,'psf data',targetWavelength);

% Extract PSF support in arcmin and plot a mesh of the PSF
psfSupportMicrons = opticsGet(optics,'psf support','um');
xGridMinutes = 60*psfSupportMicrons{1}/micronsPerDegree;
yGridMinutes = 60*psfSupportMicrons{2}/micronsPerDegree;
if (max(abs(wvfXGridMinutes(:)-xGridMinutes(:))) > 1e-6 | max(abs(wvfYGridMinutes(:)-yGridMinutes(:))) > 1e-6)
    error('Do not get same psf samples from wvf and optics objects');
end

xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
[xSfGridCyclesDeg,ySfGridCyclesDeg,otfraw] = PsfToOtf(xGridMinutes,yGridMinutes,wavePSF);
otf = psfCircularlyAverage(abs(otfraw));

waveMTF = abs(otf);
[~,idx] = min(abs(ySfCyclesDeg));
mtfSlice = waveMTF(idx,:);

if plotPsfOtf
figure(psfMeshFig); clf;
mesh(xGridMinutes,yGridMinutes,wavePSF);
figure(otfMeshFig); clf;
subplot(1,2,1);
mesh(xSfCyclesDeg,ySfCyclesDeg,abs(otfraw));
subplot(1,2,2);
mesh(xSfCyclesDeg,ySfCyclesDeg,abs(otf));
zlim([0 1]);
end

%% Get comparison data

% Several options in spreadsheet and Austin there might be some
% uncertainty about which is the right comparison. Choose one,
% read, and use.
whichCompare = 'maxStrehlWAstig';
switch (whichCompare)
    case 'maxContrastWAstig'
        whichWorksheet = 'MTFs at max Contrast w Astig';
    case 'maxStrehlWAstig'
        whichWorksheet = 'MTFs at max Strehl w Astig';
    case 'maxStrehlNoAstig'
        whichWorksheet = 'MTFs at max Strehl no Astig';
    otherwise
        error('Unknown comparison option chosen');
end
MTFcols = 'B':'L';
extraData.sf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet',whichWorksheet,'Range','A4:A15'));
rowstart = [MTFcols(TSindex(1)),'4'];
rowend = [MTFcols(TSindex(end)),'15'];
extraData.csf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet',whichWorksheet,'Range',[rowstart,':',rowend]));

%% Plot MTF for each tree shrew and wavelength
if plotMTFCompare
figure; clf;
set(gcf,'Position',[1224 460 1126 1129]);
hold on;
plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 10);
set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
ylabel('modulation');
xlim([0 20]);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
ylabel('\it MTF')
MTFcols = 'B':'L';
title({ sprintf('TS# = %d, wl = %d nm, pupli %d mm',TSindex,targetWavelength,calcPupilDiameterMM) ; ...
    sprintf('Opt = %s, comp = %s',whichFigureOfMerit,whichCompare) ...
    });

% Add to plot
extraDataColors = [1 0 0];
plot(extraData.sf, extraData.csf, ...
    'rs-', 'MarkerSize', 9, ...
    'MarkerEdgeColor', squeeze(extraDataColors),  ...
    'MarkerFaceColor', squeeze(extraDataColors)*0.5+[0.5 0.5 0.5]);
ylim([0 1.05]);
ylabel('MTF');
end

%% Create presentation display
% Generate object using display measurements from match-to-sample
% experiment.
presentationDisplay = generateTreeShrewDisplayObject_eemParams();

% Set mean luminance and FOV
meanLuminanceCdPerM2 = 20;
sceneFOVdegs = 5.0;

% Generate default treeshrew cone mosaic with a 5 msec integration time.
% Set microns/deg - might need to change this?
umPerDegree = 76;
imageProps = [1 1];
integrationTime = 1/1000;
nTrialsNum = 10;
sizeConeExcit = 564;

% Set which figures should be plotted for the scene
displayContrastProfiles = false;
displayRadianceMaps = false;
displayRetinalContrastProfiles = false;

% Initialize figures for generated scene, optical images, and cone
% responses
figure;
genScene = tiledlayout('flow','TileSpacing','Compact');
figure;
opticalImage = tiledlayout('flow','TileSpacing','Compact');
% coneResponses = tiledlayout('flow','TileSpacing','Compact');

% Loop through any set of targets to display generated scene, optical
% image, and cone responses
targetidx = 1:6;

% Initialize coneExcitations. Not sure where the size of coneExcitations is
% determined to initilize the size, but the 564 can be replaced.
coneExcitations = zeros(length(targetidx),nTrialsNum,sizeConeExcit,sizeConeExcit);
for ti = 1:length(targetidx)
    
% Load the camels (targets) and realized it on presentation display. 
fileName = ['C:/Users/eemeyer/Documents/GitHub/arcaro_rotation/imagematch_test/camels/target_',num2str(targetidx(ti)-1),'.png'];
realizedStimulusScene = sceneFromFile(fileName, 'monochrome', ...
    meanLuminanceCdPerM2, presentationDisplay);

% Set FOV in scene structure
realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', sceneFOVdegs);

% Visualize generated scene and computed optical image
% Visualize different aspects of the 
tempax = nexttile(genScene);
visualizeScene(realizedStimulusScene, 'displayContrastProfiles', displayContrastProfiles,...
    'displayRadianceMaps',displayRadianceMaps,'axesHandle',tempax);

% Compute the retinal image
oi = oiCompute(oi, realizedStimulusScene);

% Visualize different aspects of the computed optical image
% nexttile(opticalImage);
tempax = nexttile(opticalImage);
visualizeOpticalImage(oi, 'displayRetinalContrastProfiles', displayRetinalContrastProfiles,...
    'displayRadianceMaps',displayRadianceMaps,'axesHandle',tempax);

% create cone mosaic structure
theConeMosaic = coneMosaicTreeShrewCreate(...
    umPerDegree, ...
    'fovDegs', sceneFOVdegs*imageProps, ... % double check these values are correct
    'integrationTime', integrationTime);
% theConeMosaic.displayInfo();
% theConeMosaic.visualizeGrid('backgroundColor', [1 1 1], 'generateNewFigure', true); 

% Generate cone responses and visualize
emPath = zeros(nTrialsNum, 1, 2);
% Compute mosaic excitation responses
coneExcitations(ti,:,:,:) = theConeMosaic.compute(oi, 'emPath', emPath);
size(coneExcitations);

% Visualize cone mosaic and its cone excitation responses
visualizeConeMosaicResponses(theConeMosaic, squeeze(coneExcitations(ti,:,:,:)), 'R*/cone/tau');

end
