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
% tabulated Zernike coefficients in an Excel spreadsheet. Will likely
% change this once we determine an ideal default shrew.
TSindex = 1;

% Can turn on/off which figures to plot
plotFigureOfMerit = false;
plotPsfOtf = false;
plotMTFCompare = true;

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
zCoeffs = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verC.xlsx','Sheet','Aberration Summaries','Range','B7:L71'));
zCoeff4 = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verC.xlsx','Sheet','Aberration Summaries','Range','B97:L97'));

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
% This section requires iterating over multiple wavelengths for interp1 to
% work. Will leave this commented for now.


% lensAbsorbanceFile = 'treeshrewLensAbsorbance.mat';
% targetWavelenth = wavelengthSupport;
% 
% % Load TreeShrew lens unit-density
% load(lensAbsorbanceFile, 'wavelength', 'data');
% 
% % TreeShrew lens absorption. Start with the human lens.
% theLens = Lens();
% 
% if (isempty(targetWavelenth))
%     targetWavelenth = wavelength;
%     unitDensity = data;
% else
%     % Interpolate to optics wavelength support
%     unitDensity = interp1(wavelength,data,targetWavelenth, 'pchip');
% end

% Update the lens
% set(theLens,'wave', targetWavelenth);
% set(theLens,'unitDensity',unitDensity);

% from oiTreeShrewCreate
% Update the oi.lens
% oi = oiSet(oi, 'lens', theLens);

% Set the same lens in the optics structure too
% oi.optics.lens = theLens;

%% Allocate space for storing an MTF slice for each shrew
mtfSlice = zeros(1,spatialsamples);

% Define figures to display if plotPsfOtf is turned on
if plotPsfOtf
psfMeshFig = figure;
otfMeshFig = figure;
end

%% Set up OI and extract PSF and OTF data from single shrew

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

% Get x,y degrees for OTF plotting, convert PSF to OTF, and circularly
% average
xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
[xSfGridCyclesDeg,ySfGridCyclesDeg,otfraw] = PsfToOtf(xGridMinutes,yGridMinutes,wavePSF);
otf = psfCircularlyAverage(abs(otfraw));

% Extract single MTF slice at y minimum
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
extraData.sf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verC.xlsx','Sheet',whichWorksheet,'Range','A4:A15'));
rowstart = [MTFcols(TSindex(1)),'4'];
rowend = [MTFcols(TSindex(end)),'15'];
extraData.csf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018_verC.xlsx','Sheet',whichWorksheet,'Range',[rowstart,':',rowend]));

%% Plot MTF for each tree shrew and wavelength
if plotMTFCompare
figure; clf;
set(gcf,'Position',[1000 580 540 300]);
hold on;
plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 10);
set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
ylabel('modulation');
xlim([0 20]);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
ylabel('\it MTF')
MTFcols = 'B':'L';
title({sprintf('TS# = %d, wl = %d nm, pupli %d mm',TSindex,targetWavelength,calcPupilDiameterMM)});

% Add to plot
extraDataColors = [1 0 0];
plot(extraData.sf, extraData.csf, ...
    'rs-', 'MarkerSize', 9, ...
    'MarkerEdgeColor', squeeze(extraDataColors),  ...
    'MarkerFaceColor', squeeze(extraDataColors)*0.5+[0.5 0.5 0.5]);
ylim([0 1.05]);
ylabel('MTF');
end

%% Compare this with default from ISETTreeShrew

theWVFOI = oiTreeShrewCreate('opticsType', 'wvf', ...
    'name', 'wvf-based optics');

% Get PSF slice at target wavelength
optics2 = oiGet(theWVFOI,'optics');

psfSupportMicrons2 = opticsGet(optics2,'psf support','um');
xGridMinutes2 = 60*psfSupportMicrons2{1}/micronsPerDegree;
yGridMinutes2 = 60*psfSupportMicrons2{2}/micronsPerDegree;

wavePSF2 = opticsGet(optics2,'psf data',targetWavelength);
xSfCyclesDeg2 = opticsGet(optics2, 'otf fx', 'cyclesperdeg');
ySfCyclesDeg2 = opticsGet(optics2, 'otf fy', 'cyclesperdeg');
[xSfGridCyclesDeg2,ySfGridCyclesDeg2,otfraw2] = PsfToOtf(xGridMinutes2,yGridMinutes2,wavePSF2);
otf2 = psfCircularlyAverage(abs(otfraw2));

% Extract single MTF slice at y minimum
waveMTF2 = abs(otf2);
[~,idx2] = min(abs(ySfCyclesDeg));
mtfSlice2 = waveMTF2(idx2,:);

figure; clf;
set(gcf,'Position',[1000 580 540 300]);
plot(xSfCyclesDeg2, mtfSlice2, 'go-', 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 10);

hold on;
plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 5);
set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
ylabel('modulation');
xlim([0 20]);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
ylabel('\it MTF')
MTFcols = 'B':'L';
title({sprintf('TS# = %d, wl = %d nm, pupli %d mm',TSindex,targetWavelength,calcPupilDiameterMM)});


%% Create presentation display with Kell images
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
nTrialsNum = 100;
sizeConeExcit = 564;

% Define any set of targets and distractors with imgidx to display  
% generated scene, optical image, and cone responses
folderName = 'C:/Users/eemeyer/Documents/GitHub/arcaro_rotation/imagematch_test/';
imgFolder{1} = dir([folderName,'camels/target*']);
imgFolder{2} = dir([folderName,'wrenches/distractor*']);
imgidx = 1; %1:length(targetFile)

% Set which figures should be plotted for the scene
displayContrastProfiles = false;
displayRadianceMaps = false;
displayRetinalContrastProfiles = false;
displayConeMosaic = false;
displayImages = true;

% Initialize figures for generated scene and optical images
if displayImages
    figure;
    genScene = tiledlayout('flow','TileSpacing','Compact');
    figure;
    opticalImage = tiledlayout('flow','TileSpacing','Compact');
end

% Initialize coneExcitations. Not sure where the size of coneExcitations is
% determined to initilize the size, but sizeConeExcit can be replaced.
coneExcitations = {};
coneExcitations{1} = zeros(length(imgidx),nTrialsNum,sizeConeExcit,sizeConeExcit);
coneExcitations{2} = zeros(length(imgidx),nTrialsNum,sizeConeExcit,sizeConeExcit);

% Loop through the number of images and loop through the targets and
% distractors (itype). Resulting cone excitations for all images are stored
% in the cell array coneExcitations with cell 1: targets, cell 2:
% distractors
for img = 1:length(imgidx)
    for itype = 1:2
        % Load the camels (targets) and realized it on presentation display. 
        imgFile = fullfile(imgFolder{itype}(img).folder,imgFolder{itype}(img).name);
        realizedStimulusScene = sceneFromFile(imgFile, 'monochrome', ...
            meanLuminanceCdPerM2, presentationDisplay);

        % Set FOV in scene structure
        realizedStimulusScene = sceneSet(realizedStimulusScene, 'fov', sceneFOVdegs);

        % Compute the retinal image
        oi = oiCompute(oi, realizedStimulusScene);
        
        if displayImages
        % Visualize generated scene and computed optical image
        % Visualize different aspects of the 
        ax_temp = nexttile(genScene);
        visualizeScene(realizedStimulusScene, 'displayContrastProfiles', displayContrastProfiles,...
            'displayRadianceMaps',displayRadianceMaps,'axesHandle',ax_temp);

        % Visualize different aspects of the computed optical image
        % nexttile(opticalImage);
        ax_temp = nexttile(opticalImage);
        visualizeOpticalImage(oi, 'displayRetinalContrastProfiles', displayRetinalContrastProfiles,...
            'displayRadianceMaps',displayRadianceMaps,'axesHandle',ax_temp);
        end
        
        % create cone mosaic structure
        theConeMosaic = coneMosaicTreeShrewCreate(...
            umPerDegree, ...
            'fovDegs', sceneFOVdegs*imageProps, ... % double check these values are correct
            'integrationTime', integrationTime);

        % Generate cone responses and visualize
        emPath = zeros(nTrialsNum, 1, 2);
        % Compute mosaic excitation responses
        coneExcitations{itype}(img,:,:,:) = theConeMosaic.compute(oi, 'emPath', emPath);

        % Visualize cone mosaic and its cone excitation responses
        if displayConeMosaic
            theConeMosaic.displayInfo();
            theConeMosaic.visualizeGrid('backgroundColor', [1 1 1], 'generateNewFigure', true); 

            visualizeConeMosaicResponses(theConeMosaic, squeeze(coneExcitations{itype}(img,:,:,:)), 'R*/cone/tau');
        end
    end
end

%% Attempt at SVM, code from ls_inferenceBinarySVM in ISETLiveScripts

% Simulate a 2AFC task 
taskIntervals = 1;
[classificationMatrix, classLabels] = generateSetUpForClassifier(theConeMosaic, ...
    squeeze(coneExcitations{1}), squeeze(coneExcitations{2}), taskIntervals);

% Find principal components of the responses
[pcVectors, ~, ~, ~,varianceExplained] = pca(classificationMatrix);

% Project the responses onto the space formed by the first 4 PC vectors
pcComponentsNumForClassification = 2;
classificationMatrixProjection = classificationMatrix * pcVectors(:,1:pcComponentsNumForClassification);

% Visualize the first 4 principal components.
visualizePrincipalComponents(pcVectors, varianceExplained, theConeMosaic);

% Train a binary SVM classifier and visualize the support vectors in 2
% dimensions
svm = fitcsvm(classificationMatrixProjection,classLabels);

% Visualize the data along with the hyperplane computed by the SVM 
visualizeSVMmodel(svm, classificationMatrixProjection, classLabels);
