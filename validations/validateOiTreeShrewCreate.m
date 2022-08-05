% ValidateOITreeShrewCreate
%
% Description:
%     Demonstrate the unpacked create of the tree shrew optics object does
%     the same thing as oiTreeShrewCreate.
%
% See Also: validateAgainstRoordaCalcs

% History:
%   07/22/2022  dhb  Move from original Nicolas/Emily live script versions.
%   07/26/2022  eem  Move from David script to simplify/add cone mosaic.
%   08/06/2022  dhb, eem Work towards validating oiTreeShrewCreate only,
%                    leaving other things in other scripts.

% Configuration.
%
% Execute the tbUseProjec command
% below to configure, if you use ToolboxToolbox.
%{
    tbUseProject('ISETTreeShrew');
%}

%% Initialize
clear; close all;

%% Things should match up across choices of these parameters.
%
% These correspond to key/value pairs that can be passed through
% oiTreeShrewCreate to opticsTreeShrewCreate.
%
% See help opticsTreeShrewCreate for information
TSindex = 2;                        % 'whichShrew'
spatialsamples = 501;               % 'spatialSamples'
psfSamplesPerMinute = 0.06;         % 'psfSamplesPerMinute'
wavelengthSupport = 410:10:690;     % 'wavelengthSupport'
measuredWavelength = 550;           % 'measuredWavelength'
calcPupilDiameterMM = 3.0;          % 'pupilDiameterMM'
focalLengthMM = 4.35;               % 'focalLengthMM'

%% Call oiTreeShrewCreate
theWVFOI = oiTreeShrewCreate('opticsType', 'wvf', ...
    'name', 'wvf-based optics', ...
    'whichShrew', TSindex, ...
    'spatialSamples', spatialsamples, ...
    'psfSamplesPerMinute', psfSamplesPerMinute, ...
    'wavelengthSupport', wavelengthSupport, ...
    'measuredWavelength', measuredWavelength, ...
    'pupilDiameterMM',calcPupilDiameterMM, ...
    'focalLengthMM', focalLengthMM);

%% Unpacked oiTreeShrewCreate
%
% Here we explicitly do much of what oiTreeShrewCreate/opticsTreeShrew
% create do, and at the end compare that they match.

% Specify wavelengths to compare on.  Choosing 2 allows us to see if
% the LCA seems sensible
targetWavelengths = [550 460];

% Pupil diameter
%
% We're trying to reproduce Figure 2,
% which is for a 4 mm pupil. So we calculate for that.
% Measurements were also for a 4 mm pupil.
%
% This is so tightly coupled to the Sadjek et al. data
% that we read from a file that we don't allow it to be
% passed to oiTreeShrewCreate.
measuredPupilDiameterMM = 4.0;

% Set up focal length and conversion factors
focalLengthMeters = focalLengthMM / 1000;
posteriorNodalDistanceMM = focalLengthMM;
micronsPerDegree = posteriorNodalDistanceMM * 1000 * tand(1);

% Some spatial parameters for visualization and calculation.
visualizedSpatialSfrequencyCPD = 10.0;

% Read in Zernicke coefficients as matrix. The original table in
% in the spreadsheet (and plot in the paper) has an uncorrected
% number for defocus. We separately read in the defocus that
% maximizes contrast from a separate row of the spreadsheet.  This is the
% updated spreadsheet provided by Austin Roorda with the corrected defocus
% coefficient. The hope is that this will allow us to match MTFs.
rootDir = isetTreeShrewRootPath;
dataDir = fullfile(rootDir,'data','TabulatedDataFromRoorda');
zernikeFile = 'Tree_Shrew_Aberrations_Remeasured_Oct2018_verC.xlsx';
zernikePath = fullfile(dataDir,zernikeFile);
zCoeffs = cell2mat(readcell(zernikePath,'Sheet','Aberration Summaries','Range','B7:L71'));
zCoeff4 = cell2mat(readcell(zernikePath,'Sheet','Aberration Summaries','Range','B97:L97'));

%% Allocate space for storing an MTF slice for each shrew
mtfSlice = zeros(1,spatialsamples);

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

% Get PSF slice at target wavelengths and do comparison for each
for ii = 1:length(targetWavelengths)
    % Set target wavelength
    targetWavelength = targetWavelengths(ii);
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
    ySfCyclesDegFromOI = opticsGet(optics, 'otf fy', 'cyclesperdeg');
    [xSfGridCyclesDegFromOI,ySfGridCyclesDegFromOI,otfraw] = PsfToOtf(xGridMinutes,yGridMinutes,wavePSF);
    otf = psfCircularlyAverage(abs(otfraw));

    % Extract single MTF slice at y minimum
    waveMTF = abs(otf);
    [~,idx] = min(abs(ySfCyclesDegFromOI));
    mtfSlice = waveMTF(idx,:);

    % Get PSF slice at target wavelength along with coordinates for
    % plotting
    opticsFromOI = oiGet(theWVFOI,'optics');
    psfSupportMicronsFromOI = opticsGet(opticsFromOI,'psf support','um');
    xGridMinutesFromOI = 60*psfSupportMicronsFromOI{1}/micronsPerDegree;
    yGridMinutesFromOI = 60*psfSupportMicronsFromOI{2}/micronsPerDegree;

    wavePSFFromOI = opticsGet(opticsFromOI,'psf data',targetWavelength);
    xSfCyclesDegFromOI = opticsGet(opticsFromOI, 'otf fx', 'cyclesperdeg');
    ySfCyclesDegFromOI = opticsGet(opticsFromOI, 'otf fy', 'cyclesperdeg');
    [xSfGridCyclesDegFromOI,ySfGridCyclesDegFromOI,otfrawFromOI] = PsfToOtf(xGridMinutesFromOI,yGridMinutesFromOI,wavePSFFromOI);
    otfFromOI = psfCircularlyAverage(abs(otfrawFromOI));

    % Extract single MTF slice at y minimum
    waveMTFFromOI = abs(otfFromOI);
    [~,idx] = min(abs(ySfCyclesDegFromOI));
    mtfSliceFromOI = waveMTFFromOI(idx,:);

    % Plot MTF slice obtained both ways
    figure; clf;
    set(gcf,'Position',[1000 580 540 300]);
    plot(xSfCyclesDegFromOI, mtfSliceFromOI, 'go-', 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 10);

    hold on;
    plot(xSfCyclesDeg, mtfSlice, 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 5);
    set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
    ylabel('modulation');
    xlim([0 20]);
    xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
    ylabel('\it MTF')
    MTFcols = 'B':'L';
    title({sprintf('TS# = %d, wl = %d nm, pupli %d mm',TSindex,targetWavelength,calcPupilDiameterMM)});
end


