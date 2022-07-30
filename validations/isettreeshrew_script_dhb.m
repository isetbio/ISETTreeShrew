% Replicate as best we can the results in
% Sajdak et al., 2019, Experimental Eye Research,
% 185, 107683.
%
% Thanks to Austin Roorda for clarifications and for
% providing tabulated data and example calculations.
%
% In this version, we unpack some of the encapsulated
% tree shrew routines to expose more parameters.  See 
% oiTreeShrewCreate for the packed up version.

% History:
%   07/22/2022  dhb  Move from original Nicolas/Emily live script versions.

% Configuration.  Execute the tbUseProjec command
% below to configure, if you use ToolboxToolbox.
%{
    tbUseProject('ISETTreeShrew');
%}

%% Initialize
clear; close all;

%% Parameters
%
% Specify wavelengths to compute on, and which individual
% tree shrews to analyze here.
targetWavelength = 550;

% Normally this is a range of wavelengths, but we are only
% trying to make sense of results at one wavelength here,
% so just do one for speed.
wavelengthSupport = targetWavelength;

% Special flags
DIFFRACTIONLIMITED = false;
ZERODEFOCUS = false;
NOASTIG = false;
ROORDA_COMPARE = false;
ROORDA_SAMPLING = false;
FLIP_DEFOCUS_SIGN = false;

% Can specify any contiguous range between 1 and 11.
%
% There were 11 shrews measured, and we have from Roorda the
% tabulated Zernike coefficients in an Excel spreadsheet.
TSindex = 1:11;

% Define files for direct comparisons with Roorda calculations
if (ROORDA_COMPARE)
    if (DIFFRACTIONLIMITED)
        roordaWvfFile = 'pr_0D_defocus_4mm_pupil_WF.csv';
        roordaPSFFile = 'pr_0D_defocus_4mm_pupil_PSF.csv';
        roordaMTFFile = 'pr_0D_defocus_4mm_pupil_MTF.csv';
    else
        switch (TSindex)
            case 1
                roordaWvfFile = '216OS_0D_defocus_4mm_pupil_WF.csv';
                roordaPSFFile = '216OS_0D_defocus_4mm_pupil_PSF.csv';
                roordaMTFFile = '216OS_0D_defocus_4mm_pupil_fullMTF.csv';
            otherwise
                error('Invalid tree shrew index for comparing with Roorda calculations');
        end
    end
end

% Figure of merit to optimize.
% Choices are:
%   'None'
%   'PSFPeak';
%   'MTFArea';
whichFigureOfMerit = 'None';
if (~strcmp(whichFigureOfMerit,'None'))
    PSFPeakFig = figure; clf;
    set(gcf,'Position',[1224 460 1126 1129]);
    MTFAreaFig = figure; clf;
    set(gcf,'Position',[1224 460 1126 1129]);
end

% Pupil diameter
%
% We're trying to reproduce Figure 2, which is for a 4 mm pupil. So we
% calculate for that. Measurements were also for a 4 mm pupil.
measuredPupilDiameterMM = 4.0;
calcPupilDiameterMM = 4.0;

% Set up focal length and conversion factors
focalLengthMM = 4.35;
focalLengthMeters = focalLengthMM / 1000;
posteriorNodalDistanceMM = focalLengthMM;
micronsPerDegree = posteriorNodalDistanceMM * 1000 * tand(1);

% Some spatial parameters for visualization and calculation.
%
% Pixels per minute for the PSF. Need to make this large enough to capture
% PSF support, but small enough to model variations in the PSF.  A bit of
% plotting and hand fussing to choose.
roordaPSFRangeArcMin = 3600/60;
if (ROORDA_SAMPLING)
    spatialsamples = 513;
    PSFMinutesPerSample = 2*roordaPSFRangeArcMin/spatialsamples;
else
    spatialsamples = 801;
    PSFMinutesPerSample = 0.05;
end
visualizedSpatialSfrequencyCPD = 10.0;

% 840 nm light used in measurements of Sajdak et al 2019 (section 2.2)
% but we will put in best focus coeffs and are just computing
% at some wavelength without adjustment.
%
% The way we are doing this, we want to use measured wavelength to
% be the target wavelength
measuredWavelength = targetWavelength;

%% Data 
%
% Get paths to data we need
rootDir = isetTreeShrewRootPath;
dataDir = fullfile(rootDir,'data','TabulatedDataFromRoorda');
zernikeFile = 'Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx';
zernikePath = fullfile(dataDir,zernikeFile);

%% Read in Zernicke coefficients as matrix.
% 
% The original table in in the spreadsheet (and plot in the paper) has an uncorrected
% number for defocus. We separately read in the defocus that
% maximizes contrast from a separate row of the spreadsheet.  This
% was provided by Austin Roorda and the hope is this would allow
% us to match MTFs.
zCoeffs = cell2mat(readcell(zernikePath,'Sheet','Aberration Summaries','Range','B7:L71'));
zCoeff4 = cell2mat(readcell(zernikePath,'Sheet','Aberration Summaries','Range','B95:L95'));

%% Optional check of sign convention on defocus term
if (FLIP_DEFOCUS_SIGN)
    zCoeff4 = -zCoeff4;
end

%% Optional zero defocus
if (ZERODEFOCUS)
    zCoeff4 = zeros(size(zCoeff4));
end

%% Optional zeroing of astigmatism coefficients
%
% Might do this
if (NOASTIG)
    zCoeffs(:,3) = 0;
    zCoeffs(:,5) = 0;
end

%% Optional diffraction limited PSF
if (DIFFRACTIONLIMITED)
    zCoeffs = zeros(size(zCoeffs));
    zCoeff4 = zeros(size(zCoeff4));
end

%% Get Roorda comparsion data if needed
if (ROORDA_COMPARE)
    roordaWvf = csvread(fullfile(dataDir,roordaWvfFile));
    roordaPSF = csvread(fullfile(dataDir,roordaPSFFile));
    roordaPSF = roordaPSF/sum(sum(roordaPSF(:)));
    roordaMTF = csvread(fullfile(dataDir,roordaMTFFile));
end

%% Compute Roorda's MTF from his PSF
%
% Spatial samples, putting 0 at the right place to match
% where we think the center of Roorda's PSF is.
if (ROORDA_COMPARE)
    [m,n] = size(roordaPSF);
    roordaPSFSamplesPos = roordaPSFRangeArcMin*linspace(0,1,m/2);
    diffSamples = diff(roordaPSFSamplesPos);
    roordaPSFSamplesArcMin = -roordaPSFRangeArcMin-diffSamples:diffSamples:roordaPSFRangeArcMin;
    [roordaXGridMinutes,roordaYGridMinutes] = meshgrid(roordaPSFSamplesArcMin,roordaPSFSamplesArcMin);

    % Find location of PSF max.  For diffraction limited, this should be at
    % m/2+1.
    maxPSF = -Inf;
    for ii = 1:m
        for jj=1:n
            if (roordaPSF(ii,jj) > maxPSF)
                maxI = ii;
                maxJ = jj;
                maxPSF = roordaPSF(ii,jj);
            end
        end
    end
    fprintf('Maximum PSF at %d, %d\n',maxI,maxJ);
    fprintf('Sf at PSF max, %f\n',roordaPSFSamplesArcMin(maxI));

    % Find OTF/MTF from Roorda PSF
    [ourRoordaXSfGridCyclesDeg,ourRoordaYSfGridCyclesDeg,ourRoordaOTF] = PsfToOtf(roordaXGridMinutes,roordaYGridMinutes,roordaPSF);
    ourRoordaMTF = abs(ourRoordaOTF);
    ourRoordaMTFCirc = psfCircularlyAverage(ourRoordaMTF);
    ourRoordaMTFSlice = ourRoordaMTFCirc(m/2+1,m/2+1:end);
    ourRoordaMTFSliceSfsCyclesDeg = ourRoordaXSfGridCyclesDeg(m/2+1,m/2+1:end);
    roordaMTFCirc = psfCircularlyAverage(roordaMTF);

    % Plot PSF and MTF slices
    %
    % Roorda PSF and MTF in red
    % MTF computed here from PSF in black
    roordaFig = figure; clf;
    set(gcf,'Position',[120 747 1794 682]);
    subplot(1,2,1); hold on;
    plot(roordaPSFSamplesArcMin,roordaPSF(m/2+1,:),'ro','MarkerFaceColor','r','MarkerSize',10);
    plot(roordaPSFSamplesArcMin,roordaPSF(m/2+1,:),'r','LineWidth',5);
    if (DIFFRACTIONLIMITED)
        sliceXlim = 3;
    else
        sliceXlim = 10;
    end
    xlim([-sliceXlim sliceXlim]);
    xlabel('Postion (arcmin)'); ylabel('PSF');
    subplot(1,2,2); hold on;
    plot(ourRoordaMTFSliceSfsCyclesDeg,roordaMTFCirc(m/2+1,m/2+1:end),'r','LineWidth',5);
    plot(ourRoordaMTFSliceSfsCyclesDeg,roordaMTFCirc(m/2+1,m/2+1:end),'ro','MarkerFaceColor','r','MarkerSize',10);
    plot(ourRoordaMTFSliceSfsCyclesDeg,ourRoordaMTFSlice,'k','LineWidth',4);
    xlabel('Spatial Freq (c/deg)'); ylabel('Circ Avg MTF');
    drawnow;
end

%% Set up tree shrew lens absorption
lensAbsorbanceFile = 'treeshrewLensAbsorbance.mat';
targetWavelenth = wavelengthSupport;

%% Need to put in lens density for real lens

%% Allocate space for storing an MTF slice for each shrew
MTFSlice = zeros(length(TSindex),spatialsamples);

% Define figures to cycle through
PSFMeshFig = figure;
OTFMeshFig = figure;

%% Loop over shrews
for ts = 1:length(TSindex)

    % Set up Zernikes. These get offset by 1 from the array, as we
    % understand the reporting and storage conventions. The first OSA
    % coefficient is piston and that's not included. Then we have zero
    % for tip and tilt, also effectively not included. First reported
    % coefficient is astigmatism, leading to defocus dropping into
    % element 5 of our array.
    zCoeffs_TreeShrew = zeros(1,size(zCoeffs,1)+1);
    zCoeffs_TreeShrew(2:end) = zCoeffs(:,TSindex(ts));
    zCoeffs_TreeShrew(5) = zCoeff4(TSindex(ts));

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
    if (defocus0 ~= zCoeff4(TSindex(ts)))
        error('Stored defocus not as intended')
    end

    % Other wvf properties
    wvfP = wvfSet(wvfP, 'measured pupil size', measuredPupilDiameterMM);
    wvfP = wvfSet(wvfP, 'calc pupil size', calcPupilDiameterMM);
    wvfP = wvfSet(wvfP, 'ref PSF sample interval',PSFMinutesPerSample);

    % Get PSF angular samples in minutes
    wvfPSFAngularSamplesMinutes = wvfGet(wvfP,'PSF angular samples');
    [wvfXGridMinutes,wvfYGridMinutes] = meshgrid(wvfPSFAngularSamplesMinutes,wvfPSFAngularSamplesMinutes);

    % Loop over a list of defocus values and find one that maximizes the
    % Streh ratio. This is the same as maximizing the peak of the PSF.
    % Also find area under MTF, which is another figure of merit we could
    % maximize.
    if (~strcmp(whichFigureOfMerit,'None'))
        minDefocusDelta = -1.5;
        maxDefoucusDelta = 1.5;
        nDefocusDeltas = 40;
        defocusDeltas = linspace(minDefocusDelta,maxDefoucusDelta,nDefocusDeltas);
        defocusZCoeffs = zCoeffs_TreeShrew;
        for dd = 1:length(defocusDeltas)
            % Set adjusted defocus
            defocusZ1 = defocusZCoeffs;
            defocusZ1(5) = defocus0+defocusDeltas(dd);
            wvfP1 = wvfSet(wvfP,'zcoeffs',defocusZ1);

            % Compute PSF, find and store peak
            wvfP1 = wvfComputePSF(wvfP1);
            PSF1 = wvfGet(wvfP1,'psf',targetWavelength);
            peakPSF1(dd) = max(PSF1(:));

            % Compute MTF, find and store area under it
            [xSfGridCyclesDeg,ySfGridCyclesDeg,OTFraw] = PsfToOtf(wvfXGridMinutes,wvfYGridMinutes,PSF1);
            MTF = psfCircularlyAverage(abs(OTFraw));
            MTFArea1(dd) = sum(MTF(:));
        end

        % Plot peak of PSF and MTF area versus defocus
        figure(PSFPeakFig); 
        subplot(4,3,ts); hold on
        plot(defocusDeltas,peakPSF1,'ro','MarkerFaceColor','r','MarkerSize',12);
        xlabel('Defocus Delta'); ylabel('Peak PSF');
        title(sprintf('TS = %d, wavelength = %s',TSindex(ts),num2str(targetWavelength)));
        figure(MTFAreaFig); 
        subplot(4,3,ts); hold on
        plot(defocusDeltas,MTFArea1,'ro','MarkerFaceColor','r','MarkerSize',12);
        xlabel('Defocus Delta'); ylabel('MTF area');
        title(sprintf('TS = %d, wavelength = %s',TSindex(ts),num2str(targetWavelength)));
        drawnow;

        % Choose and set optimal defocus based on calcs above
        switch whichFigureOfMerit
            case 'PSFPeak'
                [~,idx] = max(peakPSF1);
                bestDefocus(ts) = defocus0+defocusDeltas(idx);
            case 'MTFArea'
                [~,idx] = max(MTFArea1);
                bestDefocus(ts) = defocus0+defocusDeltas(idx);
            case 'None'
                bestDefocus(ts) = defocus0;
            otherwise
                error('Unknown figure of merit chosen');
        end
        defocusZ1(5) = bestDefocus(ts);
        wvfP = wvfSet(wvfP,'zcoeffs',defocusZ1);
    end

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
    PSFSupportMicrons = opticsGet(optics,'psf support','um');
    xGridMinutes = 60*PSFSupportMicrons{1}/micronsPerDegree;
    yGridMinutes = 60*PSFSupportMicrons{2}/micronsPerDegree;
    if (max(abs(wvfXGridMinutes(:)-xGridMinutes(:))) > 1e-6 | max(abs(wvfYGridMinutes(:)-yGridMinutes(:))) > 1e-6)
        error('Do not get same PSF samples from wvf and optics objects');
    end
    figure(PSFMeshFig); clf;
    mesh(xGridMinutes,yGridMinutes,wavePSF);
    xlabel('x (minutes)'); ylabel('y (minutes)'); zlabel('PSF');
    axis('square');

    if (ROORDA_COMPARE)
        figure(roordaFig);
        subplot(1,2,1); hold on
        plot(xGridMinutes(floor(spatialsamples/2)+1,:),max(roordaPSF(m/2+1,:))*wavePSF(floor(spatialsamples/2)+1,:)/max(wavePSF(floor(spatialsamples/2)+1,:)), 'g-', 'LineWidth', 2);
        xlim([-sliceXlim sliceXlim]);
    end

    xSfCyclesDeg = opticsGet(optics, 'OTF fx', 'cyclesperdeg');
    ySfCyclesDeg = opticsGet(optics, 'OTF fy', 'cyclesperdeg');
    [xSfGridCyclesDeg,ySfGridCyclesDeg,OTFraw] = PsfToOtf(xGridMinutes,yGridMinutes,wavePSF);
    figure(OTFMeshFig); clf;
    subplot(1,2,1);
    mesh(xSfCyclesDeg,ySfCyclesDeg,abs(OTFraw));
    xlabel('x sf (cyc/deg'); ylabel('y sf (cyc/deg)'); zlabel('MTF');
    axis('square');
    zlim([0 1]);
    MTF = psfCircularlyAverage(abs(OTFraw));
    subplot(1,2,2);
    mesh(xSfCyclesDeg,ySfCyclesDeg,abs(MTF));
    xlabel('x sf (cyc/deg'); ylabel('y sf (cyc/deg)'); zlabel('Circ Avg MTF');
    axis('square');
    zlim([0 1]);

    [~,idx] = min(abs(ySfCyclesDeg));
    MTFSlice(ts,:) = MTF(idx,:);

end

%% Get comparison data
%
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
extraData.sf = cell2mat(readcell(zernikePath,'Sheet','MTFs at max Strehl w Astig','Range','A4:A15'));
rowstart = [MTFcols(TSindex(1)),'4'];
rowend = [MTFcols(TSindex(end)),'15'];
extraData.csf = cell2mat(readcell(zernikePath,'Sheet','MTFs at max Strehl w Astig','Range',[rowstart,':',rowend]));

%% Plot MTF for each tree shrew and wavelength
figure; clf;
set(gcf,'Position',[1224 460 1126 1129]);
for ts = 1:length(TSindex)
    subplot(4,3,ts); hold on;
    plot(xSfCyclesDeg, MTFSlice(ts,:), 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 10);
    set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
    ylabel('modulation');
    xlim([0 20]);
    xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
    ylabel('\it MTF')
    MTFcols = 'B':'L';
    title({ sprintf('TS# = %d, wl = %d nm, pupli %d mm',TSindex(ts),targetWavelength,calcPupilDiameterMM) ; ...
        sprintf('Opt = %s, comp = %s',whichFigureOfMerit,whichCompare) ...
        });

    % Add to plot
    extraDataColors = [1 0 0];
    plot(extraData.sf, extraData.csf(:,ts), ...
        'rs-', 'MarkerSize', 9, ...
        'MarkerEdgeColor', squeeze(extraDataColors),  ...
        'MarkerFaceColor', squeeze(extraDataColors)*0.5+[0.5 0.5 0.5]);
    ylim([0 1.05]);
    ylabel('MTF');

    % Add our MTF to figure of what Roorda sent
    if (ROORDA_COMPARE)
        figure(roordaFig);
        subplot(1,2,2); hold on
        plot(xSfCyclesDeg, MTFSlice(ts,:), 'g-', 'LineWidth', 2);
        if (DIFFRACTIONLIMITED)
            xlim([0 130]);
        else
            xlim([0 80]);
        end
    end
end

function lcaDiopters = treeShrewLCA(wl1NM, wl2NM)
% We dont have a model LCA for tree shrew yet.
% Here we model is as the human LCA x 5
% This creates an LCA difference between 840nm and 550nm of -4.5D
% as per Sadjak et al (2019), "Noninvasive imaging of the tree shrew eye:
% Wavefront analysis and retinal imaging with correlative histology"

constant = 1.8859 - (0.63346 ./ (0.001 .* wl1NM - 0.2141));
lcaDiopters = 1.8859 - constant - (0.63346 ./ (0.001 * wl2NM - 0.2141));
lcaDiopters = 5.15 * lcaDiopters;
end