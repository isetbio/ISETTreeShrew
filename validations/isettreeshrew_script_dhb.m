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

% Configuration.  Execute the tbUseProjec command
% below to configure, if you use ToolboxToolbox.
%{
    tbUseProject('ISETTreeShrew');
%}

%% Initialize
clear; close all;

% this code contains the original functions from iset bio used to make
% figures
% [theWVFOI,wvfP2] = oiTreeShrewCreate('opticsType', 'wvf', ...
%     'name', 'wvf-based optics');
%
% psfRangeArcMin = 50;
% visualizedSpatialSfrequencyCPD = 10.0;
% targetWavelength = 500;
% visualizeOptics(theWVFOI, targetWavelength, ...
%     psfRangeArcMin, ...
%     visualizedSpatialSfrequencyCPD, ...
%     'extraOTFData',  mtfTreeShrewFromPaper('SaidakEtAl_2019'));

%% Unpacked oiTreeShrewCreate

% Specify wavelengths to compute on, and which individual
% tree shrews to analyze here.
%
% There were 11 shrews measured, and we have from Roorda the
% tabulated Zernike coefficients in an Excel spreadsheet.
targetWavelength = 840;

% 11 shrews measured.  Can specify a contiguous range.
TSindex = 3; % 11 total shrews

% Wavelength support for calculations
wavelengthSupport = 450:10:900;

% Pupil diameter
%
% We're trying to reproduce Figure 2,
% which is for a 4 mm pupil. So we calculate for that.
% Measurements were also for a 4 mm pupil.
calcPupilDiameterMM = 4.0;
measuredPupilDiameterMM = 4.0;

% Set up focal length and conversion factors
focalLengthMM = 4.35;
focalLengthMeters = focalLengthMM / 1000;
posteriorNodalDistanceMM = focalLengthMM;
micronsPerDegree = posteriorNodalDistanceMM * 1000 * tand(1);

% Some spatial parameters for visualization and calculation.
spatialsamples = 801;
psfRangeArcMin = 20;
visualizedSpatialSfrequencyCPD = 10.0;

% Pixels per minute for the PSF. Need to make this
% large enough to capture PSF support, but small enough
% to model variations in the PSF.  A bit of plotting and
% hand fussing to chose.
psfSamplesPerMinute = 0.05;

% 840 nm light used in measurements of Sajdak et al 2019 (section 2.2)
% but we will put in best focus coeffs and are calling this 550.
% The way we are doing this, we want to use measured wavelength to
% be the wavelength we think the shrews accommodate to.
measuredWavelength = 840;

% Read in Zernicke coefficients as matrix. The original table in
% in the spreadsheet (and plot in the paper) has an uncorrected
% number for defocus. We separately read in the defocus that
% maximizes contrast from a separate row of the spreadsheet.  This
% was provided by Austin Roorda and the hope is this would allow
% us to match MTFs.
zCoeffs = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet','Aberration Summaries','Range','B7:L71'));
zCoeff4 = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet','Aberration Summaries','Range','B95:L95'));

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

% Default lens info
theLens = Lens();

% Load tree shrew lens unit-density and spline to our wavelengths
load(lensAbsorbanceFile, 'wavelength', 'data');
unitDensity = interp1(wavelength,data,targetWavelenth, 'pchip');

% Update the lens object with wavelength and density info
set(theLens,'wave', targetWavelenth);
set(theLens,'unitDensity',unitDensity);

%% Allocate space for storing an MTF slice for each shrew and calculated
% wavelength
mtfSlice = zeros(length(targetWavelength),length(TSindex),spatialsamples);

% Define figures to cycle through
psfMeshFig = figure;
otfMeshFig = figure;

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
    wvfP = wvfSet(wvfP, 'ref psf sample interval',psfSamplesPerMinute);

    % Loop over a list of defocus values and find one that maximizes the
    % Streh ratio. This is the same as maximizing the peak of the PSF.
    defocusDeltas = linspace(-1,1,40);
    defocusZ = zCoeffs_TreeShrew;
    for dd = 1:length(defocusDeltas)
        defocusZ1 = defocusZ;
        defocusZ1(5) = defocus0+defocusDeltas(dd);
        wvfP1 = wvfSet(wvfP,'zcoeffs',defocusZ1);
        wvfP1 = wvfComputePSF(wvfP1);
        psf1 = wvfGet(wvfP1,'psf',targetWavelength(1));
        peakPsf1(dd) = max(psf1(:));
    end
    figure; clf; hold on
    plot(defocusDeltas,peakPsf1,'ro','MarkerFaceColor','r','MarkerSize',12);
    title(sprintf('TS = %d, wavelength = %s',ts,num2str(targetWavelength(1))));

    [~,idx] = max(peakPsf1);
    maxStrehlDefocus(ts) = defocus0+defocusDeltas(idx);
    defocusZ1(5) = maxStrehlDefocus(ts);
    wvfP = wvfSet(wvfP,'zcoeffs',defocusZ1);
    wvfP = wvfComputePSF(wvfP);

    % Create optics structure with desired wvf object, and set other optics
    % properties.
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

    % Push lens density into oi and optics.
    % Update the oi.lens
    oi = oiSet(oi, 'lens', theLens);
    oi.optics.lens = theLens;

    % Loop over target wavelengths
    for ww = 1:length(targetWavelength)
        % Find wavelength in optics closest to target wavelength
        optics = oiGet(oi, 'optics');
        wavelengthSupport = opticsGet(optics, 'wave');
        [~,idx] = min(abs(wavelengthSupport-targetWavelength(ww)));
        targetWave = wavelengthSupport(idx);

        % Get PSF slice at target wavelength
        wavePSF = opticsGet(optics,'psf data',targetWave);

        % Extract PSF support in arcmin and plot a mesh of the PSF
        psfSupportMicrons = opticsGet(optics,'psf support','um');
        xGridMinutes = 60*psfSupportMicrons{1}/micronsPerDegree;
        yGridMinutes = 60*psfSupportMicrons{2}/micronsPerDegree;
        figure(psfMeshFig); clf;
        mesh(xGridMinutes,yGridMinutes,wavePSF);

        xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
        ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');
        [xSfGridCyclesDeg,ySfGridCyclesDeg,otfraw] = PsfToOtf(xGridMinutes,yGridMinutes,wavePSF);
        figure(otfMeshFig); clf;
        subplot(1,2,1);
        mesh(xSfCyclesDeg,ySfCyclesDeg,abs(otfraw));
        otf = psfCircularlyAverage(abs(otfraw));
        subplot(1,2,2);
        mesh(xSfCyclesDeg,ySfCyclesDeg,abs(otf));
        zlim([0 1]);

        waveMTF = abs(otf);
        [~,idx] = min(abs(ySfCyclesDeg));
        mtfSlice(ww,ts,:) = waveMTF(idx,:);
    end

end

%% Plot Avg MTF for each tree shrew and wavelength
MTFrows = 'B':'L';
for ts = 1:length(TSindex)
    for ww = 1:length(targetWavelength)

        figure; hold on;
        plot(xSfCyclesDeg, squeeze(mean(mtfSlice(ww,ts,:),2)), 'bo-', 'MarkerFaceColor', [0 0.8 1.0], 'MarkerSize', 10);
        set(gca, 'YTickLabel', 0:0.1:1, 'YTick', 0:0.1:1.0, 'YLim', [0 1.05]);
        ylabel('modulation');
        xlim([0 20]);
        xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'normal');
        title(sprintf('TS = %d, wavelength = %s',ts,num2str(targetWavelength(ww))));

        % Comparison data
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
        extraData = struct;
        extraData.sf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet','MTFs at max Strehl w Astig','Range','A4:A15'));
        rowstart = [MTFrows(TSindex(1)),'4'];
        rowend = [MTFrows(TSindex(end)),'15'];
        extraData.csf = cell2mat(readcell('Tree_Shrew_Aberrations_Remeasured_Oct2018.xlsx','Sheet','MTFs at max Strehl w Astig','Range',[rowstart,':',rowend]));
        extraData.legend = 'Saidak et al (2019)';
        extraData.ylabel = 'mtf';

        % Add to plot
        extraDataColors = [1 0 0];
        plot(extraData.sf, mean(extraData.csf,2), ...
            'rs-', 'MarkerSize', 10, ...
            'MarkerEdgeColor', squeeze(extraDataColors),  ...
            'MarkerFaceColor', squeeze(extraDataColors)*0.5+[0.5 0.5 0.5]);
        ylim([0 1.05]);
        ylabel(extraData.ylabel);
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