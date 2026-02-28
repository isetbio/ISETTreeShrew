function t_treeShrewWavefrontDerivedOptics()
% Generate optics for all 11 trewshrew subjects
%
% Description:
%   Generates wavefront-derived optics for all 11 treeshrew subjects in SaidakEtAl_2019
%   and visualizes the corresponding PSFs and MTFs

% History:
%    02/28/26  NPC  Wrote it.

    mtfData = mtfTreeShrewFromPaper('SaidakEtAl_2019');
    mtfData = mtfData{1};
   

    % Generate trewshrew optics for all subjects
    treeshrewSubjectsNum = 11;

    % Visualize PSF and OTF at the in-focus wavelength (550 nm)
    theVisualizedWavelength = 550;

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1900 900], 'Color', [1 1 1]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 3, ...
       'colsNum', treeshrewSubjectsNum, ...
       'heightMargin',  0.03, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.05, ...
       'topMargin',      0.00);

    for iSubject = 1:treeshrewSubjectsNum
        theOI = oiTreeShrewCreate( ...
            'opticsType', 'wvf', ...
            'pupilDiameterMM', 4.0, ...
            'whichShrew', iSubject, ...
            'name', 'wvf-based optics');

        
        axPSF = subplot('Position', subplotPosVectors(1, iSubject).v);
        axOTF = subplot('Position', subplotPosVectors(2, iSubject).v);
        axOTFslice = subplot('Position', subplotPosVectors(3, iSubject).v);

        [mtfSupportCyclesPerDeg, theRadiallyAveragedMTFslice(iSubject,:)] = ...
        visualizePSFandOTFaAtSingleWavelength(...
            theOI, theVisualizedWavelength, axPSF, axOTF, axOTFslice, iSubject);
    end

    % Visualize the average (across all subjects) radially-averaged MTF
    hFig = figure(44); clf;
    set(hFig, 'Position', [10 10 1000 850], 'Color', [1 1 1]);
    axMTF = subplot(1,1,1);
    
    
    plot(mtfData.sf, mtfData.csf, 'bs--', ...
        'LineWidth', 1.0, ...
        'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'b', 'MarkerSize', 12);
    hold(axMTF, 'on');

    meanMTF = mean(theRadiallyAveragedMTFslice,1);
    targetSF = 20;
    [~, idx] = min(abs(mtfSupportCyclesPerDeg-targetSF));
    theMTFsfToMatch = mtfSupportCyclesPerDeg(idx);

    [~, iidx] = min(abs(mtfData.sf - theMTFsfToMatch));
    mtfData.sf (iidx)
    scalarToMatchCSF =  mtfData.csf(iidx) / meanMTF(idx);

    plot(axMTF, mtfSupportCyclesPerDeg, meanMTF*scalarToMatchCSF, ...
        'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);

    legend({'average radially-symmetric treeshrew MTF', 'treeshrew CSF (measured by SaidakEtAl)'});
    xlabel(axMTF, 'spatial frequency (cpd)');
    ylabel(axMTF, 'CSF');
    set(axMTF, 'XLim', [0 50], 'YLim', [0 100], 'FontSize', 16);
    set(axMTF, 'XTick', 0:5:50, 'YTick', [0:20:100]);
    grid(axMTF, 'on')

end

function [otfSupportCyclesPerDeg, theVisualizeMTFslice] = visualizePSFandOTFaAtSingleWavelength(...
    theOI, theVisualizedWavelength, axPSF, axOTF, axOTFslice, iSubject)


    % Get the optics data
    optics = oiGet(theOI, 'optics');

    % Get the focal length
    focalLengthMeters = opticsGet(optics, 'focal length');
    micronsPerDegree = focalLengthMeters*tand(1)*1e6;

    % Get the wavelength support
    wavelengthSupport = opticsGet(optics, 'otf wave');
    [~,theVisualizedWavelengthIndex] = min(abs(wavelengthSupport-theVisualizedWavelength));

    % Get the psf data
    psf = opticsGet(optics, 'psf data');
    theVisualizedPSF = squeeze(psf.psf(:,:,theVisualizedWavelengthIndex));

    % Get the psf support
    psfSupportXmicrons = squeeze(psf.xy(1,:,1));
    psfSupportYmicrons = squeeze(psf.xy(:,2,2));
   
     % Get the OTF
    theOTF = opticsGet(optics,'otf');
    theVisualizedOTF = squeeze(theOTF(:,:,theVisualizedWavelengthIndex));
    theVisualizedMTF = fftshift(abs(theVisualizedOTF));

    % Circularly symmetric MTF
    theVisualizedMTF = psfCircularlyAverage(theVisualizedMTF);

    % Get the OTF support
    otfSupportCyclesPerMM = opticsGet(optics,'otf fx', 'mm');
    otfSupportCyclesPerMicron = otfSupportCyclesPerMM * 1e-3;
    otfSupportCyclesPerDeg = otfSupportCyclesPerMicron * micronsPerDegree;

    % Whether to show the y-axis label
    showYaxisLabel = iSubject==1;

    % PSF visualization range
    visualizedPSFrangeMicrons = 20*[-1 1];

    % Visualize the PSF
    imagesc(axPSF, psfSupportXmicrons, psfSupportYmicrons, theVisualizedPSF);
    axis(axPSF, 'image'); axis(axPSF, 'xy')
    set(axPSF, 'XLim', visualizedPSFrangeMicrons, 'YLim', visualizedPSFrangeMicrons);
    set(axPSF, 'XTick', -20:10:20, 'yTick', -20:10:20, 'FontSize', 12);
    xlabel(axPSF, 'space, x (microns)');
    if (showYaxisLabel)
        ylabel(axPSF, 'space, y (microns)');
    else
        set(axPSF, 'yTickLabel', {});
    end
    grid(axPSF, 'on');
    title(axPSF, sprintf('PSF (subject %d)', iSubject));
    colormap(axPSF, 1-gray(1024));
    drawnow;


    sfSupportRangeCPD = [-15 15];
    imagesc(axOTF, otfSupportCyclesPerDeg, otfSupportCyclesPerDeg, theVisualizedMTF);
    axis(axOTF, 'image');  axis(axOTF, 'xy')
    set(axOTF, 'XLim', sfSupportRangeCPD, 'YLim', sfSupportRangeCPD);
    set(axOTF, 'XTick', -20:5:20, 'yTick', -20:5:20, 'FontSize', 12);
    xlabel(axOTF, 'spatial frequency, x (cpd)');
    if (showYaxisLabel)
        ylabel(axOTF, 'spatial frequency, y (cpd)');
    else
        set(axOTF, 'yTickLabel', {});
    end
    xtickangle(axOTF, 90)
    title(axOTF, 'MTF')
    grid(axOTF, 'on');
    colormap(axOTF, gray(1024));
    drawnow;

    [~,idx] = max(theVisualizedMTF(:));
    [peakRow, peakCol] = ind2sub(size(theVisualizedMTF), idx);
    theVisualizeMTFslice = squeeze(theVisualizedMTF(peakRow,:));
    
    plot(axOTFslice, otfSupportCyclesPerDeg, theVisualizeMTFslice, 'k-', 'LineWidth', 1.5);
    set(axOTFslice, 'XLim', [0 50], 'YLim', [0 1]);
    set(axOTFslice, 'XTick', 0:10:50, 'yTick', 0:0.2:1, 'FontSize', 12);
    grid(axOTFslice, 'on')
    xlabel(axOTFslice, 'spatial frequency, x (cpd)');
    if (showYaxisLabel)
        ylabel(axOTFslice, 'MTF');
    else
        set(axOTFslice, 'yTickLabel', {});
    end
    title(axOTFslice, 'radially-averaged MTF')
    drawnow;
end


function visualizePSFaAtWavelengths(figNo, theOI, displayedWavelengths, mtfData)

    % Get the optics data
    optics = oiGet(theOI, 'optics');
    
    % Get the OTF
    otf = opticsGet(optics,'otf');
    sfSupportXCyclesPerMM = opticsGet(optics,'otf fx', 'mm');
    sfSupportYCyclesPerMM = opticsGet(optics,'otf fy', 'mm');


    % Get the psf 
    psf = opticsGet(optics, 'psf data');
    xSupport = squeeze(psf.xy(1,:,1));
    ySupport = squeeze(psf.xy(:,2,2));
    
    thePSF = psf.psf;
    
    % Extract wavevelength support
    psfWavelengths = opticsGet(optics, 'otf wave');
    
    % Compute spatial support
    psfSampleSpacing = opticsGet(optics, 'psf spacing', 'um');
    psfSupport = (1:size(thePSF,1))*psfSampleSpacing;
    psfSupportMicrons = (psfSupport - mean(psfSupport));
    
    % Convert spatial support in minutes of arc
    focalLength = opticsGet(optics, 'focal length');
    micronsPerDegree = focalLength*tand(1)*1e6;
    psfSupportArcMin = psfSupportMicrons/micronsPerDegree*60;

    sfSupportXCyclesPerMicron = sfSupportXCyclesPerMM * 1e-3;
    sfSupportXCyclesPerDegree = sfSupportXCyclesPerMicron * micronsPerDegree;


    hFig = figure(figNo+1); clf;
    set(hFig, 'Position', [10 10 1000 900], 'Color', [1 1 1]);
    psfRangeArcMin = 8*[-1 1];

    %% Display the PSFs
    for k = 1:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-psfWavelengths));
        imagesc(psfSupportArcMin, psfSupportArcMin, squeeze(thePSF(:,:,idx))); hold on;
        plot([0 0], psfRangeArcMin, 'g-'); plot(psfRangeArcMin, [0 0],'g-'); axis 'xy';  axis 'image';
        set(gca, 'XLim', psfRangeArcMin, 'YLim', psfRangeArcMin, 'XTick', [-2:0.2:2], 'YTick', [-2:0.2:2], 'FontSize', 14);
        if (k == 1)
            xlabel('retinal space (arc min.)');
            ylabel('retinal space (arc min.)');
        else
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('%d nm', psfWavelengths(idx)));
    end
    colormap(gray);
    drawnow;

    sfRange = [0.1 30];
    xTicks = [0.1 0.3 1 3 10 30];
    hFig = figure(figNo+1); clf;
    set(hFig, 'Position', [10 10 1000 900], 'Color', [1 1 1]);

    for k = 1:numel(displayedWavelengths)
        subplot(3,4,k);
        [~,idx] = min(abs(displayedWavelengths(k)-psfWavelengths));

        theOTFatThisWave = fftshift(abs(squeeze(otf(:,:,idx))));

        theOTFatThisWave = psfCircularlyAverage(theOTFatThisWave);

        [~,iidx] = max(theOTFatThisWave(:));
        [row,col] = ind2sub(size(theOTFatThisWave), iidx);

        displayedSFindices = find(sfSupportXCyclesPerDegree>=0);
        theSpatialFrequencies = sfSupportXCyclesPerDegree(displayedSFindices);
        theSpatialFrequencies(1) = 0.1;
        theOTFslice = squeeze(theOTFatThisWave(row,displayedSFindices));

        plot(theSpatialFrequencies, theOTFslice, 'k-', 'LineWidth', 1.5);
        hold on;
        plot(mtfData.sf, mtfData.csf/max(mtfData.csf), 'ro');

        if (k == 1)
            xlabel('spatial frequency (c/deg)');
        else
            set(gca, 'XTick', [], 'YTick', []);
        end
        %set(gca, 'XScale', 'Log', 'XLim', sfRange, 'YScale', 'log', 'YLim', [3 500], 'XTick', xTicks, 'YTick', [1 3 10 30 100 300 1000], 'FontSize', 14);
        set(gca, 'XLim', [0 30], 'YLim', [0 1]);
        grid 'on'
        title(sprintf('%d nm', psfWavelengths(idx)));
    end
end
