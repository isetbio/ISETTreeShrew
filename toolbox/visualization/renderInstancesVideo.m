function renderInstancesVideo(figNo, imIndex, oiArray, coneExcitationArray, theConeMosaic, timeAxis)
% Render a video of the response (first trial only) over time.
%
% Syntax:
%   renderInstancesVideo(figNo, imIndex, oiArray, coneExcitationArray, theConeMosaic)
%
% Description:
%    Render a video of the response (first trial only) over time.
%
%
% Inputs:
%    figNo                      Figure number, 1, 2 etc.
%    imIndex
%    oiArray                    Array with the ois to be visualized
%    coneExcitationArray        Array with the cone mosaic responses to be visualized
%    demosaicedExcitationArray  Array with the demosaiced cone responses to be visualized
%    theConeMosaic              The employed cone mosaic
%    timeAxis                   The response time axis
%
% Outputs:
%    None
%

% History:
%    1/24/2019  NPC   Wrote it


    videoFileName = sprintf('LeMem%d_ResponseInstances', imIndex);
    videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 30;
    videoOBJ.Quality = 100;
    videoOBJ.open();
        
    hFig = figure(figNo); 
    clf;
    set(hFig, 'Position', [10 10 1450 700], 'Color', [1 1 1]);
    
    subplotRows = 1;
    subplotCols = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.04, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.08, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);
    
    oiRGB = oiGet(oiArray{imIndex}, 'rgb');
    spatialSupportMM = oiGet(oiArray{imIndex}, 'spatial support', 'mm');
    optics = oiGet(oiArray{imIndex}, 'optics');
    focalLength = opticsGet(optics, 'focal length');
    mmPerDegree = focalLength*tand(1)*1e3;
    spatialSupportDegs = spatialSupportMM/mmPerDegree;
    spatialSupportXDegs = spatialSupportDegs(1,:,1);
    spatialSupportYDegs = spatialSupportDegs(:,1,2);
    spatialRangeDegs = max(theConeMosaic.fov)*0.5*[-1 1];
    ticksDegs = [-5:1:5];
    
    subplot('Position', subplotPosVectors(1,1).v);
    image(spatialSupportXDegs,spatialSupportYDegs, oiRGB);
    axis 'image'; axis 'ij'; grid on
    set(gca, 'FontSize', 18, ...
        'XLim', spatialRangeDegs, 'YLim', spatialRangeDegs, ...
        'XTick', ticksDegs, 'XTickLabels', sprintf('%2.1f\n', ticksDegs), ...
        'YTick', ticksDegs, 'YTickLabels', {} ...
        );
    xlabel('\it space (degs on retina)');
    title('retinal image');

    axHandle = subplot('Position', subplotPosVectors(1,2).v);
    theConeExcitationsAllInstances = coneExcitationArray{imIndex};
    
    %responseRange = prctile(theConeExcitationsAllInstances(:), [5 95]);
    responseRange = [min(theConeExcitationsAllInstances(:)) max(theConeExcitationsAllInstances(:))];
    
    iTrial = 1;
    tmp = theConeMosaic.reshapeHex2DmapToHex3Dmap(squeeze(theConeExcitationsAllInstances(iTrial,:,:)));
    
    for tBin = 1:numel(timeAxis)
        theConeMosaic.renderActivationMap(axHandle, squeeze(tmp(:,:,tBin)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 responseRange(2)], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', theConeMosaic.integrationTime*1000));

        ticksMeters = ticksDegs * mmPerDegree * 1e-3;
        axis 'image'; axis 'ij'; 
        set(gca, 'FontSize', 18, 'YTickLabel', {}, 'XTick', ticksMeters, 'XTickLabels', ticksDegs);

        xlabel('\it space (degs on retina)');
        ylabel('');
        title(sprintf('mosaic activation (%2.0f ms)', 1000*timeAxis(tBin)));
        drawnow;
        
        % Add video frame
        videoOBJ.writeVideo(getframe(hFig));
    end % tBin
    
    % Close video stream
    videoOBJ.close();
    fprintf('File saved in %s\n', videoFileName);
end