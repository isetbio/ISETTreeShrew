function visualizeConeMosaicResponses(coneMosaic, coneExcitations)
    
    figure(); clf;
    subplotRows = 2;
    subplotCols = 3;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.01, ...
       'widthMargin',    0.02, ...
       'leftMargin',     0.02, ...
       'rightMargin',    0.12, ...
       'bottomMargin',   0.04, ...
       'topMargin',      0.01);

    instancesNum = size(coneExcitations,1);
    visualizedTrialsNum = min([subplotRows*subplotCols, instancesNum]);
    
    axHandle = subplot('Position', subplotPosVectors(1,1).v);
    coneMosaic.visualizeGrid(...
        'axesHandle', axHandle, ...
        'displayVisualDegs', true);
    set(gca, 'FontSize', 12, 'YTickLabel', {});
    xlabel('space (degs)')
    ylabel('');
    
    for k = 1:2
        axHandle = subplot('Position', subplotPosVectors(1,k+1).v);
        if (k == 2)
            coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', true, ...
                'labelColorBarTicks', true, ...
                'titleForColorBar', sprintf('R*/cone/%2.0fms', coneMosaic.integrationTime*1000));
        else 
            coneMosaic.renderActivationMap(axHandle, squeeze(coneExcitations(k,:,:)), ...
                'mapType', 'modulated disks', ...
                'signalRange', [0 max(coneExcitations(:))], ...
                'showColorBar', ~true, ...
                'labelColorBarTicks', ~true);
        end
        ylabel(''); set(gca, 'YTickLabel', {});
        
        xlabel('');
        set(gca, 'XTickLabel', {});
        
        set(gca, 'FontSize', 12);
        title(sprintf('trial #%d', k));
    end
    
    % Retrieve indices of cones along horizontal meridian
    [indicesOfConesAlongXaxis, xCoordsOfConesAlongXaxis, typesOfConesAlongXaxis] = indicesOfConesAlongHorizontalMeridian(coneMosaic);

    % Extract the excitations of cones along horizontal meridian
    coneExcitations = reshape(coneExcitations, [instancesNum  size(coneExcitations,2)* size(coneExcitations,3)]);
    coneExcitationsNxXY = zeros(instancesNum, numel(indicesOfConesAlongXaxis));
    for k = 1:numel(indicesOfConesAlongXaxis)
        coneExcitationsNxXY(:,k) = coneExcitations(:,indicesOfConesAlongXaxis(k));
    end
    
    % Plot the excitations separately for L-,M- and S-cones
    subplot('Position', [0.15 0.1 0.7 0.33]);
    idx = find(typesOfConesAlongXaxis == 2);
    LconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'r.');
    hold on;
    idx = find(typesOfConesAlongXaxis == 3);
    MconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'g.');
    idx = find(typesOfConesAlongXaxis == 4);
    SconesNum = numel(idx);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'b.');
    grid on
    set(gca, 'FontSize', 14, 'YTick', [0:25:200], 'XTick', [-0.3:0.1:0.3]);
    set(gca, 'YLim', [0 max(coneExcitations(:))]);
    xlabel('\it space (degs)');
    ylabel(sprintf('\\it R*/cone/%2.0fms', coneMosaic.integrationTime*1000));
    title(sprintf('%d trials, responses of %d L- %d M- and %d-S cones (along the horiz. meridian)', ...
        instancesNum, LconesNum, MconesNum, SconesNum), 'FontWeight', 'Normal', 'FontSize', 10);
end

function [indicesOfConesAlongXaxis,coneXcoordsAlongXaxis, theConeTypes] = indicesOfConesAlongHorizontalMeridian(theMosaic)
    
    conesXcoords = squeeze(theMosaic.patternSupport(1, :, 1));
    conesYcoords = squeeze(theMosaic.patternSupport(:, 1, 2));
    [coneXcoordsMesh, coneYcoordsMesh] = meshgrid(conesXcoords, conesYcoords);
    coneXcoordsDegsMesh = coneXcoordsMesh * 1e6 / theMosaic.micronsPerDegree;
    coneYcoordsDegsMesh = coneYcoordsMesh * 1e6 / theMosaic.micronsPerDegree;
    idx = find(theMosaic.pattern > 1);
    coneXcoordsDegsMesh = coneXcoordsDegsMesh(idx);
    coneYcoordsDegsMesh = coneYcoordsDegsMesh(idx);
    dx = diameterForCircularApertureFromWidthForSquareAperture(...
            theMosaic.pigment.width) * 1e6 / theMosaic.micronsPerDegree;
    indicesOfConesAlongXaxis = find(abs(coneYcoordsDegsMesh) < dx);
    coneXcoordsAlongXaxis = coneXcoordsDegsMesh(indicesOfConesAlongXaxis);
    theConeTypes = theMosaic.pattern(idx(indicesOfConesAlongXaxis));
    indicesOfConesAlongXaxis = idx(indicesOfConesAlongXaxis);
end
