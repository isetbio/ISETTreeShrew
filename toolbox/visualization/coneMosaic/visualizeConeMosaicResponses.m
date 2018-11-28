function visualizeConeMosaicResponses(coneMosaic, coneExcitations)
    
    subplotRows = 2;
    subplotCols = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', subplotRows, ...
       'colsNum', subplotCols, ...
       'heightMargin',  0.06, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.07, ...
       'rightMargin',    0.03, ...
       'bottomMargin',   0.15, ...
       'topMargin',      0.1);

    instancesNum = size(coneExcitations,1);
    visualizedTrialsNum = min([subplotRows*subplotCols, instancesNum]);
    
    for k = 1:visualizedTrialsNum-1
        r = floor((k-1)/subplotCols)+1;
        c = mod(k-1,subplotCols)+1;
        axHandle = subplot('Position', subplotPosVectors(r,c).v);
        if (k == 1)
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
        if (c>1)
        	set(gca, 'YTickLabel', {});
            ylabel('');
        end
        xlabel(''); set(gca, 'XTickLabel', {});
        set(gca, 'FontSize', 14);
        title(sprintf('trial #%d', k));
    end
    
    % Rertieve indices of cones along horizontal meridian
    [indicesOfConesAlongXaxis, xCoordsOfConesAlongXaxis, typesOfConesAlongXaxis] = indicesOfConesAlongHorizontalMeridian(coneMosaic);

    % Extract the excitations of cones along horizontal meridian
    coneExcitations = reshape(coneExcitations, [instancesNum  size(coneExcitations,2)* size(coneExcitations,3)]);
    coneExcitationsNxXY = zeros(instancesNum, numel(indicesOfConesAlongXaxis));
    for k = 1:numel(indicesOfConesAlongXaxis)
        coneExcitationsNxXY(:,k) = coneExcitations(:,indicesOfConesAlongXaxis(k));
    end
    
    % Plot the excitations separately for L-,M- and S-cones
    subplot('Position', subplotPosVectors(subplotRows,subplotCols).v);
    idx = find(typesOfConesAlongXaxis == 2);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'r.');
    hold on;
    idx = find(typesOfConesAlongXaxis == 3);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'g.');
    idx = find(typesOfConesAlongXaxis == 4);
    plot(xCoordsOfConesAlongXaxis(idx), coneExcitationsNxXY(:,idx), 'b.');
    axis 'square'
    grid on
    set(gca, 'FontSize', 14, 'YTick', [0:25:200], 'XTick', [-0.3:0.1:0.3]);
    xlabel('\it space (degs)');
    ylabel(sprintf('\\it R*/cone/%2.0fms', coneMosaic.integrationTime*1000));
    title(sprintf('%d trials', instancesNum));
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
