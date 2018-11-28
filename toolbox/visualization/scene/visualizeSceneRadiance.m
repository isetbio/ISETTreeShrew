function visualizeSceneRadiance(spatialSupport, spatialSupportUnits, ...
    photons, wavelengthSupport, wavelengthBandsToVisualize)
% Method to visualize the scene radiance at certain wavelength bands
spatialSupportX = squeeze(spatialSupport(1,:,1));
spatialSupportY = squeeze(spatialSupport(:,1,2));
photonRange = [min(photons(:)) max(photons(:))];

subplotRows = 3;
subplotCols = 3;
subplotPos = NicePlot.getSubPlotPosVectors(...
    'rowsNum', subplotRows, 'colsNum', subplotCols, ...
    'heightMargin',  0.07, 'widthMargin',    0.01, ...
    'leftMargin',     0.01, 'rightMargin',    0.01, ...
    'bottomMargin',   0.01, 'topMargin',      0.04);

figure(); clf;
cmap = brewermap(1024,'*Spectral');
colormap(cmap);
hProgressBar = waitbar(0,'Getting photon rates ...');
for iBand = 1:numel(wavelengthBandsToVisualize)
    if (iBand > subplotRows*subplotCols)
        continue
    end
    waitbar(iBand/numel(wavelengthBandsToVisualize), hProgressBar, ...
        sprintf('Getting photon rates at %2.0f nm', wavelengthBandsToVisualize(iBand)));
    row = floor((iBand-1)/subplotCols)+1; col = mod(iBand-1,subplotCols)+1;
    subplot('Position', subplotPos(row,col).v);
    % find index of spectral slices requested
    [~,visualizedWavelengthIndex] = ...
        min(abs(wavelengthSupport-wavelengthBandsToVisualize(iBand)));
    imagesc(spatialSupportX, spatialSupportY, ...
        squeeze(photons(:,:,visualizedWavelengthIndex)), photonRange);
    axis 'xy'; axis 'image';
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'FontSize', 16);
    title(sprintf('%2.0fnm', wavelengthSupport(visualizedWavelengthIndex)));
    colorbar();
end
close(hProgressBar);
drawnow;

end