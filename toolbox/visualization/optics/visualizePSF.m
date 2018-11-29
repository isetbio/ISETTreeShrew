function visualizePSF(theOI, targetWavelength)
figure(); clf;
optics = oiGet(theOI, 'optics');
wavelengthSupport = opticsGet(optics, 'wave');
wavePSF = opticsGet(optics,'psf data',targetWavelength );

psfSupportMicrons = opticsGet(optics,'psf support','um');
xGridMinutes = 60*psfSupportMicrons{1}/optics.micronsPerDegree;
yGridMinutes = 60*psfSupportMicrons{2}/optics.micronsPerDegree;
xSupportMinutes = xGridMinutes(1,:);
ySupportMinutes = yGridMinutes(:,1);

psfRange = 70;
contourLevels = 0:0.05:1.0;
contourf(xSupportMinutes, ySupportMinutes, wavePSF/max(wavePSF(:)), contourLevels);
hold on;
[~,idx] = min(abs(ySupportMinutes));
psfSlice = wavePSF(idx,:)/max(wavePSF(:));
plot(xSupportMinutes, psfSlice*psfRange*2-psfRange, 'c-', 'LineWidth', 4.0);
plot(xSupportMinutes, psfSlice*psfRange*2-psfRange, 'b-', 'LineWidth', 1.0);
axis 'image'; axis 'xy';
grid on
set(gca, 'XLim', psfRange*1.05*[-1 1], 'YLim', psfRange*1.05*[-1 1], 'CLim', [0 1], ...
            'XTick', [-100:20:100], 'YTick', [-100:20:100]);
set(gca, 'FontSize', 16, 'XColor', [0 0 0], 'YColor', [0 0 0], 'Color', [0 0 0]);
xlabel('\it space (microns)');
set(gca, 'FontSize', 16);
cmap = brewermap(1024, 'greys');
colormap(cmap);
title(sprintf('PSF (%2.0f nm)', targetWavelength));
end
