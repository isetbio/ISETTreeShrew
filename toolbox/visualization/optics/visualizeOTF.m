function visualizeOTF(theOI, targetWavelength)
figure(); clf;

optics = oiGet(theOI, 'optics');
waveOTF = fftshift(opticsGet(optics,'otf data',targetWavelength));
xSfCyclesDeg = opticsGet(optics, 'otf fx', 'cyclesperdeg');
ySfCyclesDeg = opticsGet(optics, 'otf fy', 'cyclesperdeg');

otfRangeCyclesPerDeg = 3;
theTicks = [-10:0.5:10];
xx = find(abs(xSfCyclesDeg) <= otfRangeCyclesPerDeg);
yy = find(abs(ySfCyclesDeg) <= otfRangeCyclesPerDeg);
    
waveOTF = abs(waveOTF(yy,xx));
xSfCyclesDeg = xSfCyclesDeg(xx);
ySfCyclesDeg = ySfCyclesDeg(yy);
  
waveOTF = waveOTF / max(abs(waveOTF(:)));
contourLevels = 0:0.05:1.0;
contourf(xSfCyclesDeg, ySfCyclesDeg, waveOTF, contourLevels);
hold('on');
axis('image'); axis('xy');
set(gca, 'ZLim', [0 1], 'CLim', [0 1], 'XLim', otfRangeCyclesPerDeg *1.05*[-1 1], 'YLim', otfRangeCyclesPerDeg *1.05*[-1 1], 'FontSize', 14);
set(gca, 'XTick', theTicks, 'YTick', theTicks);
set(gca, 'YTickLabel', {});
grid('on'); box('on');
cmap = brewermap(1024, 'greys');
set(gca, 'FontSize', 16);
colormap(cmap);
xlabel('\it spatial frequency (c/deg)', 'FontWeight', 'bold');   
title(sprintf('MTF (%2.0f nm)', targetWavelength));
end

