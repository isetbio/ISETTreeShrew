function visualizeDisplayGamut(primariesXYZ)
% Method to visualize the display's gamut
% Extract the maximum luminance for each primary (Y tristimulus value)
maxLuminanceCdPerM2 = primariesXYZ(:,2);
% Extract the (x,y) chromaticity coordinates, e.g. x = X / (X+Y+Z)
xChroma = primariesXYZ(:,1) ./ sum(primariesXYZ,2);
yChroma = primariesXYZ(:,2) ./ sum(primariesXYZ,2);

figure(); clf;
plot(xChroma(1), yChroma(1), 'rs', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.5 0.5]); hold on
plot(xChroma(2), yChroma(2), 'gs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 1 0.5]);
plot(xChroma(3), yChroma(3), 'bs', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.5 1]);
xx = [xChroma; xChroma(1)]; yy = [yChroma; yChroma(1)];
plot(xx,yy, 'k--', 'LineWidth', 1.5); hold on;
for k = 1:3
text(xChroma(k)+0.03, yChroma(k), sprintf('%2.1f cd/m2', maxLuminanceCdPerM2(k)), 'FontSize', 14);
end
set(gca, 'XLim', [0 0.8], 'YLim', [0 0.8], 'XTick', 0:0.1:1, 'YTick', 0:0.1:1, 'FontSize', 16);
axis 'square'; grid on; 
xlabel('\it x-chroma'); ylabel('\it y-chroma');
legend({'R primary', 'G primary', 'B primary'})
end