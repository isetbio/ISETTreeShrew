function visualizeDisplayLUT(gammaTable)
% Method to visualize the display's LUT
figure(); clf;
% generate linearRGBvalues
lutLength = size(gammaTable,1);
linearRGBvalue = (0:(lutLength-1))/lutLength;
% plot the RGB LUTs
plot(linearRGBvalue, gammaTable(:,1), 'r-', 'LineWidth', 1.5); hold on;
plot(linearRGBvalue, gammaTable(:,2), 'g-', 'LineWidth', 1.5);
plot(linearRGBvalue, gammaTable(:,3), 'b-', 'LineWidth', 1.5);
% label plot
xlabel('\it input RGB (settings)'); 
ylabel('\it normalized output intensity, rgb (primaries)'); 
axis 'square'
set(gca, 'FontSize', 16);
legend({'R primary', 'G primary', 'B primary'}, 'Location', 'NorthWest');
title('LUTs')
end