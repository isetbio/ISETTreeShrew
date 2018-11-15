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
xlabel('\it settings value (linear value, RGB)'); 
ylabel('\it primary value (light intensity, rgb)'); 
axis 'square'
set(gca, 'FontSize', 16);
legend({'R primary', 'G primary', 'B primary'}, 'Location', 'NorthWest')
end