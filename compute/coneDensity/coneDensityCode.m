%% Comparing human and treeshrew cone densities

% Driving question: Sure, tree shrews have worse vision than humans. But
% let's compare cones per visual angle, a unit that can be compared across
% different eye sizes, as a function of degrees from the retina. Humans
% will likely have better acuity in our retina, and given how poor we are
% in the far periphery, tree shrews will likely outperform humans there.
% But how far out in our visual periphery do we have to go before humans
% and tree shrews have similar visual acuity?

% An important note: For the comparison, we are considering the visual
% acuity to a spot in reference to the front of the head. For human vision,
% this would correspond roughly to our fovea. For tree shrews, because
% their eyes are at an angle due to their skull shape, we have incorporated
% an angle shift into the calculations.

% For the human cone density data, we use data from Curcio, 1990. For the
% tree shrew data, the cone density is approximated from Muller, 1989.

% drafted by jsc

%% Parameters

% Tree Shrew:
%   Which retina (from Muller, 1989) are you interested in analyzing? 
%   'A', B', or 'C'.
    retinaChoice = 'B';

%   Eye geometry:
    tsFocalLength = 5.81; %mm
    tsFOV = 175; %deg
    
% Human:
%   Eye geometry:
    humanFocalLength = 17; %mm
    humanFOV = 160; %deg



% Get the cone density as a function of visual angle for tree shrews and
% humans (right eye only)
[rightAngTS,rightConeDensityTS] = getRightTreeshrewCD('focalLength',tsFocalLength,'rightEyeFOV',tsFOV,'retinaChoice',retinaChoice);
treeshrewData = [rightAngTS ; rightConeDensityTS];

[rightAngHuman,rightConeDensityHuman] = getRightHumanCD('focalLength',humanFocalLength,'rightEyeFOV',humanFOV);
humanData = [rightAngHuman ; rightConeDensityHuman];

%Determine the visual angle where the cone density is the same across
%species. Uses InterX function from
T = InterX(treeshrewData,humanData);
I = T(:,size(T,2));

fprintf('At %.1f degrees eccentricity, both humans and tree shrews have approximately %.0f cones per visual angle.\n',I)

%% Plotting

yMin = min([rightConeDensityHuman rightConeDensityTS]);
yMax = max([rightConeDensityHuman rightConeDensityTS]);
xLims = [-100 150];

figure()
plot(rightAngHuman,rightConeDensityHuman, 'r-')
hold on
set(gca, ...
    'YScale','log',...
    'YLim',[yMin yMax],...
    'XLim',xLims)
plot(rightAngTS,rightConeDensityTS,'b-')
plot(I(1),I(2),'ko')

line([I(1),I(1)],[yMin,I(2)],'Color','black','LineStyle',':','LineWidth',1.5)
line([xLims(1),I(1)],[I(2),I(2)],'Color','black','LineStyle',':','LineWidth',1.5)

set(gca, 'XTick', unique([I(1), get(gca, 'XTick')]));
xt = get(gca, 'XTick');
xtl = cell([1,length(xt)]);
for i = 1:length(xt)
    xtl{i} = num2str(round(xt(i),2));
end
set(gca, 'XTickLabels', xtl);

set(gca, 'YTick', unique([yMin, I(2), yMax, get(gca, 'YTick')]));
yt = get(gca, 'YTick');
ytl = cell([1,length(yt)]);
for i = 1:length(yt)
    ytl{i} = num2str(round(yt(i)));
end
set(gca, 'YTickLabels', ytl);
legend('Human Cone Density','Treeshrew Cone Density')
title(['Visualizing Cone Density as a Function' newline 'of Visual Angle for Treeshrew ' retinaChoice])
xlabel(['Visual Angle - ',char(176)])
ylabel(['Cone Density - Log(Cones/',char(176),')'])