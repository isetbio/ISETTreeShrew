%% Code for plotting human and treeshrew cone densities
% -Simplified to a single eye-

% Driving question: Sure, tree shrews have worse vision than humans. But
% let's compare cones per visual angle, a unit that can be compared
% across different eye sizes, as a function of degrees from the retina.
% Humans will likely have better acuity in our retina, but how far out in
% our visual periphery do we have to go before humans and tree shrews have
% similar visual acuity?


% drafted by jsc

% Which retina (from Muller, 1989) are you interested in analyzing? 'A',
% 'B', or 'C'.

retinaChoose = 'B';

% Get the cone density as a function of visual angle for tree shrews and
% humans (right eye only)
[rightAngTS,rightConeDensityTS] = getRightTreeshrewCD(retinaChoose);
treeshrewData = [rightAngTS ; rightConeDensityTS];
[rightAngHuman,rightConeDensityHuman] = getRightHumanCD;
humanData = [rightAngHuman ; rightConeDensityHuman];

%Determine the visual angle where the cone density is the same across
%species. Uses InterX function from 
T = InterX(treeshrewData,humanData);
I = T(:,size(T,2));

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
title(['Visualizing Cone Density as a Function' newline 'of Visual Angle for Treeshrew ' retinaChoose])
xlabel(['Visual Angle (',char(176),')'])
ylabel('Combined Cone Density (cones/degree)')