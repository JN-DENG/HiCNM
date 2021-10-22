function plotProbMatrix(Data, ColorMap)

% ------------------------------------- %
% --- plot cluster transition matrix -- %
% ----@created 2021-09-09 DN ---------- %
% ------------------------------------- %

nCluster   = size(Data,1);
FigureBox  = utils.Parameters.instance.parameters.CTM_BoxSize*(nCluster/10);
TextSize   = utils.Parameters.instance.parameters.TextSize+2*(nCluster/10);
units      = utils.Parameters.instance.parameters.units;
LineWidth_Box   = utils.Parameters.instance.parameters.LineWidthBox*(nCluster/10);

set(gca, 'position', [0.07 0.05 0.8 0.8])
% ----------------------------------------------------------------------------------------------- %
% --------    Plot of Matrix -------------------------------------------------------------------- %
% ----------------------------------------------------------------------------------------------- %

%% START figure
box on
hold on
% colormap(flipud(hot));
colormap(ColorMap);

% PLOT matrix
imagesc(Data);
Data_min = 0;
Data_max = 1;
caxis([Data_min Data_max])
h = colorbar;
h.Ticks=[0, 0.5, 1];
h.Position = [0.88 0.15 0.04 0.6];
set(get(h,'title'),'string','$P_{jk}$','interpreter','latex','FontSize',TextSize,'Rotation',0)
axis([0.5 nCluster+0.5 0.5 nCluster+0.5])
daspect([1 1 1])
set(gca,'YDir','reverse')    

% LINES
for irow = 1:nCluster+1
    plot([irow-0.5,irow-0.5], [0.5,nCluster+0.5], '-k', 'LineWidth',0.3) % horizontal up
    plot([0.5,nCluster+0.5], [irow-0.5,irow-0.5], '-k', 'LineWidth',0.3) % horizontal below
end
% Border
plot([0.5,0.5], [0.5,nCluster+0.5], '-k', 'LineWidth',1) % vertical left
plot([nCluster+0.5,nCluster+0.5], [0.5,nCluster+0.5], '-k', 'LineWidth',1) % vertical right
plot([0.5,nCluster+0.5], [0.5,0.5], '-k', 'LineWidth',1) % horizontal up
plot([0.5,nCluster+0.5], [nCluster+0.5,nCluster+0.5], '-k', 'LineWidth',1) % horizontal down

% Label
if nCluster<15
    set(gca,'xtick',1:nCluster,'ytick',1:nCluster,'XAxisLocation','top')
else  
    set(gca,'xtick',[1,round(([1:5]/5)*nCluster)],'ytick',[1,round(([1:5]/5)*nCluster)],'XAxisLocation','top')
end
xlabel('$\mathcal{C}_k$','interpreter','latex','FontSize',TextSize)
ylabel('$\mathcal{C}_j$','interpreter','latex','FontSize',TextSize,'Rotation',0)

set(gca, 'Fontsize', TextSize,'LineWidth', LineWidth_Box);
set(gcf, 'PaperUnits', units, 'PaperPosition', [0 0 1.1*FigureBox FigureBox]);

hold off

%% Finished
disp(['Finished: Plot of CTM.P'])
