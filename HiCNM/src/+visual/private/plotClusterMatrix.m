function plotClusterMatrix(D, ColorMap,caxislim)
% ------------------------------------- %
% --- plot cluster transition matrix -- %
% ----@created 2021-09-09 DN ---------- %
% ------------------------------------- %

nCluster   = utils.Parameters.instance.parameters.nClusters;
FigureBox  = utils.Parameters.instance.parameters.CTM_BoxSize;
TextSize   = utils.Parameters.instance.parameters.TextSize;
units      = utils.Parameters.instance.parameters.units;
LineWidth_Box   = utils.Parameters.instance.parameters.LineWidthBox;

if size(D,1) ~= nCluster
    disp('WARNING: Size of matrix ~= Nclusters ... using size of matrix!')
    nCluster = size(D,1);
end

% ColorMap      = vertcat([1 1 1],ColorMap);
ColorMap  = flipud(gray(50));
% ColorMap   = jet(50);
% ColorMap = getColormapForCTM();
% ColorMap = vertcat([1 1 1],ColorMap);
set(gca, 'position', [0.15 0.1 0.7 0.7])

%% PLOT figure
hold on
box on
colormap(ColorMap);

% PLOT matrix
imagesc(D);
Data_min = caxislim(1);
Data_max = caxislim(2);    
caxis([Data_min Data_max])
h = colorbar;
% h.Ticks=[0, 0.5, 1];
h.Position = [0.88 0.15 0.04 0.6];
set(get(h,'title'),'string','$D_{jk}$','interpreter','latex','FontSize',TextSize,'Rotation',0)
axis([0.5 nCluster+0.5 0.5 nCluster+0.5])
daspect([1 1 1])
set(gca,'YDir','reverse') % axis reverse      

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

set(gca, 'Fontsize', TextSize,'LineWidth',LineWidth_Box);
set(gcf, 'PaperUnits', units, 'PaperPosition', [0 0 1.05*FigureBox FigureBox]);

hold off

end
