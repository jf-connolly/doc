%---------------------------- Functions themselves ----------------------------%
%------------------------------------------------------------------------------%
%-- mpbGraph: 2D graphic of the multipeak benchmark and PSO swarm.
%-- Works with util_getName
%------------------------------------------------------------------------------%
%-- DEFINE
temp = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res = get(0,'ScreenPixelsPerInch')/2.56;
%-- END DEFINE

%------------------------------------------------------------------------------%
%-------------------------- User-defined parameters ---------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
nameAl    = 'mopso';
nameDb    = 'custom';
iteration = -1;
pbest     = 1;

%-- What to do
loaddata    = 0;
paretoFront = 1;
gArchive    = 0;

%-- Objective function global parameters
resolution = 60;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%------------------------------------------------------------------------------%
%------------------------------ Reading the data ------------------------------%
if loaddata
    %--------------------------------- Reading --------------------------------%
    fprintf('/*-- Load \n');

    %-- Swarm
    nameFile = sprintf('../savedStuff/%s_mopso.pso', nameDb);
    pso = util_readMopso(nameFile, 1, 1);

    nameFile = sprintf('../savedStuff/%s_mopso.arc', nameDb);
    archive  = util_readArchiveMopso(nameFile, 1, 1, pso.nIterations);
    
    nObjectives = pso.nObjectives;
    sizeSwarm   = pso.sizeSwarm;
    
    %--------------------- Formating the data for graphics --------------------%
    if   iteration < 0,   it = pso.nIterations;
    else                  it = iteration;
    end

    %-------------------- Particle position in both spaces --------------------%
    nFilled = archive.nFilled(it);
    x = [ archive.s(1:nFilled, 1, it),  archive(1:nFilled, 1, it)  ];
    y = [ archive.s(1:nFilled, 2, it),  archive(1:nFilled, 2, it)  ];
    z = [ 1-archive.sPm(1:nFilled, it), archive.sSz(1:nFilled, it) ];
else
    nObjectives = 2;
end

%------------------------------------------------------------------------------%
%----------------------------- Objective functions ----------------------------%
fprintf('/*-- Process\n');
xObj       = 0 : 1/resolution : 1;
[Y X]      = meshgrid(xObj,xObj);
nDimensions = 2;
counter     = 0;
x           = zeros((resolution+1)^2, 2);
f           = zeros(resolution+1, resolution+1, 2);
for i = 1: size(xObj,2);
    for j = 1: size(xObj,2);
        counter      = counter + 1;
        x(counter,:) = [xObj(i) xObj(j)];
    end
end

%------------------------- Kursawe - Original -------------------------%
if strcmp(nameDb,'custom')
    %-- Objective 1 ("y")
    std1 = 0.2;
    c1 = [0.6 0.6];
    f1 = 1 - exp( -(x(:,1)-c1(1)).^2 / std1  -(x(:,2)-c1(2)).^2 / std1 );

    %-- Objective 2 ("x")
%     std2 = 0.08;
%     c2 = [0.2 0.4];
%     c3 = [0.6 0.8];
    std2 = 0.1;
    c2 = [0.2 0.2];
    c3 = [0.8 0.8];
    f2 = 1.001 - exp(-((x(:,1)-c2(1)).^2/std2) - ((x(:,2)-c2(2)).^2/std2)) + ...
               - exp(-((x(:,1)-c3(1)).^2/std2) - ((x(:,2)-c3(2)).^2/std2));
else
    fprintf('/*-- No-name is not only a brand!\n');
    return;
end

f(:,:,1) = transpose(reshape(f1, resolution+1, resolution+1));
f(:,:,2) = transpose(reshape(f2, resolution+1, resolution+1));

%------------------------------------------------------------------------------%
%-------------------------------- Pareto Front --------------------------------%
if paretoFront
    nElements  = 0;
    nSolutions = size(f1,1);
    front = zeros(nSolutions, 2);
    frons = zeros(nSolutions, nDimensions);
    for n = 1:nSolutions;
        tested       = [f1(n) f2(n)];
        testex       = x(n,:);
        nonDominated = 1;

        for m = nElements:-1:1
            if tested(1) == front(m,1) && tested(2) == front(m,2)
                nonDominated = 1;
            elseif tested(1) >= front(m,1) && tested(2) >= front(m,2)
                nonDominated = 0;
            elseif tested(1) <= front(m,1) && tested(2) <= front(m,2)
                nElements = nElements - 1;
                front(m, :) = [];
                frons(m, :) = [];
            end
        end
        
        if nonDominated
            nElements           = nElements + 1;
            front(nElements, :) = tested;
            frons(nElements, :) = testex;
        end
    end
    
    front( nElements+1:size(front,1),: ) = [];
    frons( nElements+1:size(frons,1),: ) = [];
    
    front( 1,: ) = [];
    frons( 1,: ) = [];
    
    %------------------------ Other local pareto front ------------------------%
    %-- Filter
    nSolutions = size(f1,1);
    xo = x;    f1o = f1;    f2o = f2;
    for n = nSolutions:-1:1
        if x(n,1) >= c1(1) || x(n,2) >= c1(2)
            xo(n,:) = [];    f1o(n) = [];    f2o(n) = [];
        end
    end
    
    %-- Local Pareto front
    nElements  = 0;
    nSolutions = size(f1o,1);
    fronto = zeros(nSolutions, 2);
    fronso = zeros(nSolutions, nDimensions);
    for n = 1:nSolutions;
        tested       = [f1o(n) f2o(n)];
        testex       = xo(n,:);
        nonDominated = 1;

        for m = nElements:-1:1
            if tested(1) == fronto(m,1) && tested(2) == fronto(m,2)
                nonDominated = 1;
            elseif tested(1) >= fronto(m,1) && tested(2) >= fronto(m,2)
                nonDominated = 0;
            elseif tested(1) <= fronto(m,1) && tested(2) <= fronto(m,2)
                nElements = nElements - 1;
                fronto(m, :) = [];
                fronso(m, :) = [];
            end
        end
        
        if nonDominated
            nElements           = nElements + 1;
            fronto(nElements, :) = tested;
            fronso(nElements, :) = testex;
        end
    end
    
    fronto( nElements+1:size(fronto,1),: ) = [];
    fronso( nElements+1:size(fronso,1),: ) = [];
    
    for n = 1:size(fronto,1)
        if fronso(n,1) == c1(1) && fronso(n,2) == c1(2)
            fronso(n,:) = [];
            fronto(n,:) = [];
            break;
        end
    end
end

%------------------------------------------------------------------------------%
%----------------------------- Graphic parameters -----------------------------%
%-- General
graphWidth = 16.64 / 3; %-- centimeters
ratio      = 1/1;
nLevels    = 50;
ax         = [0 65; 0.3 1; -0.6 0.5; 0 0.9];

%-- Marker & font
mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
sizeMkr = 5;   sizeFont = 10;   widthLine = 0.3;
clrMin  = 0.3;   clrMax   = 1;    pow       = 1;

%-- Figure initialization
close all;
for o = 1:nObjectives+1
    width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    %     posScrn  = [1300 sizeScrn(2)-104-550*(o-1) width height];
    posScrn  = [1200 sizeScrn(2)-104-400*(o-1) width height];
    posPaper = [0 0 graphWidth graphWidth/ratio];

    fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
    if o <= nObjectives, ax = axes('Position', [0 0 1 1]); end
end

%------------------------------------------------------------------------------%
%-------------------------------- Search space --------------------------------%
fprintf('/*-- Figure\n');
for o = 1:nObjectives
    
    set(0,'CurrentFigure',o);
    hold on
    
    %-- Figure - objective function
    [tmp ct] = contourf(X, Y, f(:,:,o), nLevels);
    set(ct, 'LineStyle', '-', 'LineWidth', 0.2, 'LineColor', [0.1 0.1 0.25]);
    colormap bone;

    %--------------------------------- Archive --------------------------------%
    if gArchive
        clr = [0.5417    0.6250    0.6667];
        for m = 1:archive.nFilled(it)
            pl = plot(x, y, 'o', 'Color', [0 0 0], ...
              'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
              'LineWidth', widthLine);
        end
    end

    %------------------------ Other local Pareto front ------------------------%
    clr = [0.9566    0.9722    0.9722];
    opf = plot(fronso(:,1), fronso(:,2), 'o', 'Color', [0 0 0], ...
         'markerfacecolor', clr, 'MarkerSize',sizeMkr,...
         'LineWidth', widthLine);
     
    %------------------------------ Pareto front ------------------------------%
    clr = [0.32    0.32    0.44];
    pf = plot(frons(:,1), frons(:,2), 'o', 'Color', [0 0 0], ...
         'markerfacecolor', clr, 'MarkerSize',sizeMkr,...
         'LineWidth', widthLine);
        
    hold off;
    
    box  on;
%     set(gca, 'XTick', [], 'YTick', [], 'XLim', [0 1], 'YLim', [0 1]);

    set(gca, 'Position', [0.09 0.09 0.9 0.9], 'XLim', [0 1], 'YLim', [0 1], ...
        'FontName', 'Times New Roman');
    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

    xlabel('$h_1$', 'FontName', 'Times New Roman',...
     'FontSize', sizeFont, 'Interpreter','latex');
    ylabel('$h_2$', 'FontName', 'Times New Roman',...
     'FontSize', sizeFont, 'Interpreter','latex');
    

    %-- Saving the figure
    nameFile = sprintf('../figures/%s_mopso_obj%d', nameDb, o);
    nameFig  = sprintf('-f%d', o);
    saveas(fig(o), nameFile, 'pdf');
%     print(nameFig, '-dpng', '-r300', nameFile);
end

%------------------------------------------------------------------------------%
%------------------------------- Objective space ------------------------------%
set(0,'CurrentFigure',o+1);
hold on

%------------------------------- Existing space -------------------------------%
clr = [0.7179    0.8194    0.8194];
plot(f1, f2, 'o', 'Color', clr, 'markerfacecolor', clr, ...
    'MarkerSize',sizeMkr-3);

%-------------------------------- Pareto front --------------------------------%
clr = [0.32    0.32    0.44];
plot(front(:,1), front(:,2), 'o', 'Color', [0 0 0], ...
  'markerfacecolor', clr, 'MarkerSize',sizeMkr-2,...
  'LineWidth', widthLine);

%-------------------------- Other local pareto front --------------------------%
clr = [0.9566    0.9722    0.9722];
plot(fronto(:,1), fronto(:,2), 'o', 'Color', [0 0 0], ...
     'markerfacecolor', clr, 'MarkerSize',sizeMkr-2,...
     'LineWidth', widthLine);


% %------------------------- Ensemble for solutions -------------------------%
% clr = [0.8481    0.9028    0.9028];
% for i = 1:resolution+1
%     for j = 1:resolution+1
%         if f(i,j,2) < 0.2
%             plot(f(i,j,1), f(i,j,2), 'o', 'Color', clr, ...
%                  'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
%                  'LineWidth', widthLine);
%              pb=1;
%         end
%     end
% end

%----------------------------------- Archive ----------------------------------%
% clr = [0.8481    0.9028    0.9028];
% for n = 1:archive.nFilled(it)
%     pl = plot(zArc(n,1), zArc(n,2), 'o', 'Color', [0 0 0], ...
%       'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
%       'LineWidth', widthLine);
% end

hold off
 
set(gca, 'Position', [0.11 0.10 0.88 0.89], 'XLim', [0 1], 'YLim', [0 1], ...
    'FontName', 'Times New Roman');
set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

xlabel('$f_{1}(\textbf{h})$', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','latex');

ylabel('$f_{2}(\textbf{h})$', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','latex');

%------------------------------ Saving the figure -----------------------------%
nameFile = sprintf('../figures/%s_mopso_pf', nameDb );
nameFig  = sprintf('-f%d', o+1);

saveas(fig(3), nameFile, 'pdf');
print(nameFig, '-dpng', '-r400', nameFile);


%------------------------------------------------------------------------------%
%----------------------------------- Legend -----------------------------------%
o = 4;
graphWidth = 16.64 * 0.6; %-- centimeters
ratio = 13;
width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
%     posScrn  = [1300 sizeScrn(2)-104-550*(o-1) width height];
posScrn  = [1000 sizeScrn(2)-104-200*(o-1) width height];
posPaper = [0 0 graphWidth graphWidth/ratio];

fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
            'PaperSize', [posPaper(3) posPaper(4)],...
            'PaperPosition', posPaper, 'Color', [1,1,1]);
axes('Position', [0 0 1 1]);

hold on
%-- Marker & text - Pareto front
clr = [0.3056    0.3056    0.4253];
plot(0, 0, 'o', 'Color', [0 0 0], 'markerfacecolor', clr, ...
     'MarkerSize',sizeMkr+2, 'LineWidth', widthLine);

text(0.25, 0, 'Pareto front', 'FontName','Times New Roman', ...
      'FontSize', sizeFont, 'HorizontalAlignment', 'left', ...
      'Interpreter', 'Latex');

%-- Marker & text - Pareto front
clr = [0.9566    0.9722    0.9722];
plot(3, 0, 'o', 'Color', [0 0 0], 'markerfacecolor', clr, ...
     'MarkerSize',sizeMkr+2, 'LineWidth', widthLine);

text(3.25, 0, 'Other local Pareto front', 'FontName','Times New Roman', ...
      'FontSize', sizeFont, 'Interpreter', 'Latex');
hold off
grid off
box on

set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

set(gca, 'Position', [0.01 0.08 0.98 0.9], ...
          'XLim', [-0.35, 7], 'YLim', [-1, 1]);

nameFile = sprintf('../figures/%s_mopso_leg', nameDb );
nameFig  = sprintf('-f%d', 4);

saveas(fig(4), nameFile, 'pdf');
print(nameFig, '-dpng', '-r300', nameFile);
      
fprintf('/*---------------------- End ----------------------*/\n');