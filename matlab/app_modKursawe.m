%------------------------------------------------------------------------------%
%-- mpbGraph: 2D graphic of the multipeak benchmark and PSO swarm.
%-- Works with util_getName
%------------------------------------------------------------------------------%
%-- DEFINE
temp     = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res      = get(0,'ScreenPixelsPerInch')/2.56;
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
ax          = [0 65; 0.3 1; -0.6 0.5; 0 0.9];
nLevels     = 40;
resolution  = 100;
nObjectives = 2;

modified    = 0;
compute     = 1;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%----------------------------- Objective functions ----------------------------%
if compute
    xObj  = -5 :10/resolution :5;
    [Y X] = meshgrid(xObj,xObj);

    f    = zeros(resolution, resolution, nObjectives);

    %-- Objective function
    x = zeros(resolution+1, resolution+1, resolution+1, 2);
    for i = 1: size(xObj,2);
        for j = 1: size(xObj,2);

            x1 = xObj(i);
            x2 = xObj(j);

            x(i,j,1) = x1;
            x(i,j,2) = x2;

            if modified
                %-------------------------- Modified --------------------------%
                f(i,j,1) = - ( 10*exp(-0.2*sqrt( x1^2 + x2^2 )) );
                f(i,j,2) = - ( abs(x1)^0.8 + 5*(sin(x1)^3) + ...
                               abs(x2)^0.8 + 5*(sin(x2)^3) );
            else
                %-------------------------- Kursawe ---------------------------%
                f(i,j,1) = - ( 10*exp(-0.2*sqrt( x1^2 + x2^2 )));
                f(i,j,2) =     abs(x1)^0.8 + 5*(sin(x1)^3) +...
                               abs(x2)^0.8 + 5*(sin(x2)^3);
            end

        end
    end

    x1 = x(:,:,1);  x1 = x1(:);
    x2 = x(:,:,2);  x2 = x2(:);
    
    f1 = f(:,:,1);  f1 = f1(:);
    f2 = f(:,:,2);  f2 = f2(:);
    
%     [min(f1) max(f1); min(f2) max(f2)]

    %------------------------------ Pareto Front ------------------------------%
    front = zeros(resolution^2,2);
    frons = zeros(resolution^2,2);
    nElements = 0;
    nSolutions = resolution^2;
    for n = 1:nSolutions;
        tested       = [f1(n) f2(n)];
        testes       = [x1(n) x2(n)];
        nonDominated = 1;

        for m = nElements:-1:1
            if tested(1) == front(m,1) && tested(2) == front(m,2)
                nonDominated = 1;
            elseif tested(1) > front(m,1) && tested(2) > front(m,2)
                nonDominated = 0;
            elseif tested(1) < front(m,1) && tested(2) < front(m,2)
                nElements = nElements - 1;
                front(m, :) = [];
            end
        end

        if nonDominated
            nElements           = nElements + 1;
            front(nElements, :) = tested;
            frons(nElements, :) = testes;
        end
        bp=1;
    end
    
    front( nElements+1:size(front,1),: ) = [];
    frons( nElements+1:size(frons,1),: ) = [];
end

%--------------------------------- Save'n'load --------------------------------%
% if savePF
%     nameFile = sprintf('../savedStuff/%s_pf.mat', nameDb);
%     save(nameFile, 'front', 'frons', 'x', 'f', 'f1', 'f2', '-mat');
% end
% 
% if loadPF
%     nameFile = sprintf('../savedStuff/%s_pf.mat', nameDb);
%     load(nameFile, 'front', 'frons', 'x', 'f', 'f1', 'f2', '-mat');
% end

%------------------------------------------------------------------------------%
%----------------------------- Graphic parameters -----------------------------%
%-- General
graphWidth = 8;
ratio      = 1/1;

%-- Marker & font
mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
sizeMkr = 7;   sizeFont = 10;   widthLine = 0.5;
clrMin  = 0;   clrMax   = 1;    pow        = 1;

%-- Figure initialization
close all;

% for o=1:nObjectives+1
%     width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
%     posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
%     posPaper = [0 0 graphWidth graphWidth/ratio];
% 
%     fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
%     set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
%                 'PaperSize', [posPaper(3) posPaper(4)],...
%                 'PaperPosition', posPaper, 'Color', [1,1,1]);
%     ax = axes('Position', [0 0 1 1]);
%     
% end


%-------------------------------- Search space --------------------------------%
fprintf('/*-- Search space\n');
for o = 1:nObjectives
    
%     set(0,'CurrentFigure', o);
    width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
    posPaper = [0 0 graphWidth graphWidth/ratio];

    fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
%     ax = axes('Position', [0 0 1 1]);

    hold on
        %-- Figure - objective function
        [tmp ct] = contourf(X, Y, f(:,:,o), nLevels);
        set(ct, 'LineStyle', '-', 'LineWidth', 0.2, ...
                'LineColor', [0.1 0.1 0.25]);
        colormap bone;

        %-- Pareto front
        clr = [0.3056    0.3056    0.4253];
        plot(frons(:,1), frons(:,2), 'o', ...
             'markeredgecolor', clr-0.3, 'markerfacecolor', clr,        ...
             'MarkerSize',2, 'LineWidth', widthLine);

        plot(frons(:,2), frons(:,1), 'o', ...
             'markeredgecolor', clr-0.3, 'markerfacecolor', clr,        ...
             'MarkerSize',2, 'LineWidth', widthLine);
         
    hold off
    
    box off;
%     set(gca, 'XTick', [], 'YTick', []);
    
    if modified,   nameFile = sprintf('../figures/kurMod_%d', o);
    else           nameFile = sprintf('../figures/kurOri_%d', o);
    end
%     saveas(fig(o), nameFile, 'pdf');
end

%---------------------------- Pareto Optimal front ----------------------------%
fprintf('/*-- Objective space\n');
% set(0,'CurrentFigure', nObjectives+1);

o = 3;
width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
posPaper = [0 0 graphWidth graphWidth/ratio];
fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
            'PaperSize', [posPaper(3) posPaper(4)],...
            'PaperPosition', posPaper, 'Color', [1,1,1]);
ax = axes('Position', [0 0 1 1]);

hold on 

    %-- Objective space
    clr = [0.7179    0.8194    0.8194];
    plot(f1, f2, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
        'MarkerSize',1);
    
    %-- Pareto front
    clr = [0.3056    0.3056    0.4253];
    plot(front(:,1), front(:,2), 'o', 'Color', [0 0 0], ...
        'markerfacecolor', clr, 'MarkerSize', 2, 'LineWidth', widthLine);
hold off

set(gca, 'Position', [0.15 0.13 0.82 0.85], ...
'FontName', 'Times New Roman');%, 'XLim', [0 1], 'YLim', [0 1]);

xlabel('Objective 1', 'FontSize', sizeFont, 'Interpreter','latex');
ylabel('Objective 2', 'FontSize', sizeFont, 'Interpreter','latex');

if modified,   nameFile = sprintf('../figures/kurMod_pf');
else           nameFile = sprintf('../figures/kurOri_pf');
end
saveas(fig(o), nameFile, 'pdf');
% print('-f3', '-dtiff', nameFile, '-r300');

fprintf('/*---------------------- End ----------------------*/\n');
