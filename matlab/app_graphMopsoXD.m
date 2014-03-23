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
nameDb    = 'kur3';
iteration = -1;
pbest     = 1;

%-- What to do
loaddata  = 0;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%-- Objective function global parameters
nSolutions = 10000;

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
    if   iteration < 0,   it = pso.nIterations(1, 1);
    else                  it = iteration;
    end

    %----------------------------- Objective space ----------------------------%
    %-- Particles
    if pbest
        z = [ 1-pso.pPm(:, it), pso.pSz(:, it) ];
    else
        z = [ 1-pso.sPm(:, it), pso.pSz(:, it) ] ;
    end

    %-- Archive
    nFilled = archive.nFilled(it);
    zArc = [ 1-archive.sPm(1:nFilled, it), archive.sSz(1:nFilled, it) ];
else
    nObjectives = 2;
end

%------------------------------------------------------------------------------%
%----------------------------- Objective functions ----------------------------%
fprintf('/*-- Process\n');
%------------------------- Kursawe - Original -------------------------%
if strcmp(nameDb,'kur3')
    %-- Random samples with adjusted domains
    nDimensions = 3;
%     x = rand(nSolutions,nDimensions)*5.1-5; 
%     x = rand(nSolutions,nDimensions)*3.5-1.75; 
    x = rand(nSolutions,nDimensions)*10-5; 

    f1 = -10 * sum(exp(-0.2*sqrt( x(:,1:nDimensions-1).^2 + ...
                                  x(:,2:nDimensions)  .^2) ), 2);
    f2 = sum( abs(x.^0.8) + 5*( sin(x).^3 ) , 2);
    bp=1;
%--------------------------------- Deb --------------------------------%
elseif strcmp(nameDb,'deb')
    %-- Adjustments
    x1 = x1*0.9 +0.1;
    x2 = x2*0.9 +0.1;

    %-- Objective 1 ("y")
    f1       = (1-x1)/0.9;
    f(i,j,1) = f1;

    %-- Objective 2 ("x")
    f2       = ( 18.32311- ( 2 -       exp( -((x2-0.2)/0.004)^2 ) ...
                 - 0.8 * exp( -((x2-0.6)/0.4  )^2 ) ) / x1 )/17.123 ;
    f(i,j,2) = f2;
else
    fprintf('/*-- No-name is not only a brand!\n');
    return;
end

%------------------------------------------------------------------------------%
%----------------------------- Graphic parameters -----------------------------%
%-- General
graphWidth = 8; %-- centimeters
ratio      = 1/1;

%-- Marker & font
mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
sizeMkr = 7;   sizeFont = 10;   widthLine = 1;
clrMin  = 0.3;   clrMax   = 1;    pow       = 1;

%-- Figure initialization
o=1;
close all;
width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
%     posScrn  = [1300 sizeScrn(2)-104-550*(o-1) width height];
posScrn  = [530 sizeScrn(2)-104-400*(o-1) width height];
posPaper = [0 0 graphWidth graphWidth/ratio];

fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
            'PaperSize', [posPaper(3) posPaper(4)],...
            'PaperPosition', posPaper, 'Color', [1,1,1]);
ax = axes('Position', [0 0 1 1]);

hold on
%------------------------------- Objective space ------------------------------%
clr = [0 0 0];
% clr = [0.3056    0.3056    0.4253]
clr = [0.7179    0.8194    0.8194];
plot(f1, f2, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
    'MarkerSize',sizeMkr-3.5);

%------------------------------- Explored space -------------------------------%
% clr = [0.2222    0.2222    0.3108];
% 
% z = [ reshape(1-pso.sPm(:,1:it), szSwarm*it,1) ...
%       reshape(  pso.sSz(:,1:it), szSwarm*it,1) ];
% for n = size(z,1):-1:1
%    if z(n,2) > 9000
%        z(n,:) = [];
%    end
% end
% 
% pl = plot(z(:,1), z(:,2), 'o', 'Color', [0 0 0], ...
%           'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
%           'LineWidth', widthLine);
clr = bone(it);

for i = 1:it
    for n = 1:szSwarm
        z = [ 1-pso.sPm(:,i) pso.sSz(:,i) ];
    end

    for n = size(z,1):-1:1
       if z(n,2) > 9000
           z(n,:) = [];
       end
    end

    plot(z(:,1), z(:,2), 'o', 'Color', [0 0 0], ...
              'markerfacecolor', clr(i,:), 'MarkerSize',sizeMkr-3,...
              'LineWidth', widthLine);
end

%----------------------------------- Archive ----------------------------------%
clr = [0.8481    0.9028    0.9028];
pl = plot(zArc(:,1), zArc(:,2), 'o', 'Color', [0 0 0], ...
          'markerfacecolor', [1 1 1], 'MarkerSize',sizeMkr-1,...
          'LineWidth', widthLine);
hold off

set(gca, 'Position', [0.15 0.13 0.82 0.85]);

xlabel('Objective 1', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','tex');
ylabel('Objective 2', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','tex');


%------------------------------ Saving the figure -----------------------------%
nameFile = sprintf('../figures/%s_mopso', nameDb );
saveas(fig(1), nameFile, 'pdf');

fprintf('/*---------------------- End ----------------------*/\n');
