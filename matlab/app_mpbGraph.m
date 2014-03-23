%-- DEFINE
RIEN = 1;     BAR = 2;     NORMAL = 3;
temp = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res = get(0,'ScreenPixelsPerInch')/2.56;
%-- END DEFINE

%------------------------------------------------------------------------------%
%---------------------------- Param�tres variables ----------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
nameDb  = 'mpb';
nShifts = 4;
t       = 4;
pbest   = 1;
nameHp  = 'Hdnc';

D       = 4;
bar     = 1;
ax      = [0 65; 0.3 1; -0.6 0.5; 0 0.9];
nDeniv  = 10;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%--------------------- Choix des param�tres du graphique ----------------------%
%-- Marker & font
mkrSubswarms = 'os^dpv><'; mkrFree = '*';    sizeMkr = 7;
clrMin       = 0;          clrMax     = 1;   pow     = 1;
sizeFont     = 10;         witdthLine = 1;

%-- General
graphWidth     = 8;
ratio          = 1/1;

%-- Figure initialization
close all;
width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
posScrn = [2400 sizeScrn(2)-104 width height];
posPaper = [0 0 graphWidth graphWidth/ratio];

fig(1) = figure(1);   clf(fig(1));   set(fig(1),'Color',[1,1,1]);
set(fig(1), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
            'PaperSize', [posPaper(3) posPaper(4)],...
            'PaperPosition', posPaper, 'Color', [1,1,1]);
ax = axes('Position', [0 0 1 1]);

clear frames;

%-- Objective function global parameters
resolution = 50;
xObj       = 0 :1/resolution :1;
[Y X]      = meshgrid(xObj,xObj);

%------------------------------ Reading the data ------------------------------%
fprintf('/*-- Load \n');

%-- Objective function - initialization
nameFile = sprintf('../database/%s/%dshifts/shift_1.db', nameDb, nShifts);
dbase    = util_readHeader(nameFile);

nPeaks = dbase.nPatterns;
nDims  = dbase.nFeatures - 2;

posPeaks = zeros(nPeaks, nDims);
height   = zeros(nPeaks);
width    = zeros(nPeaks);

%-- Objective function - peaks parameters
nameFile = sprintf('../database/%s/%dshifts/shift_%d.db',nameDb,nShifts,t);
dbase    = util_readDb(nameFile);

for p = 1:nPeaks
    posPeaks(p,:) = dbase.data(p, 1:nDims);
    height(p)     = dbase.data(p, nDims +1);
    width (p)     = dbase.data(p, nDims +2);
end

%-- Swarm
nameFile = sprintf('../savedStuff/%s_%s_%dBlocs_2D.pso', ...
                  nameDb, nameHp, nShifts);
pso = util_readPso(nameFile, nShifts, 1);

%------------------------------------------------------------------------------%
%-------------------------------- Movie frames --------------------------------%
fprintf('/*-- Process\n');

currentframe = 0;
subswarms = zeros(pso.sizeSwarm, pso.nIterations(t));
masters   = zeros(10,1);
flagMstrs = zeros(10,1);

nSS       = 0;

%--------------------------- Objective function ---------------------------%
%-- Initialization
f = zeros(resolution, resolution);
normS = zeros(1, nPeaks);
fpeak = zeros(1, nPeaks);

%-- Objective function
for i = 1: size(xObj,2);
    for j = 1: size(xObj,2);
        posCurr = [xObj(i) xObj(j)];

        for p = 1:nPeaks
            posPeaks(p,:);
            normS(p) = sum( (posCurr - posPeaks(p,:)).^2 );
            fpeak(p) = height(p) * exp( - normS(p) / width(p)^2 );
        end
        f(i,j) = max(fpeak);   %-- Because maximization
    end
end

%-- Figure - objective function
[tmp c(t)] = contourf(X,Y,f, 50);
set(c(t), 'LineStyle', 'none');
colormap bone;
box off;
set(gca, 'XTick',[], 'YTick', []);

%-- Swarm's position and fitness
x = [];   y = [];  z = [];

nIt = pso.nIterations(1, t);
if pbest
    x(:,:) = pso.p  (:, 1, 1:nIt, 1, t);
    y(:,:) = pso.p  (:, 2, 1:nIt, 1, t);
    z(:,:) = pso.pPm(:,    1:nIt, 1, t);
else
    x(:,:) = pso.s  (:, 1, 1:nIt, 1, t);
    y(:,:) = pso.s  (:, 2, 1:nIt, 1, t);
    z(:,:) = pso.sPm(:,    1:nIt, 1, t);
end

%-- Normalizing z (for the color)
z = (z - min(min(z))) / (max(max(z)) - min(min(z)));

%---------------- Identify the subswarms & free particles -----------------%
it = pso.nIterations(t);

%-- Identify the masters
for p = 1:pso.sizeSwarm

    %-- For masters & followers only
    if pso.lbest(p,it,1,t) > 0

        %-- Check for existing masters
        flagMaster = 1;
        for i = 1:nSS
            if pso.lbest(p,it,1,t) == masters(i),   flagMaster = 0;   end
        end

        %-- New masters
        if flagMaster
            nSS = nSS +1;
            masters(nSS) = pso.lbest(p,it,1,t);
        end
    end
end

%-- Assigne the subswarms
for p = 1:pso.sizeSwarm

    %-- Masters & followers
    if pso.lbest(p,it,1,t) > 0

        for i = 1:nSS
            if pso.lbest(p,it,1,t) == masters(i);
                subswarms(p,it) = i;
            end
        end
    %-- Free particles    
    else
        subswarms(p,it) = 0;
    end
end

%--------------------------------- Swarm ----------------------------------%
fprintf('/*-- Figure (time = %d)\n', t);

%-- Which particle to plot
graph = ones(1, pso.sizeSwarm);
if ~pbest   for i = 1:pso.sizeSwarm   for d = 1:pso.nDimensions
    if pso.s(i, d, it, 1, t) < 0 || pso.s(i, d, it, 1, t) > 1
        graph(i) = 0;
    end
end,  end,  end

%-- Figure - objective function
[tmp c(t)] = contourf(X,Y,f, 50);
set(c(t), 'LineStyle', 'none');
colormap bone;
box off;
set(gca, 'XTick', [], 'YTick', []);

%-- Figure - swarm
gbest = pso.gbest(it, 1, t);
hold on
    for n = 1:pso.sizeSwarm
        if graph(n)

            %-- If master or follower, identify the master
            if subswarms(n, it)
                marker = mkrSubswarms(subswarms(n, it));
            else
                marker = mkrFree;
            end

            %-- Plot the graph
            tonte =[z(n,it) z(n,it) z(n,it)] .^ pow ...
                      * (clrMax-clrMin) + clrMin;
            sizeMarker  = sizeMkr;%z(n,it) ^ pow * sizeMkr + 1;

            %-- Global best marker
            if n == gbest,   mkg = marker;   end
            
            pl = plot(x(n,it), y(n,it), marker, 'Color', [0 0 0], ...
                  'markerfacecolor', tonte, 'MarkerSize',sizeMarker+1, ...
                  'LineWidth', witdthLine);
%             text( x(n,it)+0.006, y(n,it)-0.006, int2str(n), ...
%                   'FontName', 'Times New Roman',...
%                   'FontSize', sizeFont-3, 'FontWeight','Bold');
        end
    end
    
    %-- Area around masters
    radius = 0.1;
    for n = 1:nSS
        lb(n) = annotation('ellipse', [0.5, 0.5,  0.6, 0.6], 'LineStyle', ':');
                 
        set(lb(n), 'Position', ...
                   [x(masters(n),it) - radius, y(masters(n),it) - radius, ...
                    2*radius, 2*radius] );
    end
                    


    %-- Global best
    n  = pso.gbest(it, 1, t);
    xg = pso.p(n, 1, it, 1, t);
    yg = pso.p(n, 2, it, 1, t);
    tonte = [z(n,it) z(n,it) z(n,it)].^(pow) * (clrMax-clrMin) + clrMin;
    sizeMarker  = z(n)^pow*sizeMkr + 1;
    plot(xg, yg, mkg, 'Color', tonte, ...
        'markerfacecolor', tonte, 'MarkerSize',sizeMarker);
    plot(xg,yg, [mkg, 'k'], 'MarkerSize', sizeMkr + 2, 'LineWidth', 1.5);
    plot(xg,yg,      '+k' , 'MarkerSize', sizeMkr - 3, 'LineWidth', 1.5);
    plot(xg,yg, [mkg, 'r'], 'MarkerSize', sizeMkr + 2, 'LineWidth', 0.7);
    plot(xg,yg,      '+r' , 'MarkerSize', sizeMkr - 4, 'LineWidth', 0.7);

    set(gca, 'XLim',[0 1], 'YLim',[0 1]);
hold off

%-- Saving the figure
nameFile = sprintf('../figures/mpb_shift_%d', t);
saveas(fig(1), nameFile, 'pdf');

fprintf('/*---------------------- End ----------------------*/\n');
