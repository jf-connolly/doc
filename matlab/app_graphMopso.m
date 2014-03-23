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
nShifts   = 1;
t         = 1;
iteration = -1;
pbest     = 1;
radius    = 0.1;

%-- What to do
loaddata  = 0;
%-- What to plot
archi       = 0;
swarm       = 0;
searchSpace = 0;

ax        = [0 65; 0.3 1; -0.6 0.5; 0 0.9];
nLevels   = 50;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%-- Objective function global parameters
nSolutions = 50000;
minD       = -5;
maxD       =  5;

%------------------------------------------------------------------------------%
%------------------------------ Reading the data ------------------------------%
if loaddata
    %--------------------------------- Reading --------------------------------%
    fprintf('/*-- Load \n');

    %-- Swarm
    nameFile = sprintf('../savedStuff/%s_%s_%dBlocks.pso', ...
                        nameDb, nameAl, nShifts);
    pso = util_readMopso(nameFile, nShifts, 1);

    nameFile = sprintf('../savedStuff/%s_%s_%dBlocks.arc', ...
                        nameDb, nameAl, nShifts);
    archive  = util_readArchiveMopso(nameFile, nShifts, 1, pso.nIterations);

    nObjectives = pso.nObjectives;
    sizeSwarm   = pso.sizeSwarm;
    
    %--------------------- Formating the data for graphics --------------------%
    subswarms = zeros(pso.sizeSwarm, nObjectives);
    lbests    = zeros(pso.sizeSwarm, nObjectives);
    flagMstrs = zeros(pso.sizeSwarm, nObjectives);
    nSS       = zeros(1,             nObjectives);

    if   iteration < 0,   it = pso.nIterations(1, t);
    else                  it = iteration;
    end

    xp = []; yp = []; zp = []; x = []; y = []; z = [];

    xp(:,:) =   pso.p  (:, 1, it, 1, t, :);
    yp(:,:) =   pso.p  (:, 2, it, 1, t, :);
    zp(:,:) = [ 1-pso.pPm(:,    it, 1, t, 1), pso.pSz(:,    it, 1, t, 2) ];


    x = [ pso.s(:, 1, it, 1, t),     pso.s(:, 1, it, 1, t) ];
    y = [ pso.s(:, 2, it, 1, t),     pso.s(:, 2, it, 1, t) ];
    z = [ 1-pso.sPm(:,    it, 1, t), pso.sSz(:,    it, 1, t) ];

    if pbest,   x = xp;  y = yp;  z = zp;   end
    

    %----------------------------- Objective space ----------------------------%
    %-- Particles
    if pbest
        z = [ 1-pso.pPm(:,    it, 1, t, 1), pso.pSz(:,    it, 1, t, 1) ;
              1-pso.pPm(:,    it, 1, t, 2), pso.pSz(:,    it, 1, t, 2) ];
    else
        z = [ 1-pso.sPm(:,    it, 1, t), pso.pSz(:,    it, 1, t) ] ;
    end

    %-- Archive
    zArc = [ 1-archive.sPm(:,pso.nIterations), archive.sSz(:,pso.nIterations) ];
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
%------------------------------- Kursawe ------------------------------%
elseif strcmp(nameDb,'kur')
    %-- Random samples with adjusted domains
    nDimensions = 2;
%     x = rand(nSolutions,nDimensions)*5.1-5; 
%     x = rand(nSolutions,nDimensions)*3.5-1.75; 
    x = rand(nSolutions,nDimensions)*10-5; 

    %-- Objective 1 ("y")
    f1 = 1-(10*exp(-0.2*sqrt( x1^2 + x2^2 )));

    %-- Objective 2 ("x")
    f2 = abs(x(1,:))^0.8 + 5*(sin(x(1,:)))^3 + ...
         abs(x(2,:))^0.8 + 5*(sin(x2))^3;

    f(i,j,1) = f1;
    f(i,j,2) = f2;
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
sizeMkr = 7;   sizeFont = 12;   widthLine = 1;
clrMin  = 0.3;   clrMax   = 1;    pow       = 1;

%-- Figure initialization
close all;
for o=1:nObjectives+1
    width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
%     posScrn  = [1300 sizeScrn(2)-104-550*(o-1) width height];
    posScrn  = [1130 sizeScrn(2)-104-400*(o-1) width height];
    posPaper = [0 0 graphWidth graphWidth/ratio];

    fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
%     ax = axes('Position', [0 0 1 1]);
end

%-------------------------------- Search space --------------------------------%
if searchSpace
    for o = 1:nObjectives

        set(0,'CurrentFigure',o);
        hold on

        %-- Figure - objective function
        [tmp ct] = contourf(X, Y, f(:,:,o), nLevels);
        set(ct, 'LineStyle', '-', 'LineWidth', 0.2, 'LineColor', [0.1 0.1 0.25]);
        colormap bone;
        box off;
        set(gca, 'XTick', [], 'YTick', []);

        %---------------------------------- Swarm ---------------------------------%
        if swarm
            %-- Which particle to plot
            graph = ones(1, pso.sizeSwarm);
            if ~pbest   for i = 1:pso.sizeSwarm   for d = 1:pso.nDimensions
                if pso.s(i, d, it, 1, t) < 0 || pso.s(i, d, it, 1, t) > 1
                    graph(i) = 0;
                end
            end,  end,  end

            %-- Figure - swarm
            gbest = pso.gbest(it, 1, t, o);

            %------------------------------ Particles -----------------------------%
            for n = 1:pso.sizeSwarm
                if graph(n)

                    %-- If master or follower, identify the master
                    if subswarms(n,o)
                        marker = mkrSubswarms(subswarms(n,o));
                    else
                        marker = mkrFree;
                    end

                    %-- Plot the graph
                    tonte      = [1 1 1] * z(n,o)^pow * (clrMax-clrMin) +clrMin;
                    sizeMarker = sizeMkr;

                    %-- Global best marker
                    if n == gbest,   mkg = marker;   end

                    pl = plot(x(n,o), y(n,o), marker, 'Color', [0 0 0], ...
                          'markerfacecolor',tonte, 'MarkerSize',sizeMarker+1,...
                          'LineWidth', widthLine);
                    if numbers
                        text( x(n,o)+0.01, y(n,o)-0.01, ...
                              ['\textbf{',int2str(n), '}'], ...
                              'Interpreter', 'latex', 'FontSize', sizeFont);
                    end
                end
            end

            %----------------------------- Global best ----------------------------%
            n  = pso.gbest(it, 1, t, o);
            xg = pso.p(n, 1, it, 1, t, o);
            yg = pso.p(n, 2, it, 1, t, o);
            tonte = [1 1 1] * z(n, o)^pow * (clrMax-clrMin) + clrMin;
            sizeMarker  = sizeMkr;
            plot(xg, yg, mkg, 'Color', tonte, ...
                'markerfacecolor', tonte, 'MarkerSize',sizeMarker);
            plot(xg,yg,[mkg, 'k'], 'MarkerSize', sizeMkr + 2, 'LineWidth', 1.5);
            plot(xg,yg,     '+k' , 'MarkerSize', sizeMkr - 3, 'LineWidth', 1.5);
            plot(xg,yg,[mkg, 'r'], 'MarkerSize', sizeMkr + 2, 'LineWidth', 0.7);
            plot(xg,yg,     '+r' , 'MarkerSize', sizeMkr - 4, 'LineWidth', 0.7);
        end

        %--------------------------------- Archive --------------------------------%
        if archi
            clr = [0.5417    0.6250    0.6667];
            xArc = archive.s(:,:, pso.nIterations, t, 1);

            for m = 1:archive.size
                if archive.filled(m, pso.nIterations, t, 1)
                    pl = plot(xArc(m,1), xArc(m,2), 'o', 'Color', [0 0 0], ...
                      'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
                      'LineWidth', widthLine);
                end
            end
        end

        set(gca, 'XLim',[0 1], 'YLim',[0 1]);
        hold off;

        %-- Saving the figure
        if pbest,    namePb  = 'p'  ;   else namePb  = ''   ;   end
        if o == 1,   nameObj = 'Clr';   else nameObj = 'Cpn';   end
        if archi || swarm 
            nameFile = sprintf('../figures/%s%s_%s%s_shift_%d', ...
                                nameDb, nameObj, nameAl, namePb, t);
        else
            nameFile = sprintf('../figures/%s%s_shift_%d', nameDb, nameObj, t);
        end

        saveas(fig(o), nameFile, 'pdf');
    end
end

%------------------------------------------------------------------------------%
%------------------------------- Objective space ------------------------------%
set(0,'CurrentFigure',nObjectives+1);
hold on

clr = [0 0 0];
% clr = [0.3056    0.3056    0.4253]
clr = [0.7179    0.8194    0.8194];
% clr = [0.2222    0.2222    0.3108];
plot(f1, f2, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
    'MarkerSize', 1);

%-- Boundaries
if archi
    x = [0, archive.boundaries'];
    for m = 1:archive.nMemetics+1
        plot([x(m) x(m)], [0 1], 'k--','LineWidth', widthLine);
    end
end
%---------------------------------- Particles ---------------------------------%
if swarm
    clr = [0.5417    0.6250    0.6667]
    pl = plot(z(:,2), z(:,1), 'o', 'Color', [0 0 0], ...
          'markerfacecolor', clr, 'MarkerSize',sizeMarker-1,...
          'LineWidth', widthLine);

    if numbers
        for n = 1:sizeSwarm
            text( z(n,2)+0.006, z(n,1)-0.006, int2str(n), ...
                  'FontName', 'Times New Roman', 'Interpreter', 'latex', ...
                  'FontSize', sizeFont-3, 'FontWeight','Bold' );
        end
        for n = sizeSwarm+1:sizeSwarm*2
            np = n-sizeSwarm;
            text( z(np,2)+0.006, z(np,1)-0.006, int2str(np), ...
                  'FontName', 'Times New Roman', 'Interpreter', 'latex', ...
                  'FontSize', sizeFont-3, 'FontWeight','Bold' );
        end
    end
end

%----------------------------------- Archive ----------------------------------%
if archi
    clr = [0.8481    0.9028    0.9028];
    for n = 1:archive.size
        if archive.filled(n, pso.nIterations, t, 1)
            pl = plot(zArc(n,2), zArc(n,1), 'o', 'Color', [0 0 0], ...
              'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
              'LineWidth', widthLine);
        end
    end
end
hold off

set(gca, 'Position', [0.08 0.08 0.91 0.91], 'FontName', 'Times New Roman', ...
    'XLim', [-20 -5], 'YLim', [-10 25], 'XTick', [], 'YTick', [] );



xlabel('$f_2(\textbf{h})$', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','latex');
ylabel('$f_1(\textbf{h})$', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','latex');


%------------------------------ Saving the figure -----------------------------%
if pbest,    namePb  = 'p'  ;   else namePb  = ''   ;   end
if o == 1,   nameObj = 'Clr';   else nameObj = 'Cpn';   end
if archi || swarm 
    nameFile = sprintf('../figures/%sPar_%s%s_shift_%d_%dmemetics', ...
                        nameDb, nameAl, namePb, t, archive.nMemetics(it) );
else
    nameFile = sprintf('../figures/%sPar_shift_%d', nameDb, t );
end
saveas(fig(o), nameFile, 'pdf');
% print(fig(o), '-djpeg', '-r600', nameFile);

fprintf('/*---------------------- End ----------------------*/\n');
