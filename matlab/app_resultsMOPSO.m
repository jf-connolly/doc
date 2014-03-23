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
nameDb      = 'cnrc64';
nameAl      = 'mopso';
nShifts     = 1;
nMemetics   = 1;
sizeHood    = 6;
nIterations = 100;
nObjectives = 2;

t           = 1;
ax          = [0 65; 0.3 1; -0.6 0.5; 0 0.9];
nLevels     = 40;

loadFile = 0;
compute  = 1;
maximiz  = 0;
savePF   = 0;
loadPF   = 0;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%-- Objective function global parameters
resolution = 80;

%------------------------------ Load the archive ------------------------------%
if loadFile
    nameFile = sprintf('../savedStuff/%s_%s_%dBlocks_%dmem_%dnd.arc', ...
                        nameDb, nameAl, nShifts, nMemetics, sizeHood);
    archive  = util_readArchive(nameFile, nShifts, 1, nIterations);
    
    sizeArc  = archive.size;
    nMembers = archive.nMembers(nIterations);
    zArc = [ archive.sPm(:,nIterations), archive.sCn(:,nIterations) ];
    sArc = [ archive.s(:,:,nIterations) ];
    
    
    zRed = zArc;
    sRed = sArc;
    for m = sizeArc:-1:1
        if ~archive.filled(m,nIterations)
            zRed(m,:) = [];
            sRed(m,:) = [];
        end
    end

end

%----------------------------- Objective functions ----------------------------%
if compute
    xObj  = 0 :1/resolution :1;
    [Y X] = meshgrid(xObj,xObj);
%     [Y X] = meshgrid(( xObj*(0.55-0.3) ) + 0.3, ( xObj*(0.55-0.3) ) + 0.3);
    f     = zeros(resolution, resolution, nObjectives);

    %-- Load multipeaks benchmark problem - if needed
    if strcmp(nameDb, 'mpb')

        %-- Objective function - initialization
        nameFile = sprintf('../database/mpbClr/%dshifts/shift_1.db', nShifts);
        dbase    = util_readHeader(nameFile);

        nPeaks = dbase.nPatterns - 1;
        nDims  = dbase.nFeatures - 2;

        posPeaks = zeros(nPeaks, nDims, nObjectives);
        heights  = zeros(nPeaks,        nObjectives);
        widths   = zeros(nPeaks,        nObjectives);
        minf     = zeros(1,             nObjectives);
        maxf     = zeros(1,             nObjectives);

        normS = zeros(1, nPeaks);
        fpeak = zeros(1, nPeaks);

        %-- Objective function - peaks parameters
        nameFile = sprintf('../database/mpbClr/%dshifts/shift_1.db', nShifts);
        dbase(1) = util_readDb(nameFile);

        nameFile = sprintf('../database/mpbCpn/%dshifts/shift_1.db', nShifts);
        dbase(2) = util_readDb(nameFile);

        for o = 1:nObjectives
            for p = 1:nPeaks
                posPeaks(p,:,o) = dbase(o).data(p, 1:nDims);
                heights (p,  o) = dbase(o).data(p, nDims +1);
                widths  (p,  o) = dbase(o).data(p, nDims +2);
            end
            minf(o) = dbase(o).data(nPeaks+1, 1);
            maxf(o) = dbase(o).data(nPeaks+1, 2);
        end
    end

    %-- Objective function
%     x = zeros(resolution+1, resolution+1, 2);
    x = zeros(resolution+1, resolution+1, resolution+1, 3);
    for i = 1: size(xObj,2);
        for j = 1: size(xObj,2);

%             x(i,j,1) = ( xObj(i)*(0.55-0.3) ) + 0.3;
%             x(i,j,2) = ( xObj(j)*(0.55-0.3) ) + 0.3;
            x(i,j,1) = xObj(i);
            x(i,j,2) = xObj(j);

            %----------------------- Multipeak benchmark ----------------------%
            if strcmp(nameDb,'mpb')
                x1 = x(i,j,1);   x2 = x(i,j,2);
                
                for o = 1: nObjectives
                    for p = 1:nPeaks
                        posPeaks(p,:,o);
                        normS(p) = sum( ([x1 x2] - posPeaks(p,:,o)).^2 );
                        fpeak(p) = heights(p,o) * exp(- normS(p)/widths(p,o)^2);
                    end
                    f(i,j,o) = (max(fpeak) - minf(o)) / (maxf(o)- minf(o));


                end
            %----------------------------- Kursawe ----------------------------%
            elseif strcmp(nameDb,'kur')
                %-- Adjustments with min-max
                x1 = x(i,j,1)*10 - 5;   x2 = x(i,j,2)*10 - 5;
%                 x1 = x(i,j,1)*2.5 - 2.001;  x2 = x(i,j,2)*2.5 - 2.002;

                %-- Objective 1 ("y")
                f1 = 1-( abs(x1)^0.8 + 5*(sin(x1)^3) + ...
                         abs(x2)^0.8 + 5*(sin(x2)^3) + ...
                         7.1656 ) / (16.9352 + 7.1656);

                %-- Objective 2 ("x")
                f2 = 1-(10*exp(-0.2*sqrt( x1^2 + x2^2 )) -2.4312) / (10-2.4312);
                
                f(i,j,1) = f1;
                f(i,j,2) = f2;

            %----------------------------- Deb - 3 ----------------------------%
            elseif strcmp(nameDb,'de3')
                %-- Adjustments
                x1 = x(i,j,1);   x2 = x(i,j,2)*60 -30;

                %-- Objective 1 ("y")

                g = 11 + x2^2 - 10*cos(2*pi*x2);
                h = 0;   if f1 <= g,   h = 1-sqrt(x1/g);   end
                
                f1 = x1;        f(i,j,1) = 1-f1;
                f2 = g*h;       f(i,j,2) = 1-f2/901;
            %----------------------------- Deb - 4 ----------------------------%
            elseif strcmp(nameDb,'de4')
                %-- Adjustments
                x1 = x(i,j,1)*0.9 + 0.1;   x2 = x(i,j,2)*0.9 + 0.1;

                %-- Objective 1 ("y")
                f1       = (1-x1)/0.9;
                f(i,j,1) = f1;

                %-- Objective 2 ("x")
                 f2       = ( 18.32311- ( 2 -       exp( -((x2-0.2)/0.004)^2 ) ...
                     - 0.8 * exp( -((x2-0.6)/0.4  )^2 ) ) / x1 )/17.123;
                f(i,j,2) = f2;
                
                
            %------------------------------ ZDT3 ------------------------------%
            elseif strcmp(nameDb,'zd3')
                x1 = x(i,j,1);   x2 = x(i,j,2);
                
                %-- Objective 1 ("y")
                f1 = x1;
                
                %-- Objective 1 ("y")
                g  = 1 + 9/3 + x2;
                h  = 1 - sqrt(f1/g) - (f1/g) * sin(10*pi * f1);
                f2 = g * h;
                
                f(i,j,1) = 1-f1;
                f(i,j,2) = 1-(f2-1.3061)/(5-1.3061);
                f(i,j,1) = f1;
                f(i,j,2) = f2;
            %------------------------------ ZDT4 ------------------------------%
            elseif strcmp(nameDb,'zd4')
                x1 = x(i,j,1);   x2 = x(i,j,2)*10 - 5;
                
                %-- Objective 1 ("y")
                f1 = x1;
                
                %-- Objective 1 ("y")
                g  = 1 + 10 + x2^2 - 10*cos(4*pi*x2);
                h  = 1 - sqrt(f1/g);
                f2 = g * h;
                
%                 f(i,j,1) = 1-f1;
%                 f(i,j,2) = 1-f2/43.5625;
                f(i,j,1) = f1;
                f(i,j,2) = f2;
            else
                fprintf('/*-- No-name is not only a brand!\n');
                return;
            end
            
            
        end
    end

    x1 = x(:,:,1);  x1 = x1(:);
    x2 = x(:,:,2);  x2 = x2(:);
    f1 = f(:,:,1);  f1 = f1(:);
    f2 = f(:,:,2);  f2 = f2(:);
    [min(f1) max(f1); min(f2) max(f2)]

    %------------------------------ Pareto Front ------------------------------%
    front = zeros(resolution^2,2);
    frons = zeros(resolution^2,2);
    nElements = 0;
    nSolutions = resolution^2;
    for n = 1:nSolutions;
        tested       = [f1(n) f2(n)];
        testes       = [x1(n) x2(n)];
        nonDominated = 1;

        if maximiz
            for m = nElements:-1:1
                if tested(1) == front(m,1) && tested(2) == front(m,2)
                    nonDominated = 1;
                elseif tested(1) <= front(m,1) && tested(2) <= front(m,2)
                    nonDominated = 0;
                elseif tested(1) >= front(m,1) && tested(2) >= front(m,2)
                    nElements = nElements - 1;
                    front(m, :) = [];
                end
            end
        else
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
if savePF
    nameFile = sprintf('../savedStuff/%s_pf.mat', nameDb);
    save(nameFile, 'front', 'frons', 'x', 'f', 'f1', 'f2', '-mat');
end

if loadPF
    nameFile = sprintf('../savedStuff/%s_pf.mat', nameDb);
    load(nameFile, 'front', 'frons', 'x', 'f', 'f1', 'f2', '-mat');
end

%------------------------------------------------------------------------------%
%---------------------------- Performance indicator ---------------------------%
% %-- Generational distance
% gd = 0;
% for m = 1:nMembers
%     zTested   = repmat(zRed(m,:), nElements, 1);
%     distances = sum((zTested-front).^2, 2);
% 
%     gd = gd + min(distances);
% end
% 
% gd = sqrt(gd) / nMembers;
% 
% %-- Spacing distance
% di = zeros(1,nMembers);
% i  = 0;
% for m = 1:nMembers
%     zTested   = repmat(zRed(m,:), nMembers-1, 1);
%     zRed2     = zRed;
%     zRed2(m,:) = [];
%     distances = sum( abs(zTested-zRed2), 2);
%     di(m)     = min(distances);
% end
% 
% sp = sqrt( sum( (mean(di)-di).^2 ) / (nMembers-1) );
% 
% %-- Error ratio
% ei = zeros(1,nMembers);
% for m = 1:nMembers
%     zTested   = repmat(zRed(m,:), nElements, 1);
%     distances = sum((zTested-front).^2, 2);
% 
%     if min(distances) > 0.001,   ei(m) = 1;   end
% end
% 
% er = sum(ei)/nMembers;
% 
% fprintf('/*-- gd    sp    er\n');
% [gd sp er]
% 
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

for o=1:nObjectives+1
    width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
    posPaper = [0 0 graphWidth graphWidth/ratio];

    fig(o) = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
    ax = axes('Position', [0 0 1 1]);
end


%-------------------------------- Search space --------------------------------%
fprintf('/*-- Search space\n');
for o = 1:nObjectives
    
    set(0,'CurrentFigure', o);

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
             'MarkerSize',sizeMkr-3, 'LineWidth', widthLine);

        plot(frons(:,2), frons(:,1), 'o', ...
             'markeredgecolor', clr-0.3, 'markerfacecolor', clr,        ...
             'MarkerSize',sizeMkr-3, 'LineWidth', widthLine);
         
    hold off
    
    box off;
    set(gca, 'XTick', [], 'YTick', []);
    
    nameFile = sprintf('../figures/%s_%d_TheoricPar_%s', nameDb, o, nameAl);
                  
    saveas(fig(o), nameFile, 'pdf');
end

%---------------------------- Pareto Optimal front ----------------------------%
fprintf('/*-- Objective space\n');
set(0,'CurrentFigure', nObjectives+1);

hold on 

    %-- Objective space
    clr = [0.7179    0.8194    0.8194];
    plot(f1, f2, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
        'MarkerSize',sizeMkr-3.5);
    
    %-- Pareto front
    clr = [0.3056    0.3056    0.4253];
    plot(front(:,1), front(:,2), 'o', 'Color', [0 0 0], ...
        'markerfacecolor', clr, 'MarkerSize',sizeMkr-1, 'LineWidth', widthLine);
hold off

set(gca, 'Position', [0.15 0.13 0.82 0.85], 'XLim', [0 1], 'YLim', [0 1]);

xlabel('Objective 1', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','tex');
ylabel('Objective 2', 'FontName', 'Times New Roman',...
 'FontSize', sizeFont, 'Interpreter','tex');

nameFile = sprintf('../figures/%sTheoricParFront', nameDb);
saveas(fig(o+1), nameFile, 'pdf');

fprintf('/*---------------------- End ----------------------*/\n');
