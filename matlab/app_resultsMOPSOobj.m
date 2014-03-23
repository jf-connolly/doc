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
nameDb        = 'kur3';
nameAl        = 'multi';
nReplications = 2;
nObjectives   = 2;

t           = 1;
ax          = [0 65; 0.3 1; -0.6 0.5; 0 0.9];

objectives  = 0;
paretoFront = 0;
maximiz     = 0;
savePF      = 0;
loadPF      = 1;
indicators  = 1;
graph       = 0;

fprintf('\n/*--------------------- Start ---------------------*/\n');

%-- Objective function global parameters
nSolutions = 10;

%----------------------------- Objective functions ----------------------------%
if objectives
    nEval = 0;
    %--------------------------- Kursawe --------------------------%
    if strcmp(nameDb,'kur3')
        %-- Random samples with adjusted domains
        nDimensions = 3;
%         x           = rand(nSolutions,nDimensions)*5.1-5; 
%         x           = rand(nSolutions,nDimensions)*3.5-1.75; 
        x           = rand(nSolutions,nDimensions)*10-5; 

        f1 = -10 * sum(exp(-0.2*sqrt( x(:,1:nDimensions-1).^2 + ...
                                      x(:,2:nDimensions)  .^2) ), 2);
        f2 = sum( abs(x.^0.8) + 5*( sin(x).^3 ) , 2);
        bp=1;
        
    %---------------------------------- ZDT3 ----------------------------------%
    elseif strcmp(nameDb,'zd3')
        %-- Random samples
        nDimensions = 30;
        x1          = rand(nSolutions,1);
%         xOthers     = rand(nSolutions,nDimensions-1);
        xOthers     = sin(rand(nSolutions,nDimensions-1));
        
        f1 = x1;
    
        g  = 1 + 9/(nDimensions-1) + sum(xOthers, 2);
        h  = 1 - sqrt(f1./g) - (f1./g) .* sin(10*pi*f1);
        f2 = g.*h;
        
        x = [x1, xOthers];

    %---------------------------------- ZDT4 ----------------------------------%
    elseif strcmp(nameDb,'zd4')
        nDimensions = 10;                                   %- Nb. of dimensions
        x1          = rand(nSolutions,1);                   %- Random samples
        xOthers     = rand(nSolutions,nDimensions-1)*10 -5; %- Adjust domains
        
        f1 = x1;
    
%         g  = 1 + 10*(nDimensions-1) + sum(xOthers.^2 - 10*cos(4*pi*xOthers),2);
        g = zeros(nSolutions,1);
        for d = 1:nDimensions-1
           g = g + xOthers(:,d).^2 - 10*cos(4*pi*xOthers(:,d));
        end
        g = g + 1 + 10*(nDimensions-1);
        h  = 1 - sqrt(f1./g);
        f2 = g.*h;

        x = [x1, xOthers];
        
%         %-- Objective 1 ("y")
%         f1 = x1;
% 
%         %-- Objective 1 ("y")
%         g  = 1 + 10 + x2^2 - 10*cos(4*pi*x2);
%         h  = 1 - sqrt(f1/g);
%         f2 = g * h;
    else
        fprintf('/*-- No-name is not only a brand!\n');
        return;
    end
    
    [min(f1) max(f1); min(f2) max(f2)]
end

%-------------------------------- Pareto Front --------------------------------%
if paretoFront
    nElements = 0;
    front = zeros(nSolutions, 2);
    frons = zeros(nSolutions, nDimensions);
    for n = 1:nSolutions;
        tested       = [f1(n) f2(n)];
        testex       = x(n,:);
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
                    frons(m, :) = [];
                end
            end
        else
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
        end
        
        if nonDominated
            nElements           = nElements + 1;
            front(nElements, :) = tested;
            frons(nElements, :) = testex;
        end
    end
    
    front( nElements+1:size(front,1),: ) = [];
    frons( nElements+1:size(frons,1),: ) = [];
end

%--------------------------------- Save'n'load --------------------------------%
if savePF
    nameFile = sprintf('../savedStuff/kursawe3d_pf.mat');
    save(nameFile, 'front', 'frons', 'x', 'f1', 'f2', '-mat');
end

if loadPF
    nameFile = sprintf('../savedStuff/kursawe3d_pf.mat');
    load(nameFile, 'front', 'frons', 'x', 'f1', 'f2', '-mat');
end

%------------------------------------------------------------------------------%
%---------------------------- Performance indicator ---------------------------%
%-- Load the archive 
if indicators
    nameFile = sprintf('../savedStuff/%s_%s.pso', nameDb, nameAl);
    pso      = util_readPsoMulti(nameFile, 1, nReplications);
    
    nameFile = sprintf('../savedStuff/%s_%s.arc', nameDb, nameAl);
    archive  = util_readArchive(nameFile, 1, nReplications, pso.nIterations);
    
    sizeArc     = archive.size;

    %-- Indicators
    gd = zeros(1,nReplications);
    sp = zeros(1,nReplications);
    er = zeros(1,nReplications);

    for r = 1:nReplications
        %-- Preprocessing
        nIts     = pso.nIterations(r);
        nMembers = archive.nMembers(1, nIts, 1, r);
        zArc     = [ archive.sPm(:,nIts, 1, r), archive.sSz(:,nIts, 1, r) ];
        sArc     = [ archive.s(:, :, nIts, 1, r) ];

        zRed = zArc;
        sRed = sArc;
        for m = sizeArc:-1:1
            if ~archive.filled(m, nIts, 1, r)
                zRed(m,:) = [];
                sRed(m,:) = [];
            end
        end
        
        nElements = size(front,1);

        %-- Generational distance
        gd(r) = 0;
        for m = 1:nMembers
            zTested   = repmat(zRed(m,:), nElements, 1);
            distances = sum((zTested-front).^2, 2);

            gd(r) = gd(r) + min(distances);
        end

        gd(r) = sqrt( gd(r) ) / nMembers;

        %-- Spacing distance
        di = zeros(1,nMembers);
        i  = 0;
        for m = 1:nMembers
            zTested    = repmat(zRed(m,:), nMembers-1, 1);
            zRed2      = zRed;
            zRed2(m,:) = [];
            distances  = sum( abs(zTested-zRed2), 2 );
            di(m)      = min(distances);
        end
%-- Indicators
    gd = zeros(1,nReplications);
    sp = zeros(1,nReplications);
    er = zeros(1,nReplications);

    for r = 1:nReplications
        %-- Preprocessing
        nIts     = pso.nIterations(r);
        nMembers = archive.nMembers(1, nIts, 1, r);
        zArc     = [ 1-archive.sPm(:,nIts, 1, r), archive.sSz(:,nIts, 1, r) ];
        sArc     = [ archive.s(:, :, nIts, 1, r) ];

        zRed = zArc;
        sRed = sArc;
        for m = sizeArc:-1:1
            if ~archive.filled(m, nIts, 1, r)
                zRed(m,:) = [];
                sRed(m,:) = [];
            end
        end
        
        zRed(:,1) = zRed(:,1) - 1;

        %-- Generational distance
        gd(r) = 0;
        for m = 1:nMembers
            zTested   = repmat(zRed(m,:), nElements, 1);
            distances = sum((zTested-front).^2, 2);

            gd(r) = gd(r) + min(distances);
        end

        gd(r) = sqrt( gd(r) ) / nMembers;

        %-- Spacing distance
        di = zeros(1,nMembers);
        i  = 0;
        for m = 1:nMembers
            zTested   = repmat(zRed(m,:), nMembers-1, 1);
            zRed2     = zRed;
            zRed2(m,:) = [];
            distances = sum( abs(zTested-zRed2), 2);
            di(m)     = min(distances);
        end

        sp(r) = sqrt( sum( (mean(di)-di).^2 ) / (nMembers-1) );

        %-- Error ratio
        ei = zeros(1,nMembers);
        for m = 1:nMembers
            zTested   = repmat(zRed(m,:), nElements, 1);
            distances = sum((zTested-front).^2, 2);

            if min(distances) > 0.001,   ei(m) = 1;   end
        end

        er(r) = sum(ei)/nMembers;
    end

    confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
    pValue       = norminv(confInterval) / sqrt(nReplications);

    fprintf('/*--\n   gd        sp        er\n');
    moResults(1,:) = [mean(gd) mean(sp) mean(er)];
    moResults(2,:) = [std(gd) std(sp) std(er)] * pValue;
    moResults

        sp(r) = sqrt( sum( (mean(di)-di).^2 ) / (nMembers-1) );

        %-- Error ratio
        ei = zeros(1,nMembers);
        for m = 1:nMembers
            zTested   = repmat(zRed(m,:), nElements, 1);
            distances = sum((zTested-front).^2, 2);

            if min(distances) > 0.001,   ei(m) = 1;   end
        end

        er(r) = sum(ei)/nMembers;
    end

    confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
    pValue       = norminv(confInterval) / sqrt(nReplications);

    fprintf('/*--\n   gd        sp        er\n');
    moResults(1,:) = [mean(gd) mean(sp) mean(er)];
    moResults(2,:) = [std(gd) std(sp) std(er)] * pValue;
    moResults
end
% 
%------------------------------------------------------------------------------%
%----------------------------- Graphic parameters -----------------------------%
%-- General
graphWidth = 8;
ratio      = 1/1;

%-- Marker & font
mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
sizeMkr = 7;   sizeFont = 10;   widthLine = 1;
clrMin  = 0;   clrMax   = 1;    pow        = 1;

%-- Figure initialization
% close all;

width  = res*graphWidth+2;   height = res*graphWidth/ratio+2;
posScrn = [1300 sizeScrn(2)-104-550 width height];
posPaper = [0 0 graphWidth graphWidth/ratio];

fig(1) = figure(1);   clf(fig(1));   set(fig(1),'Color',[1,1,1]);
set(fig(1), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
            'PaperSize', [posPaper(3) posPaper(4)],...
            'PaperPosition', posPaper, 'Color', [1,1,1]);
ax = axes('Position', [0 0 1 1]);


%---------------------------- Pareto Optimal front ----------------------------%
if graph
    fprintf('/*-- Objective space\n');

    hold on 
        %-- Objective space
        clr = [0.7179    0.8194    0.8194];
        plot(f1, f2, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
            'MarkerSize', 1);

        %-- Pareto front
%         clr = [0.3056    0.3056    0.4253];
%         plot(front(:,1), front(:,2), 'o', 'Color', [0 0 0], ...
%             'markerfacecolor', clr, 'MarkerSize',sizeMkr-1, 'LineWidth', widthLine);

%     %-- Objective space
%     clr = [0.7179    0.8194    0.8194];
%     plot(f2, f1, 'ko', 'Color', clr, 'markerfacecolor', clr, ...
%         'MarkerSize',sizeMkr-3.5);
% 
    %-- Pareto front
    clr = [0.3056    0.3056    0.4253];
    [temp index] = sort(front(:,2));
    sFront       = front(index,:);
    plot(sFront(:,1), sFront(:,2), 'o', 'Color', [0 0 0], ...
        'markerfacecolor', clr, 'MarkerSize', 2, 'LineWidth', widthLine);

%     %-- Extra points
%     extra1 = f(1:4,32:36,1);   extra2 = f(1:4,32:36,2);
%     extra1 = extra1(:);   extra2 = extra2(:);
% 
%     plot(extra2, extra1, 'o', ...
%          'Color', [0 0 0], 'markerfacecolor', [1 1 1],    ...
%          'MarkerSize',sizeMkr-1, 'LineWidth', widthLine);

    %-- Archive
    if indicators
        clr = [1 1 1];
        plot(zRed(:,1), zRed(:,2), 's', 'Color', [0.7 0.7 0], ...
                  'markerfacecolor', clr, 'MarkerSize',sizeMkr-1,...
                  'LineWidth', widthLine);
    end
    hold off

%     set(gca, 'Position', [0.15 0.13 0.82 0.85], 'XLim', [0 1], 'YLim', [0 1]);
%     set(gca, 'Position', [0.15 0.13 0.82 0.85]);
    set(gca, 'Position', [0.15 0.13 0.82 0.85]);
%     set(gca, 'XLim', [-20 -5], 'YLim', [-11 25.5]);
    set(gca, 'XLim', [-20 -12.5], 'YLim', [-11 0.5]);

    xlabel('Objective 1', 'FontSize', sizeFont, 'Interpreter','latex');
    ylabel('Objective 2', 'FontSize', sizeFont, 'Interpreter','latex');

    nameFile = sprintf('../figures/%sTheoricPar', nameDb);
    saveas(fig(1), nameFile, 'pdf');
end
fprintf('/*---------------------- End ----------------------*/\n');
