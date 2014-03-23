%------------------------------------------------------------------------------%
%-- Results: batch -vs- inremental for different settings
%-- Works with util_getName
%------------------------------------------------------------------------------%
fprintf('\n/*---------------- Results - static ----------------*/\n');

%-- DEFINE
elsevier = 1; ieee = 2;    lnsc = 3;
classRate = 1; compression = 2; sizeEnsemble = 3;

temp = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res = get(0,'ScreenPixelsPerInch')/2.56;
%-- END DEFINE

%------------------------------------------------------------------------------%
%--------------------------- User-defined parameters --------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
nameDb         = 'cnrc64';%  {cis, p2, cnrc64, mobof64}
nameSc         = 'add';%    {add, upd}
nBlocks        =  10;
nameTest       = 'total';%   {'Ens','vid','tst'}
nTest          =  4;
legLocation    = 'EastOutside';
legOrientation = 'vertical';
nReplications  =  50;

compute    = 1;
saveResult = 0;
loadResult = 0;
graph      = 0;
show       = 1;
hypTest    = 0;

%-- Graphs
errG = 1;
cpnG = 1;
ensG = 1;

%-- Graph width & ratio
typeDoc  = elsevier; %  elsevier;
ratioErr = 3; %  {1.75, 2.5};
ratioOtr = 2;

%---------------------------- Graphics parameters -----------------------------%
%-- Graphs width
switch typeDoc
    case elsevier
        errWidth = 16.64;
    case ieee
        errWidth = 8.96;
    case lnsc
        errWidth = 12.2;
end
otherWidth = errWidth/2.1;

%------------------------------------------------------------------------------%
%------------------------------- Initialization -------------------------------%
close all;
widthErr  = res*errWidth;     heightErr = res*errWidth/ratioErr;
widthOtr  = res*otherWidth;   heightOtr = res*otherWidth/ratioOtr;

if graph
    %-- Classification rate 
    fig_1 = figure(1);     clf(fig_1);
    posScrn = [0 sizeScrn(2)-heightErr-104 widthErr heightErr];
    set(fig_1, 'Position', posScrn, 'Color', [1,1,1]);
           
    a(1) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.11 0.155 0.87 0.82]);

    %-- Compression
    fig_2 = figure(2);     clf(fig_2);
    
    posScrn = [0 sizeScrn(2)-heightErr-heightOtr-187 ...
               widthOtr heightOtr];
    set(fig_2, 'Position', posScrn, 'Color', [1,1,1]);
    a(2) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.11 0.20 0.87 0.81]);

    %-- Ensemble sizes
    fig_3 = figure(3);     clf(fig_3);
    posScrn = [0 sizeScrn(2)-heightErr-2*heightOtr-273 ...
               widthOtr heightOtr];
    set(fig_3, 'Position', posScrn, 'Color', [1,1,1]);
    a(3) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.10 0.20 0.88 0.81]);
end

%-- Matrixes initialization
nameFile=sprintf('../database/%s/%s%dBlocksExt/D1^t1.db',nameDb,nameSc,nBlocks);
dbase = util_readHeader(nameFile);

%-- Confidence interval
confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
pValue       = norminv(confInterval) / sqrt(nReplications);

%------------------------------------------------------------------------------%
%---------------------- Calcul des points du graphiques -----------------------%
%------------------------------------------------------------------------------%
if compute
    graphErr = zeros(3*nTest, nBlocks);
    graphCpn = zeros(3*nTest, nBlocks);
    graphAcn = zeros(3*nTest, nBlocks);
    cpnTest  = zeros(nBlocks, nReplications);
    graphEns = zeros(3*nTest, nBlocks);

    cpnAv = zeros(nBlocks, nReplications);
    cpn   = zeros(nBlocks, nReplications);

    textLegend = '';

    tTest = 1;
    while tTest <= nTest
        
%         if tTest == 1,   tTest = 2;   end
%         if tTest == 2,   tTest = 3;   end
%         if tTest == 3,   tTest = 4;   end
%         if tTest == 4,   tTest = 5;   end

        %-- Name parameters & legend in function of the current curve
        name = util_getName(nameTest, tTest);

		if strcmp(name.ln, 'Bth'),   nBlocksReal = 1;  nReplications = 10;
        else                         nBlocksReal = nBlocks;
        end
        
        %-- Result file loaded
        nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.result',...
                 nameDb, nameSc, name.ln, name.hp, name.al, nBlocks, name.width)
        result = util_readResults(nameFile, nBlocksReal, nReplications);

        %-- Legend
        textLegend(tTest,1:length(name.legend)) = name.legend;

        %------------ Graph points - for all blocs & replications -------------%
        %-- f1(x) - Classification rate
        graphErr(tTest*3-2, :) = 1:nBlocks;
        graphErr(tTest*3-1, :) = mean(result.clsRate')*100;
        graphErr(tTest*3,   :) = std (result.clsRate')*100 * pValue;

        %-- f2(x) - Compression
        for t = 1:nBlocksReal
            for r = 1:nReplications
				[r t];
                sizeEns  = result.sizeEns(t,r);
                sizeF2   = sum(sum(result.sizeClasses(:,1:sizeEns,t,r)));
                sizeF2av = mean(sum(result.sizeClasses(:,1:sizeEns,t,r)));
                if ~sizeF2,     sizeF2   = result.nClasses;   end
                if ~sizeF2av,   sizeF2av = result.nClasses;   end
                cpnAv(t,r) = result.nPatternsLearned(t,r) / sizeF2av;
                cpnAv(t,r) = sizeF2av;
                cpn(t,r)   = result.nPatternsLearned(t,r) / sizeF2;
                cpn(t,r)   = sizeF2;
            end
        end
        graphCpn(tTest*3-2, :) = 1:nBlocks;
        graphCpn(tTest*3-1, :) = mean(cpn');
        graphCpn(tTest*3,   :) = std(cpn') * pValue;

        graphAcn(tTest*3-2, :) = 1:nBlocks;
        graphAcn(tTest*3-1, :) = mean(cpnAv');
        graphAcn(tTest*3,   :) = std(cpnAv') * pValue;

        %-- f3(x) - Ensemble size
        graphEns(tTest*3-2, :) = 1:nBlocks;
        graphEns(tTest*3-1, :) = mean(result.sizeEns');
        graphEns(tTest*3,   :) = std (result.sizeEns') * pValue;
        
        tTest = tTest +1;
    end
end

%------------------------------------------------------------------------------%
%------------------------------ Saving & loading ------------------------------%
if saveResult
    nameFile = sprintf('../savedStuff/%s_%s%s.mat', nameDb, nameSc, nameTest);
    save(nameFile, 'graphErr', 'graphCpn', 'graphAcn', 'graphEns', '-mat');
end

if loadResult
    nameFile = sprintf('../savedStuff/%s_%s%s.mat', nameDb, nameSc, nameTest);
    load(nameFile, 'graphErr', 'graphCpn', 'graphAcn', 'graphEns', '-mat');
end

%------------------------------------------------------------------------------%
%------------------------------ Hypothesis tests ------------------------------%
if hypTest
    y       = zeros(nTest, nBlocks);
    ic      = zeros(nTest, nBlocks);
    tested  = zeros(nTest, nBlocks);
    resHyp  = zeros(nTest, nBlocks);
    df      = zeros(1,     nBlocks);
    pValues = zeros(1,     nBlocks);

    for tTest = 1:nTest
        y(tTest,:)  = graphErr(tTest*3-1, :);
        ic(tTest,:) = graphErr(tTest*3, :) / pValue;
    end

    for tTest = 1:nTest-1
        mu1 = y(tTest,:);
        mu2 = y(tTest+1,:);

        std1 = ic(tTest,:); 
        std2 = ic(tTest+1,:); 

        tested(tTest,:) = abs(mu2-mu1) ./ ...
                                 sqrt( (std1.^2+std2.^2) / nReplications );

        df(1,:) = (nReplications-1) * (std1.^2+std2.^2).^2 ./ (std1.^4+std2.^4);
        pValues = tinv( confInterval, df );

        resHyp(tTest,:) = heaviside(tested(tTest,:)-pValues);
    end

    resHyp
end

%------------------------------------------------------------------------------%
%----------------------------------- Graphs -----------------------------------%
%-- Legend settings
textLegend = '';
for tTest = 1:nTest
    %-- Name parameters & legend in function of the current curve
    name = util_getName(nameTest, tTest);

    textLegend(tTest,1:length(name.legend)) = name.legend;
end

textLegend( nTest+1, 1:size(legLocation,    2) ) = legLocation;
textLegend( nTest+2, 1:size(legOrientation, 2) ) = legOrientation;
%graphErr(1:3,:) = [];
%graphCpn(1:3,:) = [];
%graphEns(1:3,:) = [];
if graph
    %-- Classification rate
    set(0,'CurrentFigure',1);
    gErr     = util_graphBvsi(graphErr, nTest, classRate, textLegend);
    nameFile = sprintf('export_fig ../figures/%s_%s%sErr -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);

    %-- Compression
    set(0,'CurrentFigure',2);
    gCpn = util_graphBvsi(graphCpn, nTest, compression, textLegend);
    nameFile = sprintf('export_fig ../figures/%s_%s%sCpn -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);

    %-- Ensemble sizes
    set(0,'CurrentFigure',3);
    gEns = util_graphBvsi(graphEns, nTest, sizeEnsemble, textLegend);
    nameFile = sprintf('export_fig ../figures/%s_%s%sEns -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);
end
 
if show
    if strcmp(nameTest,'total')
        fprintf('\n   MA+d      LB+d      MOPSO     ');
        fprintf('gbest     knn       Batch\n');
    else
        fprintf('\n   LB+d      all       gbest\n');
    end
    [100-graphErr([2:3:nTest*3-1], nBlocks) , graphErr([3:3:nTest*3], nBlocks)]'
    [graphEns([2:3:nTest*3-1], nBlocks) , graphEns([3:3:nTest*3], nBlocks)]'
    [graphAcn([2:3:nTest*3-1], nBlocks) , graphAcn([3:3:nTest*3], nBlocks)]'
    [graphCpn([2:3:nTest*3-1], nBlocks) , graphCpn([3:3:nTest*3], nBlocks)]'
end
fprintf('/*----------------- End of results -----------------*/\n');
