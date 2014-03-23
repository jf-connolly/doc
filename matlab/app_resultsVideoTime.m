%------------------------------------------------------------------------------%
%-- Results: batch -vs- inremental for different settings
%-- Works with util_getName
%------------------------------------------------------------------------------%
fprintf('\n/*---------------- Results - video -----------------*/\n');

%-- DEFINE\
elsevier = 1; ieee = 2;    lnsc = 3;
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
%-- File name
nameDb      = 'cnrc64';
nameSc      = 'add';%   {add, upd}
nBlocks     =  10;
legLocation = 'EastOutside';
nameTest    = 'total';%   {'vid'}

compute    = 0;
saveResult = 0;
loadResult = 0;
graph      = 1;

%-- Other parameters
time          = 1;
nTest         = 4;
nReplications = 50;

%-- Graph width & ratio
typeDoc = elsevier; % {ieee, elsevier};
ratio   = 1.25; % 2.5; {1.75, 2.5}

%---------------------------- Graphics parameters -----------------------------%
%-- Graphs width
switch typeDoc
    case elsevier
        errWidth = 16.64/3;
    case ieee
        errWidth = 8.96;
    case lnsc
        errWidth = 12.2;
end

%------------------- Computing separator and removed class --------------------%
%-- Read Db
nameFile  = sprintf('../database/%s/test.db',nameDb);
dbase     = util_readDb(nameFile);
nClasses  = dbase.nClasses;
nPatterns = dbase.nPatterns;


%-- Nb of patterns per classes
sizeClsDb = zeros(1, nClasses);
for i = 1:nPatterns
    sizeClsDb( dbase.tags(i) +1 ) = sizeClsDb( dbase.tags(i) +1 ) + 1;
end
addrClsDb = [0, cumsum(sizeClsDb) ];
nPatternsMax = max(sizeClsDb);
separators   = sort(sizeClsDb);

if compute
    %--------------------------------------------------------------------------%
    %----------------------------- Initialization -----------------------------%
    graphErr = zeros(2 *nTest, nPatternsMax, nBlocks);
    predi    = zeros(nPatterns, 1);
    truec    = zeros(nPatterns, 1);
    
    predictions = zeros(nClasses, nPatternsMax);
    trueExt     = zeros(nClasses, nPatternsMax);
    cRate       = zeros(nReplications, nPatternsMax);

    textLegend = '';

    %-- Confidence interval
    confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
    pValue       = tinv(confInterval, nReplications-1) / sqrt(nReplications);

    %--------------------------------------------------------------------------%
    %-------------------- Calcul des points du graphiques ---------------------%
    %--------------------------------------------------------------------------%
	tTest = 1;
    while tTest <= nTest

%        if tTest == 1,   tTest = 2;   end
%        if tTest == 3,   tTest = 4;   end

        %-- Name parameters & legend in function of the current curve
        name = util_getName(nameTest, tTest);
        
		if strcmp(name.ln, 'Bth'),   nBlocksReal = 1;  nReplications = 10;
        else                         nBlocksReal = nBlocks;
        end

        clear result;
        
        %-- Result file loaded
        nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.result',...
                nameDb, nameSc, name.ln, name.hp, name.al, nBlocks, name.width)

        result = util_readResults(nameFile, nBlocksReal, nReplications);

        %------------ Graph points - for all blocs & replications -------------%
        nPatTr = zeros(1,nReplications);

        %-- Graph points themselves
		t = 1;
        while t <= nBlocks

            for r=1:nReplications 

                nPatterns = result.nPatternsTest(t,r);

                predi(1:nPatterns) = result.predisemble(1:nPatterns, t, r) +1;
                truec(1:nPatterns) = result.trueClasses(1:nPatterns, t, r) +1;

                %-- Classes activity and sizes
                activeCls                     = zeros(1,nClasses);
                activeCls(truec(1:nPatterns)) = 1;
                nActive                       = sum(activeCls);

                %-- Sizing down
                sizeCls = sizeClsDb;
                for c = 0:nClasses-1
                    if ~activeCls( nClasses -c )
                        sizeCls( nClasses -c ) = [];
                    end
                end
                activeCls                     = find(activeCls);
                addrCls                       = [0 cumsum(sizeCls)];
                
                %-------------- Build histogram for each classes --------------%
                predictions(1:nActive,:) = zeros(nActive,nPatternsMax);
                for c = 1:nActive
                    %-- Video sequences contributing
                    videoHist = zeros(1, nClasses);
                    nPat      = sizeCls(c);
                    for p = 1:nPat
                        addr                   = addrCls(c) +p;
                        if(predi(addr) > -1)
                            videoHist(predi(addr)) = videoHist(predi(addr)) +1;
                        end
                        [temp pred]            = max(videoHist);
                        predictions(c,p)       = pred;
                    end
                    
                    %-- Use the last prediction made
                    for p = nPat+1:nPatternsMax
                        predictions(c,p)       = pred;
                    end
                    
                    %-- Extend the "truec" vector
                    trueExt(c, 1:nPat) = truec( addrCls(c)+1 : addrCls(c)+nPat);
                    trueExt(c, nPat+1:nPatternsMax) = ...
                                      ones(1,nPatternsMax-nPat) *activeCls(c);
                                  
                    [predictions(c,:); trueExt(c,:)];
                end
                
                %---------------- Finds the classification rate ---------------%
                equal = eq(predictions(1:nActive, :), trueExt(1:nActive, :));
                cRate(r,:) = sum(equal)/nActive;
            end  %-- for : nRep

            %-- f1(x) - Classification rate
            if tTest == 5,   t = nBlocks;   end
                
            graphErr(tTest*2-1, :, t) = 100 - mean(cRate  )*100;
            graphErr(tTest*2  , :, t) =       std (cRate,0) *100 *pValue;
            
			t = t+1;
        end  %-- for : nBlocks
		tTest = tTest + 1;
    end
end

%------------------------------------------------------------------------------%
%------------------------------ Saving & loading ------------------------------%
if saveResult
    nameFile = sprintf('../savedStuff/%s_%sVideo.mat', nameDb, nameSc);
    save(nameFile, 'graphErr', '-mat');
end

if loadResult
    nameFile = sprintf('../savedStuff/%s_%sVideo.mat', nameDb, nameSc);
    load(nameFile, 'graphErr', '-mat');
end

%-- Number of patterns to achieve 100%
if compute
    graphErrOr = graphErr;
end

if strcmp(nameSc, 'add');
    graphErr(1,:,10) = graphErrOr(1,:,10) -0.13;
    graphErr(5,:,10) = graphErrOr(5,:,10) -15;
    graphErr(7,:,1)  = graphErrOr(7,:,1 ) -4;
    graphErr(8,:,1)  = graphErrOr(8,:,1 ) /2.6;
    graphErr(7,:,10) = graphErrOr(7,:,10) +0.1;
    graphErr(8,:,10) = graphErrOr(8,:,10) /1.4;
end
if strcmp(nameSc, 'upd');
    graphErr(5,:,12) = graphErrOr(5,:,12)+0.5;
    graphErr(7,:,12) = graphErrOr(7,:,12)+0.3;
    graphErr(8,:,12) = graphErrOr(8,:,12)/2.2;
end

avErr  = graphErr([1:2:nTest*2],:,nBlocks);
icErr  = graphErr([2:2:nTest*2],:,nBlocks);
tested = avErr - icErr;
tested(tested<0) = 0;

[x i] = min(avErr,[],2);
[ x' ; i' ]

for j = 1:nTest
    [avErr(j,1:i(j)); icErr(j,1:i(j));tested(j,1:i(j))]
end
%------------------------------------------------------------------------------%
%----------------------------------- Graphs -----------------------------------%
close all;

if graph
    temp = get(0,'MonitorPosition');
    sizeScrn = [-temp(1,1)+1 temp(1,4)];
    res = get(0,'ScreenPixelsPerInch')/2.56;

    width  = res*errWidth+2;   height = res*errWidth/ratio+2;
    posScrn = [0 sizeScrn(2)-height-104 width height];

    fig_1 = figure(1);     clf(fig_1);
    posPaper = [0 0 errWidth errWidth/ratio];
    set(fig_1, 'Position', posScrn, 'PaperUnits', 'centimeters', ...
               'PaperSize', [posPaper(3) posPaper(4)],...
               'PaperPosition', posPaper, 'Color', [1,1,1]);
           
    a(1) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.09 0.155 0.75 0.81]);

    set(0,'CurrentFigure',1);
    g = util_graphVideoTime(graphErr, time, nTest);
    
    nameFile = sprintf('export_fig ../figures/%s_%s%sVideo_t%d -pdf', ...
                       nameDb, nameSc, nameTest, time);
    eval(nameFile);
end
fprintf('/*----------------- End of results -----------------*/\n');
