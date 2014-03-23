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
nameDb      = 'mobof32';
nameSc      = 'upd';%   {add, upd}
nBlocks     =  26;
legLocation = 'EastOutside';
nameTest    = 'dnpso';%   {'vid'}

compute    = 1;
saveResult = 1;
loadResult = 0;
graph      = 0;

%-- Other parameters
time          = 26;
nTest         = 3;
nReplications = 20;

%-- Graph width & ratio
typeDoc = elsevier; % {ieee, elsevier};
ratio   = 1.75; % 2.5; {1.75, 2.5}

%---------------------------- Graphics parameters -----------------------------%
%-- Graphs width
switch typeDoc
    case elsevier
        errWidth = 16.64/2;
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

%-- Read Info MoBo
mobo         = util_readInfoMobo();
nPatternsMax = max(max(max(mobo.info)));

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

        %-- Result file loaded
        nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.result',...
                nameDb, nameSc, name.ln, name.hp, name.al, nBlocks, name.width)

        result = util_readResults(nameFile, nBlocksReal, nReplications);

        %------------ Graph points - for all blocs & replications -------------%
        nPatTr = zeros(1,nReplications);

        %-- Graph points themselves
		t = 1;
        addrT = 1;
        while t <= nBlocks

            equal = zeros(mobo.nSequences, nPatternsMax);
            cRate = zeros(nReplications,   nPatternsMax);
            nSeqTot = 0;
            for r=1:nReplications 

                nPatternsT          = result.nPatternsTest(t,r);
                truec(1:nPatternsT) = result.trueClasses(1:nPatternsT, t, r) +1;

                %-- Classes activity and sizes
                activeCls                      = zeros(1,nClasses);
                activeCls(truec(1:nPatternsT)) = 1;
                nActive                        = sum(activeCls);

                %-- Number of sequences
                nSequences = 0;
                for s = 1:mobo.nSequences
                    addrS    = mobo.addr(s)+1;
                    cClasses = dbase.tags(addrS)+1;
                    if activeCls(cClasses)
                        nSequences = nSequences+1;
                    end
                end

                addrT = 1;
                for s = 1:mobo.nSequences
                    addrS = mobo.addr(s)+1;
                    cClass = dbase.tags(addrS)+1;
                    
                    if activeCls(cClass && mobo.info(s))
                        
                        nSeqTot   = nSeqTot +1;
                        nPatterns = mobo.info(s);
                        addrE     = addrT + nPatterns -1;
                        
                        predi = result.predisemble(addrT:addrE, t, r) +1;
                        
                        %-------- Build histogram for the current class -------%
                        predictions = zeros(1,nPatternsMax);
                        
                        %-- Video sequences contributing
                        vHist = zeros(1, nClasses);
                        for p = 1:nPatterns
                            addr            = addrS +p -1;
                            vHist(predi(p)) = vHist(predi(p)) +1;
                            [temp pred]     = max(vHist);
                            predictions(p)  = pred;
                        end
                        
                        %-- Use the last prediction made to extend
                        for p = nPatterns+1:nPatternsMax
                            predictions(p)     = pred;
                        end

                        %-- "truec" vector
                        truec = ones(1,nPatternsMax)*cClass;

                        %---------------- Finds the classification rate ---------------%
                        equal(nSeqTot,:) = eq(predictions, truec);
                        
                        addrT = addrT + nPatterns;
                    end
                end
                cRate(r,:) = mean(equal(1:nSeqTot,:))*100;
            end  %-- for : nRep

            
            %-- f1(x) - Classification rate
            if tTest == 5,   t = nBlocks;   end
                
            graphErr(tTest*2-1, :, t) = 100 - mean(cRate  );
            graphErr(tTest*2  , :, t) =       std (cRate,0) *pValue;
            
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
    g = util_graphVideoTime(graphErr, time);
    
    nameFile = sprintf('export_fig ../figures/%s_%s%sVideo_t%d -pdf', ...
                       nameDb, nameSc, nameTest, time);
    eval(nameFile);
end
fprintf('/*----------------- End of results -----------------*/\n');
