%------------------------------------------------------------------------------%
%-- Results: batch -vs- inremental for different settings
%-- Works with util_getName
%------------------------------------------------------------------------------%
fprintf('\n/*--------------- Results - settings ---------------*/\n');

%-- DEFINE\
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
%-- File name
nameDb         = 'cnrc64';%  {cnrc, mobo}
nameSc         = 'upd';%   {add, upd}
nBlocks        =  12;
nameAl         = 'dnpso';
nameTest       = 'multi';%   {'Ens','Vid','tst'}
legLocation    = 'EastOutside';
legOrientation = 'vertical';
wdVec          = [1600];
wdVec          = [50 100 200 400 800 1600];

%-- Other parameters
typeTest = 1;
nReps    = 10;

%-- Graphs
compute    = 1;
loadResult = 0;
graph      = 1;

%-- Graph width & ratio
typeDoc  = elsevier; %  elsevier;
ratioErr = 1.2; %  {1.75, 2.5};
ratioOtr = 1.2;

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
errWidth = errWidth/2.1;
ratioErr = 1.5;

%-- Confidence interval
confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
pValue       = tinv(confInterval, nReps-1) / sqrt(nReps);

%-- Name
name = util_getName(nameTest, typeTest);

%------------------------------------------------------------------------------%
%------------------------- Reading the Test data base -------------------------%
nameFile  = sprintf('../database/%s/test.db', nameDb);
dbase     = util_readDb(nameFile);
nPatterns = dbase.nPatterns;
nFeatures = dbase.nFeatures;
A(1:nPatterns,1:nFeatures)              =     dbase.data(:,:);
A(1:nPatterns,nFeatures+1: 2*nFeatures) = 1 - dbase.data(:,:);
tags      = dbase.tags;

%------------------------------------------------------------------------------%
%---------------------------------- Computing ---------------------------------%
if compute 
    nBlks = 1;
    t     = 1;
    nTest = size(wdVec, 2);

    %-- Matrixes initialization
    Qstat    = zeros(2, nTest);
    Rho      = zeros(2, nTest);
%     divFam   = zeros(2, nTest);
    divSwarm = zeros(2, nTest);
    
    tempNormWa = zeros(1,nPatterns);
    normWalpha = zeros(nPatterns, nPatterns);
    Wj         = zeros(nPatterns, nFeatures);
    normAW     = zeros(nPatterns, nFeatures);
    T          = zeros(nPatterns, nFeatures/2);
    T_j2       = zeros(nPatterns, nFeatures/2);
    T_jStar    = zeros(nPatterns, 1);
    T_jStar2   = zeros(nPatterns, 1);
    
    for tTest = 1:nTest
        widthMc = wdVec(tTest);
        if tTest == nTest,   currentTest = 9;
        else,                currentTest = typeTest;
        end
        name    = util_getName(nameTest, currentTest);

        fprintf('/*-- Width: %d, hp: %s\n', wdVec(tTest), name.hp);
        
        %------------------------------- Loading ------------------------------%
        fprintf('     Load - pso');
        %-- Archive
        nameFile = sprintf('../savedStuff/%s_updBth%s%s_%dBlocks_wd%d.pso', ...
                           nameDb, 'Hdnc', nameAl, nBlocks, widthMc);
        pso = util_readPsoMulti(nameFile, nBlks, nReps);
        nDimensions = pso.nDimensions;
        
        %-- Archive
        fprintf('   archive');
        nameFile = sprintf('../savedStuff/%s_updBth%s%s_%dBlocks_wd%d.arc', ...
                           nameDb, 'Hdnc', nameAl, nBlocks, widthMc);
        archive  = util_readArchive(nameFile, nBlks, nReps, pso.nIterations);
        szArchive = archive.size;

        %-- Result
        fprintf('   result\n');
        nameFile = sprintf('../savedStuff/%s_upd%s%s%s_%dBlocks_wd%d.result',...
                           nameDb, 'Bth', name.hp, nameAl, nBlocks, widthMc);
        result   = util_readResults(nameFile, nBlks, nReps);

        %-- Indicator initialization
        maxSizeEns = max(result.sizeEns);
        Q_av       = zeros(nReps,1);
        rho_av     = zeros(nReps,1);

        %----------------------------------------------------------------------%
        %----------------------- FAM diversity indicator ----------------------%
%         fprintf('     Compute\n       FAM diversity');
%         nameFam = sprintf('../savedStuff/%s_updBth%s%s_%dBlocks_wd%d.artmap',...
%                            nameDb, name.hp, nameAl, nBlocks, widthMc);
%         ambiguityEns = zeros(1,nReps);
%         for r = 1:nReps,   fprintf(' %d', r);
%             
%             sizeEns = result.sizeEns(r);
%             for e = 1:sizeEns
% % e
%                 ambiguity = zeros( 1, sizeEns               );
%                 dist      = zeros( 1, sizeEns*(sizeEns-1)/2 );
% 
%                 fam = util_readFam(nameFam, nBlks, nReps, result.sizeEns,r,t,e);
%                 sizeNn = fam.sizeNn;
% 
%                 %-- 1 / (|Wij| + alpha) for all the patterns
%                 tempNormWa(1:sizeNn)  = transpose(fam.normWalpha(1:sizeNn));
%                 normWalpha(:,1:sizeNn)= ...
%                                repmat( tempNormWa(1:sizeNn), nPatterns, 1 );
% 
%                 %------ Compute T(j) = |min(A,Wij)| / (|Wij| + alpha) -----%
%                 %-- |min(A,Wij)| for all patterns and all nodes
%                 normAW = zeros(nPatterns, sizeNn);
%                 for j = 1:sizeNn
%                     %-- Temp - W
%                     Wj = repmat(fam.W(j,:), nPatterns, 1);
% 
%                     %-- |min(A,Wij)| for all the patterns
%                     normAW(:,j) = sum( min(A(:,:),Wj), 2 );
%                 end
% 
%                 %-- T(j) for all patterns
%                 T(:,1:sizeNn) = ...
%                     normAW(:,1:sizeNn).*normWalpha(:,1:sizeNn);
% 
%                 %-- T_jStar(a)
%                 [T_jStar winners] = max(T(:,1:sizeNn),[],2);
% 
%                 %-- max(T_j(a), k~=kStar)
%                 [tmp kj_one] = max( transpose( fam.Wab ) );
%                 kStar_one    = result.predictions(:,e,1,r)+1;
% 
%                 kj    = repmat( kj_one, nPatterns, 1);
%                 kStar = repmat( kStar_one, 1, sizeNn);
% 
%                 filterPred = not( eq(kStar,kj) );                
% 
%                 T_j2(:,1:sizeNn)  = T(:,1:sizeNn) .* filterPred;
%                 [T_jStar2 second] = max(T_j2(:,1:sizeNn),[],2);
% 
%                 ambiguity(e) = sum( T_jStar - T_jStar2 );
%             end
%             
%             %--------------------- Values for the ensemble --------------------%
%             for i = 1:sizeEns-1
%                 for j = i+1:sizeEns
%                     ambiguityEns(r) = abs(ambiguity(i) - ambiguity(j)) + ...
%                                       ambiguityEns(r);
%                 end
%             end
%         end
%         
%         %-- Final answer!
%         divFam(1, tTest) = mean(ambiguityEns);
%         divFam(2, tTest) = std (ambiguityEns) * pValue;

        %----------------------------------------------------------------------%
        %----------------- Q-stat and rho diversity indicators ----------------%
        fprintf('\n       Standard measures');
        for r = 1:nReps,   fprintf(' %d', r);
            sizeEnsemble = result.sizeEns(r);
            Q            = zeros(sizeEnsemble, sizeEnsemble);
            rho          = zeros(sizeEnsemble, sizeEnsemble);

            for i = 1:sizeEnsemble-1
                for j = i+1:sizeEnsemble

                    %-- Actual tested members
                    predi = result.predictions(:,i,t,r);
                    predj = result.predictions(:,j,t,r);
                    truec = result.trueClasses(:,t,r);

                    %-- Good classification by i & j
                    goodi = eq( predi, truec );
                    goodj = eq( predj, truec );

                    %-- N11, N00 - Similar good & bad classification,
                    n11 = and( goodi,      goodj      );    N11 = sum(n11);
                    n00 = and( not(goodi), not(goodj) );    N00 = sum(n00);

                    %-- N10 - Different - i,j = good, j,i = bad
                    n10 = and( goodi,      not(goodj) );    N10 = sum(n10);
                    n01 = and( not(goodi), goodj      );    N01 = sum(n01);

                    %-- Qstat & rho
                    Q  (i,j) = (N11*N00 - N10*N01) / (N11*N00 + N10*N01);
                    rho(i,j) = (N11*N00 - N10*N01) / ...
                              (sqrt( (N11+N10)*(N01+N00)*(N11+N01)*(N10+N00) ));
                end
            end

            %-- Average Qstat and rho of the ensemble
            for i = 1:sizeEnsemble-1
                for j = i+1:sizeEnsemble
                    Q_av(r)   = Q_av(r)   + Q  (i,j);
                    rho_av(r) = rho_av(r) + rho(i,j);
                end
            end

            Q_av(r)   = Q_av(r)   *2 / (sizeEnsemble *(sizeEnsemble-1));
            rho_av(r) = rho_av(r) *2 / (sizeEnsemble *(sizeEnsemble-1));
        end
        
        %-- Final answer!
        Qstat(1, tTest) = mean(Q_av)
        Qstat(2, tTest) = std (Q_av) * pValue;
        
        Rho(1, tTest) = mean(rho_av);
        Rho(2, tTest) = std (rho_av) * pValue;
        
        %----------------------------------------------------------------------%
        %------------- Diversity in the optimization environment --------------%
%         fprintf('\n       Swarm            ');
%         av_d = zeros(1, nReps);
%         for r = 1:nReps,   fprintf(' %d', r);
%             
%             %-- Particles forming the ensemble
%             nIt          = pso.nIterations(r);
%             data         = zeros(sizeEnsemble, nDimensions);
%             sizeEnsemble = result.sizeEns(r);
%             distances    = zeros(sizeEnsemble, sizeEnsemble);
%             d_e          = zeros(1, sizeEnsemble *(sizeEnsemble-1));
%             
%             for e = 1: sizeEnsemble
%                 data(e,:) = archive.s(result.members(e,t,r), :, nIt, t, r);
%             end
% 
%             %-- Computing the distances
%             for i = 1:sizeEnsemble-1
%                 for j = i+1:sizeEnsemble
%                     dst = 0;
%                     for d = 1:nDimensions
%                         s1=data(i,d);   s2=data(j,d);   dst=dst + (s2-s1)^2;
%                     end
%                     distances(i,j,r) = sqrt(dst);
%                 end
%             end
%             
%             %-- Average value of the distances
%             k = 0;
%             for i = 1:sizeEnsemble-1
%                 for j = i+1:sizeEnsemble
%                     k       = k+1;
%                     d_e(k) = distances(i,j);
%                 end
%             end
%         end
%         
%         %-- Final answer!
%         divSwarm(1, tTest) = mean(d_e);
%         divSwarm(2, tTest) = std (d_e) * pValue;
        
        fprintf('\n       Done\n');
    end
    
    %-- Saving the data
%     nameData = sprintf('../savedStuff/%s_%sBth_div.mat', nameDb, nameSc);  
%     save(nameData, 'Qstat', 'Rho', 'divFam', 'divSwarm', '-mat');
end

%------------------------------------------------------------------------------%
%------------------------------ Loading the data ------------------------------%
if loadResult
    nameData = sprintf('../savedStuff/%s_%sBth_div.mat', nameDb, nameSc);  
    load(nameData);
end

%------------------------------------------------------------------------------%
%---------------------------------- Graphics ----------------------------------%
if graph
    
    sizeLine   = 1;      sizeLineCI = 0.5;
    sizeMrkr   = 4;      sizeMrkrCI = 2;
    sizeFont   = 8;      sizeFontL  = 10;

    %-- Initialization
    close all;
    widthErr  = res*errWidth;     heightErr = res*errWidth/ratioErr;

    %-- Diversity FAM
    fig_1 = figure(1);     clf(fig_1);
    posScrn = [0 sizeScrn(2)-heightErr-104 widthErr heightErr];
    set(fig_1, 'Position', posScrn, 'Color', [1,1,1]);
           
    a(1) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.19 0.18 0.79 0.82]);

    hold on            
    %-------------------------------- Reference -------------------------------%
    graphColor = ones(1,3)*0.2;
    
    xVec = [ 0, 800, 800, 0, 0];
    yVec = [ Qstat(1,nTest)+Qstat(2,nTest), Qstat(1,nTest)+Qstat(2,nTest), ...
             Qstat(1,nTest)-Qstat(2,nTest), Qstat(1,nTest)-Qstat(2,nTest), ...
             Qstat(1,nTest)+Qstat(2,nTest) ];
         
    g = fill(xVec,yVec, graphColor+0.2, 'EdgeColor', graphColor, ...
             'LineWidth', 0.5);
         
    %-- Curve (only)                
    plot([0 800], [Qstat(1,nTest) Qstat(1,nTest)], '-', 'Color', graphColor, ...
         'LineWidth', sizeLine);

    graphColor = ones(1,3)*0.8;
    xVec = [ 0, 800, 800, 0, 0];
    yVec = [ Rho(1,nTest)+Rho(2,nTest), Rho(1,nTest)+Rho(2,nTest), ...
             Rho(1,nTest)-Rho(2,nTest), Rho(1,nTest)-Rho(2,nTest), ...
             Rho(1,nTest)+Rho(2,nTest) ];

    g = fill(xVec,yVec, graphColor+0.2, 'EdgeColor', graphColor, ...
             'LineWidth', 0.5);
         
    %-- Curve (only)                
    plot([0 800], [Rho(1,nTest) Rho(1,nTest)], '-', 'Color', graphColor, ...
         'LineWidth', sizeLine);
     
    %--------------------------------- Curves ---------------------------------%
    for i = 1:2
        
        switch i
            case 1,  data = Qstat;      graphType = 'o-';
            case 2,  data = Rho;        graphType = 's-';
            case 3,  data = divSwarm;   graphType = 'v-';
        end

        clr = 0.2 + 0.5 * (i-1)/(3-1) * ones(1,3);

        %-- Curves (incremental)
        g(i) = plot(wdVec(1:nTest-1), data(1,1:nTest-1), graphType, ...
                     'MarkerSize', sizeMrkr, 'Color', clr, ...
                     'MarkerFaceColor', clr, 'LineWidth', sizeLine);

        %-- Custumized error bar
        for tTest = 1:nTest-1
            x        = wdVec(tTest);
            negBound = data(1, tTest) - data(2, tTest);
            posBound = data(1, tTest) + data(2, tTest);

            z = plot([x,x], [negBound, posBound], graphType,     ...
                   'LineWidth', sizeLineCI, 'Color', clr, ...
                   'MarkerFaceColor', clr, 'MarkerSize', sizeMrkrCI);
        end
    end
    
    
    hold off
    set(gca,'FontSize', sizeFont, 'XLim', [0 800]);
    
    
    %-- Legend only for the classification rate
%     textLeg = ['$\mathit{EoMA}_t^\mathrm{+d}$'; '$\mathit{EoMA}_t^\mathrm{+d}$'];
    textLeg = ['$Q_{\textnormal{stat}}$'; '$\rho$                 '];
    leg = legend(g(1:2), textLeg);
    set(leg, 'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
             'Orientation', 'Vertical', 'Interpreter','Latex', ...
             'Location', 'NorthEast');
    
    %-- Labels definition
    xlabel('Width of each region', ...
           'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
           'FontSize', sizeFont, 'Interpreter','Latex');
    ylabel('Classifer correlation', 'FontName', 'Times New Roman',...
           'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

    nameFile = sprintf('export_fig ../figures/%s_div -pdf', nameDb);
    eval(nameFile);
end

fprintf('/*----------------- End of results -----------------*/\n');
