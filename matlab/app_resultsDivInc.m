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
nameAl         = 'multi';
nameTest       = 'vid';%   {'Ens','Vid','tst'}
legLocation    = 'EastOutside';
legOrientation = 'vertical';

%-- Other parameters
nTest         = 2;
nReplications = 50;

%-- Graphs
compute    = 0;
loadResult = 1;
graph      = 1;

%-- Graph width & ratio
typeDoc  = elsevier; %  elsevier;
ratioErr = 1.2; %  {1.75, 2.5};
ratioOtr = 1.2;

%---------------------------- Graphics parameters -----------------------------%
%-- Graphs width
switch typeDoc
    case elsevier
        errWidth = 16.64/3.1;
    case ieee
        errWidth = 8.96;
    case lnsc
        errWidth = 12.2;
end

%-- Confidence interval
confInterval = 0.95;  %-- 1 - alpha/2, or alpha = 10%
pValue       = tinv(confInterval, nReplications-1) / sqrt(nReplications);

%------------------------------------------------------------------------------%
%---------------------- Calcul des points du graphiques -----------------------%
%------------------------------------------------------------------------------%
if compute
    %-- Test data base
    nameFile  = sprintf('../database/%s/test.db', nameDb);
    dbase     = util_readDb(nameFile);
    nPatterns = dbase.nPatterns;
    nFeatures = dbase.nFeatures;
    A(1:nPatterns,1:nFeatures)              =     dbase.data(:,:);
    A(1:nPatterns,nFeatures+1: 2*nFeatures) = 1 - dbase.data(:,:);
    tags      = dbase.tags;

    %-- Read swarm
    nameFile = sprintf('../savedStuff/position/%s_%s%s%s%s_%dBlocks.pso', ...
                       nameDb, nameSc, 'Inc', 'Hdnc', nameAl, nBlocks);
    pso = util_readPsoMulti(nameFile, nBlocks, nReplications);
    nDimensions = pso.nDimensions;

	%-- Matrixes initialization
	graphDiv   = zeros(9, nBlocks);
	graphDivS  = zeros(9, nBlocks);
	tempNormWa = zeros(1,nPatterns);
	normWalpha = zeros(nPatterns, nPatterns);
	Wj         = zeros(nPatterns, nFeatures);
	normAW     = zeros(nPatterns, nFeatures);
	T          = zeros(nPatterns, nFeatures/2);
	T_j2       = zeros(nPatterns, nFeatures/2);
	T_jStar    = zeros(nPatterns, 1);
	T_jStar2   = zeros(nPatterns, 1);

    for typeTest = 1:nTest
        %-- Name parameters & legend in function of the current curve
        name = util_getName(nameTest, typeTest);

        %-- Result file loaded
        nameFile = sprintf('../savedStuff/position/%s_%s%s%s%s_%dBlocks.result', ...
                           nameDb, nameSc, name.ln, name.hp, nameAl, nBlocks);

        result = util_readResults(nameFile, nBlocks, nReplications);

        %-- Reading the FAMs
        nameFam = sprintf('../savedStuff/position/%s_%s%s%s%s_%dBlocks.artmap', ...
                           nameDb, nameSc, name.ln, name.hp, nameAl, nBlocks);

        %------------ Graph points - for all blocks & replications ------------%
        for t = 1:nBlocks

            ambiguityEns = zeros(1,nReplications);
            for r = 1:nReplications 
                fprintf('/-- %d, time: %d,replication: %d\n', typeTest, t,r);
			    %------------ Average margin of all the classifiers -----------%
                sizeEns = result.sizeEns(t,r);
                for e = 1:sizeEns

                    ambiguity = zeros( 1, sizeEns               );
                    dist      = zeros( 1, sizeEns*(sizeEns-1)/2 );
                    fam = util_readFam(nameFam, nBlocks, nReplications, ...
					                   result.sizeEns, r, t, e);
					sizeNn = fam.sizeNn;

                    %-- 1 / (|Wij| + alpha) for all the patterns
					tempNormWa(1:sizeNn)  = transpose(fam.normWalpha(1:sizeNn));
                    normWalpha(:,1:sizeNn)= ...
					               repmat( tempNormWa(1:sizeNn), nPatterns, 1 );

                    %------ Compute T(j) = |min(A,Wij)| / (|Wij| + alpha) -----%
                    %-- |min(A,Wij)| for all patterns and all nodes
                    normAW = zeros(nPatterns, sizeNn);
                    for j = 1:sizeNn
                        %-- Temp - W
                        Wj = repmat(fam.W(j,:), nPatterns, 1);

                        %-- |min(A,Wij)| for all the patterns
                        normAW(:,j) = sum( min(A(:,:),Wj), 2 );
                    end

                    %-- T(j) for all patterns
                    T(:,1:sizeNn) = ...
                        normAW(:,1:sizeNn).*normWalpha(:,1:sizeNn);

                    %-- T_jStar(a)
                    [T_jStar winners] = max(T(:,1:sizeNn),[],2);

                    %-- max(T_j(a), k~=kStar)
                    [tmp kj_one] = max( transpose( fam.Wab ) );
                    kStar_one    = result.predictions(:,e,t,r)+1;

                    kj    = repmat( kj_one, nPatterns, 1);
                    kStar = repmat( kStar_one, 1, sizeNn);

                    filterPred = not( eq(kStar,kj) );                

                    T_j2(:,1:sizeNn)  = T(:,1:sizeNn) .* filterPred;
                    [T_jStar2 second] = max(T_j2(:,1:sizeNn),[],2);

                    ambiguity(e) = sum( T_jStar - T_jStar2 );
                end

                %------------------- Values for the ensemble ------------------%
                for i = 1:sizeEns-1
                    for j = i+1:sizeEns
                        ambiguityEns(r) = abs(ambiguity(i) - ambiguity(j)) + ...
                                          ambiguityEns(r);
                    end
                end
                
                %--------- Diversity in the optimization environment ----------%
                %-- Particles forming the ensemble
                data      = zeros(sizeEns, nDimensions);
                distances = zeros(sizeEns, sizeEns, nBlocks, nReplications);

                for e = 1: sizeEns
                    data(e,:) = pso.p(result.members(e,t,r), :, 1, r, t);
                end

                %-- Computing the distances
                for i = 1:sizeEns-1
                    for j = i+1:sizeEns
                        dst = 0;
                        for d = 1:nDimensions
                            s1=data(i,d);   s2=data(j,d);   dst=dst + (s2-s1)^2;
                        end
                        distances(i,j,t,r) = sqrt(dst);
                    end
                end

                %-- Average value of the distances
                k = 0;
                for i = 1:sizeEns-1
                    for j = i+1:sizeEns
                        k = k+1;   dist(k) = distances(i,j,t,r);
                    end
                end
            end  %-- for : nRep

            %-- Diversity - classifiers
            graphDiv(typeTest*3-2, t) = t;
            graphDiv(typeTest*3-1, t) = mean(ambiguityEns);
            graphDiv(typeTest*3, t)   = std(ambiguityEns) * pValue;

            %-- Diversity - swarm
            graphDivS(typeTest*3-2, t) = t;
            graphDivS(typeTest*3-1, t) = mean(dist);
            graphDivS(typeTest*3, t)   = std(dist) * pValue;

        end  %-- for : nBlocks
    end

    %-- Saving the data
%     nameData = sprintf('../savedStuff/%s_%sInc_diversity.mat', nameDb, nameSc);  
%     save(nameData, 'graphDiv', 'graphDivS', '-mat');
end

if loadResult
    %-- Saving the data
    nameData = sprintf('../savedStuff/%s_%sInc_diversity.mat', nameDb, nameSc);
    load(nameData);
end

%------------------------------------------------------------------------------%
%---------------------------------- Graphics ----------------------------------%
if graph
    %------------------------- TEMPORARY MODIFICATION -------------------------%
    graphDiv(1:6, :) = graphDiv(4:9, :);
    graphDiv(7:9, :) = [];
    graphDivS(1:6,:) = graphDivS(4:9,:);
    graphDivS(7:9,:) = [];
    %------------------------- TEMPORARY MODIFICATION -------------------------%
    
    %-- Initialization
    close all;
    widthErr  = res*errWidth;     heightErr = res*errWidth/ratioErr;

    %-- Diversity FAM
    fig_1 = figure(1);     clf(fig_1);
    posScrn = [0 sizeScrn(2)-heightErr-104 widthErr heightErr];
    set(fig_1, 'Position', posScrn, 'Color', [1,1,1]);
           
    a(1) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.19 0.18 0.79 0.82]);
            
    g = util_graphDiversity(graphDiv,  nTest, 1);
    
    nameFile = sprintf('export_fig ../figures/%s_%s%sDivFAM -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);

            
    %-- Diversity Swarm
    fig_2 = figure(2);     clf(fig_2);
    posScrn = [0 sizeScrn(2)-2*heightErr-187 widthErr heightErr];
    set(fig_2, 'Position', posScrn, 'Color', [1,1,1]);
    a(2) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.18 0.18 0.79 0.82]);
    
	g = util_graphDiversity(graphDivS, nTest, 2); 
    
    nameFile = sprintf('export_fig ../figures/%s_%s%sDivSwm -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);

    %-- Correlation
    fig_3 = figure(3);     clf(fig_3);
    posScrn = [0 sizeScrn(2)-3*heightErr-273 widthErr heightErr];
    set(fig_3, 'Position', posScrn, 'Color', [1,1,1]);
    a(3) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0.21 0.21 0.79 0.79]);
    
    data = [graphDivS(2,:) ; graphDiv(2,:) ; 
            graphDivS(5,:) ; graphDiv(5,:) ];
    g = util_graphDiversity(data, nTest, 3);
    
    nameFile = sprintf('export_fig ../figures/%s_%s%sCorrDiv -pdf', ...
                       nameDb, nameSc, nameTest);
    eval(nameFile);

    %-- Legend
    fig_4 = figure(4);     clf(fig_4);
    posScrn = [sizeScrn(1)+widthErr sizeScrn(2)-heightErr-104 300 30];
    set(fig_4, 'Position', posScrn, 'Color', [1,1,1]);
    a(3) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
                'Position', [0 0 1 1]);

    %-- Text of the legend
    textLegend = [ 'LBESTS$_\mathrm{+d}$' ;
                   'SWARM               ' ];

    %-- Figure
    sizeLine   = 1.5;    sizeLineCI = 0.5;
    sizeMrkr   = 4;      sizeMrkrCI = 2;
    sizeFont   = 8;      sizeFontL  = 10;
    typeMrkr   = 'o^svd<*>';
    hold on
    for t = 1:2
        
        %-- Colors & Markers
        if(t == 1)
            TONTE = 0;
        else
            TONTE = 0.2 + 0.5 * (t-1)/(nTest-1);
        end
        graphColor = [TONTE, TONTE, TONTE]
        graphType   = [typeMrkr(t), ''];

        width = 0.25;
        %-- Line
        plot([2*t-width, 2*t+width], [0 0], '-', 'Color', graphColor, 'LineWidth', sizeLine);
        %-- Marker
        plot(2*t, 0, graphType, 'Color', graphColor, ...
            'LineWidth', sizeLine, 'MarkerSize', sizeMrkr, ...
            'MarkerFaceColor', graphColor);
        txt = text( 2*t +width+0.1, 0, textLegend(t,:), ...
              'FontName','Times New Roman', 'FontSize', sizeFont, ...
              'HorizontalAlignment', 'left', 'Interpreter', 'Latex');
    end

    hold off
    grid off
    box on

    set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', [], ...
             'XLim', [1.70 4.75], 'YLim', [-0.5 0.5]);

    nameFile = sprintf('export_fig ../figures/%s_divLeg -pdf', nameDb);
    eval(nameFile);

end

fprintf('/*----------------- End of results -----------------*/\n');
