%------------------------------------------------------------------------------%
%-- P2 synthetique problem with hold out validation (balanced update scenario)
%------------------------------------------------------------------------------%

%------------------------------------------------------------------------------%
%---------------------------------- Settings ----------------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
%
nameDb  = 'cis';
nameSc  = 'upd';
nameLn  = 'Inc';
nameHp  = 'EoLB_p'; %- EoMA, EoLB_p
nameAl  = 'multi';
widthMc = 30;
nameTp  = 't';
nBlocks = 1;
nRep    = 1;

t       = 1;
r       = 1;
e       = 1;

nameNn  = 'fam';

readData  = 1;
saveGraph = 1;

boxes      = 0;
boundaries = 1;
compute    = 1;
data       = 1;

resolution = 0.006;      %-- min resolution = 0.007
fprintf('\n/*------------ Start - graph fam ------------*/\n');

%-------------------- Read the network & define the mapping -------------------%
if readData
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.result',...
              nameDb, nameSc, nameLn, nameHp, nameAl, nBlocks, widthMc);
    result = util_readResults(nameFile, nBlocks, nRep);

    %-- Network
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.artmap',...
              nameDb, nameSc, nameLn, nameHp, nameAl, nBlocks, widthMc);
    fam = util_readFam(nameFile, nBlocks, nRep, result.sizeEns, r, t, e);

    %-- Data base
    nameFile = sprintf('../database/%s/%s%dBlocksExt/B%d^%s%d.db',...
                        nameDb, nameSc, nBlocks, t, nameTp, r);
    dbase = util_readDb(nameFile);
end

%-- Number of patterns per classes
sizeCls = zeros(1, dbase.nClasses);
for i=1:dbase.nPatterns
    sizeCls(dbase.tags(i)+1) = sizeCls(dbase.tags(i)+1) + 1;
end
addrCls = [0 cumsum(sizeCls)];

%-- Define the class mapping
if boundaries && compute
    %-- Find the mapping
    [mapFam grid] = util_mappingFam(fam, resolution);
    
    %-- Transform into vectors
    nSteps = floor(1/resolution) +1;
    graphFam = zeros(nSteps^2, 2, fam.nClasses);
    nPoints  = zeros(1,           fam.nClasses);
    
    for posX = 1:nSteps
        for posY = 1:nSteps
            
            pClass          = mapFam(posX,posY);     %- Predicted class
            nPoints(pClass) = nPoints(pClass) +1;

            graphFam(nPoints(pClass), 1, pClass) = grid(posX);
            graphFam(nPoints(pClass), 2, pClass) = grid(posY);
        end
    end
end

%---------------------------------- Graphic -----------------------------------%
nSteps = 1/resolution +1;
typeMk = 's';
sizeMk = 2;
typeMarker = 'o^svd<*>';


%-- Initialization
fig1 = figure(1);   clf(fig1);
set(fig1, 'Color', [1,1,1], 'Position', [50 200 400 400], ...
          'PaperPositionMode', 'auto');

a = axes('Position', [0 0 1 1]);
set(a, 'XTick', [], 'YTick', [], 'LineWidth', 3);
axis([0 1 0 1], 'square', 'on');
box on;

hold on

    %-- The boxes
    if boxes
        sizeNn = fam.sizeNn;
        
        %-- Ranking of the largest boxes
        sizeBoxes = zeros(1, sizeNn);
        for j = 1:sizeNn
            sizeBoxes(j) = (1-fam.W(j,3)-fam.W(j,1)) *(1-fam.W(j,4)-fam.W(j,2));
        end
        [temp index] = sort(sizeBoxes,'descend');
        
        for j = 1:sizeNn
            [tmp ck] = max(fam.Wab(index(j),:));
            colorMk  = ( 0.3 + 0.5 * (ck)/(fam.nClasses) ) *ones(1,3);
            width  = 1-fam.W(index(j),3) -fam.W(index(j),1) +0.0001;
            height = 1-fam.W(index(j),4) -fam.W(index(j),2) +0.0001;

            pos = [fam.W(index(j),1) fam.W(index(j),2) width height];
            rectangle('Position', pos, 'FaceColor', colorMk);

            text( 1-fam.W(index(j),3)-0.03, fam.W(index(j),2)+0.02, ...
                  int2str(index(j)), 'FontName', 'Times New Roman',...
                  'FontSize', 14, 'FontWeight','Bold');
        end
    end
    
    %-- The network
    if boundaries
        for c = 1:fam.nClasses
            if c == 1
                clr = [0.1750    0.1750    0.3000];
%                     clr = [0.2917    0.4167    0.4167];
            else
                clr = [0.5250    0.6500    0.6500];
%                     clr = [0.6458    0.7083    0.7083];
            end
            g = plot( graphFam(1:nPoints(c), 1, c), ...
                      graphFam(1:nPoints(c), 2, c), typeMk, 'Color', clr, ...
                      'MarkerFaceColor', clr, 'MarkerSize', sizeMk);
        end
    end
    
    %-- The data
    if data
        sizeMk = 6;
        for c = 1:dbase.nClasses
            if c == 1;
                clr = [0.3500    0.4125    0.4750];
            else
                clr = [0.7625    0.8250    0.8250];
            end
            typeMk  = typeMarker(c);

            g = plot(dbase.data(addrCls(c) +1 : addrCls(c+1), 1), ...
                     dbase.data(addrCls(c) +1 : addrCls(c+1), 2), ...
                     typeMk, 'MarkerSize', sizeMk,...
                   'MarkerFaceColor', clr, 'MarkerEdgeColor', [0 0 0]);
        end
    end

    %-- Theoritical class boundaries
    util_graphSynClasses(nameDb);
    
    %-- The hyperparameters
    textHp = sprintf('$\\textbf{h}=(%d, %1.2f, %1.2f, %1.2f)$', ...
                     round(fam.h(1)), fam.h(2), fam.h(3), fam.h(4) );
                 
%     txt = text(0.015, 0.955, textHp, 'Interpreter', 'latex', ...
%                'BackgroundColor', [1 1 1], 'EdgeColor', 'none', ...
%                'FontSize', 13, 'FontName', 'Times New Roman', ...    
%                'Margin', 1);

hold off
    
%----------------------------------- Saving -----------------------------------%
if saveGraph
    nameFile = sprintf('../figures/%s_%s%s_fam%d_%d%d', ...
                       nameDb, nameSc, nameHp, e, t, r);

    print('-dpng', '-r300', nameFile);
%     print('-dpdf', nameFile)
end
fprintf('/*------------- End - graph fam -------------*/\n');
