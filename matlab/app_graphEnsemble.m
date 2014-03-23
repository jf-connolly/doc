%------------------------------------------------------------------------------%
%-- Graph of the decision boundaries for an ensemble of fuzzy ARTMAP
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
nameHp  = 'EoMA'; %- EoMA, EoLB_p
nameAl  = 'multi';
widthMc = 30;
nameTp  = 't';
nBlocks = 1;
nRep    = 1;
t       = 1;
r       = 1;

readData   = 1;
boundaries = 1;
data       = 1;

nameNn     = 'fam';
resolution = 0.007;      %-- min resolution = 0.007
fprintf('\n/*---------- Start - graph ensemble ---------*/\n');

%-------------------- Read the network & define the mapping -------------------%
if readData
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.result',...
              nameDb, nameSc, nameLn, nameHp, nameAl, nBlocks, widthMc);
    result = util_readResults(nameFile, nBlocks, nRep);
    sizeF2 = sum(sum(result.sizeClasses(1:2,1:3)))
    
    %-- Network
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.artmap',...
              nameDb, nameSc, nameLn, nameHp, nameAl, nBlocks, widthMc);

    sizeEns = result.sizeEns;
    sizeEns = 3
    for e = 1:sizeEns
        fam(e) = util_readFam(nameFile, nBlocks, nRep, sizeEns, r, t, e);
    end

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
mapFam = zeros((1/resolution) +1, (1/resolution) +1, sizeEns);
for e = 1:sizeEns
    [temp grid]   = util_mappingFam(fam(e), resolution);
    mapFam(:,:,e) = temp;
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

    %-- The ensemble
    if boundaries
        mapEns = round(mean(mapFam,3));
        for posX = 1:nSteps
            for posY = 1:nSteps

                %-- The voting process (using the mean - good for 2 classes)
                ck      = mapEns(posX, posY);
                if ck == 1
                    clr = [0.1750    0.1750    0.3000];
%                     clr = [0.2917    0.4167    0.4167];
                else
                    clr = [0.5250    0.6500    0.6500];
%                     clr = [0.6458    0.7083    0.7083];
                end
                
                g = plot( grid(posX), grid(posY), typeMk, 'Color', clr, ...
                          'MarkerSize', sizeMk, 'MarkerFaceColor', clr);
            end
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

hold off
    
%----------------------------------- Saving -----------------------------------%
nameFile = sprintf('../figures/%s_%s%s_ens_%d%d', nameDb, nameSc, nameHp, t, r);
print('-dpng', '-r300', nameFile);
% saveas(fig1, nameFile, 'png');

fprintf('/*----------- End - graph ensemble ----------*/\n');
