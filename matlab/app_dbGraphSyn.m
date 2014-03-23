%------------------------------------------------------------------------------%
%-- P2 synthetique problem with hold out validation (balanced update scenario)
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
nameDb  = 'p2';
nameSc  = 'upd';
nBlocks =  1;
nameLn  = 'D';
t       =  1;
nameTp  = 't';
r       =  1;

fprintf('\n/*---------- Start - syn db graph -----------*/\n');

%------------------------------ Reading the data ------------------------------%
%-- Reading
if     nameLn == 'T'
    nameFile = sprintf('../database/%s/test.db', nameDb);
elseif nameLn == 'L'
    nameFile = sprintf('../database/%s/learn.db', nameDb);
else
    nameFile = sprintf('../database/%s/%s%dBlocksExt/%s%d^%s%d.db',...
                        nameDb, nameSc, nBlocks, nameLn, t, nameTp, r)
end

dbase = util_readDb(nameFile);

%-- Number of patterns per classes
sizeCls = zeros(1, dbase.nClasses);
for (i=1:dbase.nPatterns)
    sizeCls(dbase.tags(i)+1) = sizeCls(dbase.tags(i)+1) + 1;
end
addrCls = [0 cumsum(sizeCls)];

%-- Printing info
fprintf('\n/*-- Number of patterns - %d\n\n', dbase.nPatterns);
fprintf('/*-- Number of patterns per class - ');
for i = 1:dbase.nClasses
    fprintf('%d ', sizeCls(i));
end

fprintf('\n\n');

%-------------------------------- Scatter plot --------------------------------%
typeMarker = 'osvd<*>';
sizeMk     = 5;
%-- Initialization
fig1 = figure(1);
set(fig1, 'Color', [1,1,1], 'Position', [50 200 400 400],...
                 'PaperPositionMode', 'auto');
clf(fig1);

a = axes('Position', [0 0 1 1]);
set(a, 'XTick', [], 'YTick', [], 'LineWidth', 3);
axis([0 1 0 1], 'square', 'on');
box on;


hold on
    %-- The data
    for c = 1:dbase.nClasses
    %     TONTE       = 0.7 - 0.4 * (c)/(nClasses-1);
        colorMk = ( 0.3 + 0.5 * (c)/(dbase.nClasses) ) *ones(1,3);
        typeMk  = typeMarker(c);

        g = plot(dbase.data(addrCls(c) +1 : addrCls(c+1), 1),...
                 dbase.data(addrCls(c) +1 : addrCls(c+1), 2),...
                 typeMk);

        set(g, 'Color', colorMk, 'MarkerSize', sizeMk,...
                    'MarkerFaceColor', colorMk);
    end
    
    %-- Theoritical class boundaries
    util_graphSynClasses(nameDb);
    
hold off;

%----------------------------------- Saving -----------------------------------%
nameFile = sprintf('../figures/db_%s%s%dBlocks_%s%d^%s%d',...
                    nameDb, nameSc, nBlocks, nameLn, t, nameTp, r);
print(fig1,'-dpng', nameFile);

fprintf('\n/*----------- End - syn db graph ------------*/\n');
