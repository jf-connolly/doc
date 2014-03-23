function  g = util_graphCmcTime(data, nFrames, nTest)

%-- Graphic's parameters
sizeLine   = 1;      sizeLineCI = 0.5;
sizeMrkr   = 5;      sizeMrkrCI = 2;
sizeFont   = 8;      sizeFontL  = 10;
typeMrkr   = 'vs^oh+d<p*>x';

%---------------------------------- Graphic -----------------------------------%
nRank  = size(data,2);
nShown = nTest;

sp(1) = subplot(121);

%----------------------------------- Curves -----------------------------------%
maxY = 0;
minY = 100;
for tTest = nShown:-1:1
    addr     = (2*tTest)-1;
    dataTemp = data(addr:addr+1,:);
    
    %-- Colors & Markers
    TONTE       = 0.7 * (tTest-1)/(nShown-1);
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(nShown+1-tTest), '-'];

    hold on
    %-- Curve (only)                
    g(tTest) = plot(1:nRank, dataTemp(1, :), ...
         '-', 'Color', graphColor, 'LineWidth', sizeLine);

    %-- Custumized error bar
    for r = 1:nRank
        negBound = dataTemp(1, r) - dataTemp(2, r);
        posBound = dataTemp(1, r) + dataTemp(2, r);

        plot(r, dataTemp(1, r), graphType, ...
             'Color', graphColor, 'MarkerFaceColor', graphColor, ...
             'MarkerSize', sizeMrkr);

        plot([r,r], [negBound, posBound], graphType, ...
             'LineWidth', sizeLineCI, 'Color', graphColor,...
             'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
    end    
    
    minY = min([minY dataTemp(1,:)-dataTemp(2,:)]);
end

hold off     

set(gca,'YLim', [minY-1 100.1], 'XLim', [1 nRank], ...
        'FontName', 'Times New Roman', 'FontSize', sizeFont);

%-- Labels definition
xlabel('Ranking', ...
       'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
       'FontSize', sizeFont, 'Interpreter','Latex');
ylabel('Cumulative Match (\%)', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

%------------------------------- Legend - Swarm -------------------------------%
sp(2) = subplot(122);

%-- Legend only for the classification rate
textLegend = ['ADNPSO'; ' DNPSO'; ' MOPSO'; ' GBEST'];

%-- Figure
hold on
    for t = nShown:-1:1
        l = 0.3;
        %-- Colors & Markers
        TONTE       = 0.7 * (t-1)/(nShown-1);
        graphColor  = [TONTE, TONTE, TONTE];
        graphType   = [typeMrkr(nShown-t+1)];

        %-- Line
        plot([0,l], nShown+1-[t,t], '-', 'Color', graphColor, ...
                    'LineWidth', sizeLine);
        %-- Marker
        plot(l/2, nShown+1-t, graphType, 'Color', graphColor, ...
            'LineWidth', sizeLine, 'MarkerSize', sizeMrkr, ...
            'MarkerFaceColor', graphColor);
        %-- Text
        text( l + 0.92, nShown+1-t, textLegend(t,:), ...
              'FontName','Times New Roman', 'FontSize', sizeFont, ...
              'HorizontalAlignment', 'right', 'Interpreter', 'Latex');
    end
    
hold off
grid off
box on

set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

set(sp(1), 'Position', [0.18 0.18 0.81 0.82]);
set(sp(2), 'Position', [0.65 0.20 0.35 0.25], ...
    'XLim', [-0.1, l+1], 'YLim', [0.5, nShown+0.5]);
