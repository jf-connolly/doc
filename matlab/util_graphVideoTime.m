function  g = util_graphVideoTime(data, t, nTest)

%-- Graphic's parameters
sizeLine   = 1;      sizeLineCI = 0.5;
sizeMrkr   = 3;      sizeMrkrCI = 1;
sizeFont   = 8;      sizeFontL  = 10;
typeMrkr   = 'vs^oh+d<p*>x';

%---------------------------------- Graphic -----------------------------------%
s(1) = subplot(121);

nShown    = nTest;
nPatterns = size(data,2);
% nPatterns = 50;
x         = 1:nPatterns;

nMarkers = 10;
marked = round( 1:(nPatterns-1)/(nMarkers-1):nPatterns );

hold on

%----------------------------------- Curves -----------------------------------%
maxY = 0;
minY = 100;
for tTest = nShown:-1:1
    addr     = (2*tTest)-1;
    dataTemp = data(addr:addr+1,:,t);

    %-- Colors & Markers
    TONTE       = 0.7 * (tTest-1)/(nShown-1);
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(nShown+1-tTest), '-'];

    %-- Curve (only)                
    g(t) = plot(x, dataTemp(1, 1:nPatterns), ...
         '-', 'Color', graphColor, 'LineWidth', sizeLine);

    %-- Custumized error bar
    for p = 1:nMarkers
        xm = marked(p);
        
        %-- Markers
        plot(xm, dataTemp(1, xm), graphType, ...
             'Color', graphColor, 'MarkerFaceColor', graphColor, ...
             'MarkerSize', sizeMrkr);
         
        %-- Error bars
        negBound = dataTemp(1, xm) - dataTemp(2, xm);
        posBound = dataTemp(1, xm) + dataTemp(2, xm);

        plot([xm,xm], [negBound, posBound], graphType, ...
             'LineWidth', sizeLineCI, 'Color', graphColor,...
             'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
    end
    
    maxY = max([maxY dataTemp(1,:)+dataTemp(2,:)]);
    minY = min([minY dataTemp(1,:)-dataTemp(2,:)]);
    
end

hold off     
ylim = get(gca, 'YLim');
% ylim = [0 20];
set(gca,'YLim', [0 31], 'XLim', [1 nPatterns+1], ...
        'FontName', 'Times New Roman', 'FontSize', sizeFont);

%-- Labels definition
xlabel('Number of ROIs', ...
       'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
       'FontSize', sizeFont, 'Interpreter','Latex');
ylabel('Video-based error rate (\%)', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

%------------------------------- Legend - Swarm -------------------------------%
s(2) = subplot(122);

%-- Legend only for the classification rate
textLegend = ['ADNPSO'; ' DNPSO'; ' MOPSO'; ' GBEST'];

%-- Figure
hold on
    for t = nShown:-1:1
        l = 0.3;
        
        %-- Colors & Markers
        TONTE       = 0.7 * (t-1)/(nShown-1);
        graphColor  = [TONTE, TONTE, TONTE];
        graphType   = '-';
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

set(s(1), 'Position', [0.18 0.18 0.81 0.82]);
set(s(2), 'Position', [0.65 0.77 0.35 0.25], ...
    'XLim', [-0.1, l+1], 'YLim', [0.5, nShown+0.5]);
