function  g = util_graphVideoB(data, typeTest)

%-- Graphic's parameters
sizeLine   = 1;      sizeLineCI = 0.5;
sizeMrkr   = 3;      sizeMrkrCI = 1;
sizeFont   = 8;      sizeFontL  = 10;
typeMrkr   = 'o^svh+d<p*>x';

%---------------------------------- Graphic -----------------------------------%
s(1) = subplot(121);
nBlocks = size(data,3);

nShown = 4;
shown  = round( 0:(nBlocks)/(nShown-1):nBlocks );
shown(1) = 1;

nPatterns = size(data,2);
% nPatterns = 50;
x = 1:nPatterns;

nMarkers = 10;
marked = round( 1:(nPatterns-1)/(nMarkers-1):nPatterns );

hold on
%--------------------------------- Error Area ---------------------------------%
% for t = 1:nShown
%     dataTemp = data(:,:,shown(t));
% 
%     %-- Colors & Markers
%     TONTE       = 0.7 * (t-1)/(nShown-1);
%     graphColor  = [TONTE, TONTE, TONTE];
%     graphType   = [typeMrkr(t), '-'];
% 
%     %-- Custumized error area
%     negBound = dataTemp(typeTest*2-1, :) - dataTemp(typeTest*2, :);
%     posBound = dataTemp(typeTest*2-1, :) + dataTemp(typeTest*2, :);
% 
%     xVec = [ 1:nPatterns, nPatterns:-1:1,   1           ];
%     yVec = [ negBound,    fliplr(posBound), negBound(1) ];
% 
%     g = fill(xVec,yVec, graphColor+0.2, 'EdgeColor', graphColor, ...
%              'LineWidth', 0.5);
%          
%     %-- Curve (only)                
%     g(t) = plot(1:nPatterns, dataTemp(typeTest*2-1, 1:nPatterns), ...
%          '-', 'Color', graphColor, 'LineWidth', sizeLine);
% 
%     %-- Markers
%     for p = 1:nMarkers
%         xm = marked(p);
%         
%         plot(xm, dataTemp(typeTest*2-1, xm), graphType, ...
%              'Color', graphColor, 'MarkerFaceColor', graphColor, ...
%              'MarkerSize', sizeMrkr);
%     end
% end

%----------------------------------- Curves -----------------------------------%
for t = 1:nShown
    dataTemp = data(:,:,shown(t));

    %-- Colors & Markers
    TONTE       = 0.7 * (t-1)/(nShown-1);
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(t), '-'];

    %-- Curve (only)                
    g(t) = plot(1:nPatterns, dataTemp(typeTest*2-1, 1:nPatterns), ...
         '-', 'Color', graphColor, 'LineWidth', sizeLine);

    %-- Custumized error bar
    for p = 1:nMarkers
        xm = marked(p);
        
        %-- Markers
        plot(xm, dataTemp(typeTest*2-1, xm), graphType, ...
             'Color', graphColor, 'MarkerFaceColor', graphColor, ...
             'MarkerSize', sizeMrkr);
         
        %-- Error bars
        negBound = dataTemp(typeTest*2-1, xm) - dataTemp(typeTest*2, xm);
        posBound = dataTemp(typeTest*2-1, xm) + dataTemp(typeTest*2, xm);

        plot([xm,xm], [negBound, posBound], graphType, ...
             'LineWidth', sizeLineCI, 'Color', graphColor,...
             'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
    end
end

hold off     
ylim = get(gca, 'YLim');
% ylim = [0 20];
set(gca,'YLim', [0 ylim(2)], 'XLim', [1 nPatterns], ...
        'FontName', 'Times New Roman', 'FontSize', sizeFont);

%-- Labels definition
xlabel('Number of detected ROIs used for classification', ...
       'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
       'FontSize', sizeFont, 'Interpreter','Latex');
ylabel('Video-based error rate (\%)', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

%------------------------------- Legend - Swarm -------------------------------%
s(2) = subplot(122);

%-- Legend only for the classification rate
textLegend = '';
for t = 1:nShown
    if shown(t)<10,    textTemp = sprintf('$t=\\:\\;%d$', shown(t));
    else               textTemp = sprintf('$t=%d$', shown(t));
    end
    textLegend(t,1:length(textTemp)) = textTemp;
end

% leg = legend(g, textLegend);
% set(leg, 'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
%          'Orientation', 'vertical', 'Interpreter','Latex', ...
%          'Location', 'SouthEast' );

%-- Figure
hold on
    for t = 1:nShown
        l = 0.4;
        %-- Colors & Markers
        TONTE       = 0.7 * (t-1)/(nShown-1);
        graphColor  = [TONTE, TONTE, TONTE];
        graphType   = '-';
        graphType   = [typeMrkr(t)];

        %-- Line
        plot([0,l], [t,t], '-', 'Color', graphColor, 'LineWidth', sizeLine);
        %-- Marker
        plot(l/2, t, graphType, 'Color', graphColor, ...
            'LineWidth', sizeLine, 'MarkerSize', sizeMrkr, ...
            'MarkerFaceColor', graphColor);
        %-- Text
%         txt = sprintf('$%d$', shown(t));
        text( l + 0.7, t, textLegend(t,:), ...
              'FontName','Times New Roman', 'FontSize', sizeFont, ...
              'HorizontalAlignment', 'right', 'Interpreter', 'Latex');
    end
    
%     t = text( 0.1, nShown+1, '$D_t$', ...
%       'FontName','Times New Roman', 'FontSize', sizeFont, ...
%       'HorizontalAlignment', 'center', 'Interpreter', 'latex');

hold off
grid off
box on

set(gca, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);

h    = 0.6/10*nBlocks;
yPos = (1-h)/2;

set(s(1), 'Position', [0.105 0.17 0.90 0.83]);
set(s(2), 'Position', [0.81 0.77 0.2 0.25], ...
    'XLim', [-0.1, l+0.8], 'YLim', [0.5, nShown+0.5]);
