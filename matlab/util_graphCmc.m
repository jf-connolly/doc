function  g = util_graphCmc(data, widthFrames, maxFrames)

%-- Graphic's parameters
sizeLine   = 1;      sizeLineCI = 0.5;
sizeMrkr   = 3;      sizeMrkrCI = 1;
sizeFont   = 8;      sizeFontL  = 10;
typeMrkr   = 'o^svh+d<p*>x';

%---------------------------------- Graphic -----------------------------------%
nRank = size(data,2);
nShown = size(data,1)/2;

sp(1) = subplot(121);

%--------------------------------- Error Area ---------------------------------%
% hold on
% for s = 1:nShown
%     %-- Colors & Markers
%     TONTE       = 0.7 * (s-1)/(nShown-1);
%     graphColor  = [TONTE, TONTE, TONTE];
%     graphType   = [typeMrkr(s), '-'];
% 
%     %-- Custumized error area
% %         for r = 1:nRank
% %             negBound = data(2*s-1, r) - data(2*s, r);
% %             posBound = data(2*s-1, r) + data(2*s, r);
%    
%     negBound = data(2*s-1, :) - data(2*s, :);
%     posBound = data(2*s-1, :) + data(2*s, :);
%     
%     xVec = [ 1:nRank,  nRank:-1:1,       1           ];
%     yVec = [ negBound, fliplr(posBound), negBound(1) ];
% 
%     g = fill(xVec,yVec, graphColor+0.2, 'EdgeColor', graphColor, ...
%              'LineWidth', 0.5);
%          
%     %-- Curve (only)                
%     plot(1:nRank, data(2*s-1,:), '-', 'Color', graphColor, ...
%          'LineWidth', sizeLine);
% 
%     %-- Markers
%     for r = 1:nRank
%         plot(r, data(2*s-1, r), graphType, ...
%              'Color', graphColor, 'MarkerFaceColor', graphColor, ...
%              'MarkerSize', sizeMrkr);
%     end
% end

%----------------------------------- Curves -----------------------------------%
for s = 1:nShown
    
    %-- Colors & Markers
    TONTE       = 0.7 * (s-1)/(nShown-1);
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(s), '-'];

    hold on
        %-- Curve (total)
%         g(s) = plot([1:nRank], data(2*s-1,:), ...
%                     graphType, 'Color', graphColor, 'MarkerSize', sizeMrkr, ...
%                     'MarkerFaceColor', graphColor, 'LineWidth', sizeLine, ...
%                     'Visible', 'off');
                
        %-- Curve (curve only)                
        plot(1:nRank, data(2*s-1,:), ...
             '-', 'Color', graphColor, 'LineWidth', sizeLine);

        %-- Custumized error bar
        for r = 1:nRank
            negBound = data(2*s-1, r) - data(2*s, r);
            posBound = data(2*s-1, r) + data(2*s, r);

            plot(r, data(2*s-1, r), graphType, ...
                 'Color', graphColor, 'MarkerFaceColor', graphColor, ...
                 'MarkerSize', sizeMrkr);

            plot([r,r], [negBound, posBound], graphType, ...
                 'LineWidth', sizeLineCI, 'Color', graphColor,...
                 'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
        end
    hold off     
end

minY = min( data(1:2:13,1) - data(2:2:14,1) );

hold off
% set(gca,'FontSize', sizeFont, 'XLim', [1 nRank], 'YLim', [92 100.2], ...
%         'FontName', 'Times New ROman');
set(gca,'FontSize', sizeFont, 'XLim', [1 nRank], 'YLim', [minY-0.1 100.2], ...
        'FontName', 'Times New Roman');
% %         set(gca,'FontSize', sizeFont, 'XLim', [1 nRank], 'YLim', [94.8 100.2]);
% %         set(gca,'FontSize', sizeFont, 'XLim', [1 nRank], 'YLim', [93.3 100.2]);

%-- Labels definition
xlabel('Ranking', ...
       'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
       'FontSize', sizeFont, 'Interpreter','Latex');
ylabel('Cumulative Match (\%)', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

%------------------------------- Legend - Swarm -------------------------------%
sp(2) = subplot(122);

%-- Legend only for the classification rate
textLegend = '';
for s = 1:nShown
    f = s*widthFrames;
    if f<10,  textTemp = sprintf('$\\:\\;%d$', f);
    else      textTemp = sprintf('$%d$', f);
    end
    textLegend(s,1:length(textTemp)) = textTemp;
end

%-- Figure
hold on

    text( -0.3, nShown + 2.4, 'Nb. of ROIs used', ...
          'FontName','Times New Roman', 'FontSize', sizeFont, ...
          'HorizontalAlignment', 'left', 'Interpreter', 'Latex');
    text( -0.3, nShown + 1.3, 'for classification', ...
          'FontName','Times New Roman', 'FontSize', sizeFont, ...
          'HorizontalAlignment', 'left', 'Interpreter', 'Latex');

    for s = 1:nShown
        l = 0.6;
        %-- Colors & Markers
        TONTE       = 0.7 * (s-1)/(nShown-1);
        graphColor  = [TONTE, TONTE, TONTE];
        graphType   = '-';
        graphType   = [typeMrkr(s)];

        %-- Line
        plot([0,l], [s,s], '-', 'Color', graphColor, 'LineWidth', sizeLine);
        %-- Marker
        plot(l/2, s, graphType, 'Color', graphColor, ...
            'LineWidth', sizeLine, 'MarkerSize', sizeMrkr, ...
            'MarkerFaceColor', graphColor);
        %-- Text
%         txt = sprintf('$%d$', shown(t));
        text( l + 0.5, s, textLegend(s,:), ...
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

h    = 0.6/10*nShown;
yPos = (1-h)/2;

set(sp(1), 'Position', [0.12 0.17 0.88 0.83]);
set(sp(2), 'Position', [0.77 0.2 0.2 0.35], ...
    'XLim', [-0.1, l+0.6], 'YLim', [0.4, nShown+0.6]);
