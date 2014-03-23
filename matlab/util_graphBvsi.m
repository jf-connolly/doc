function  g = util_graphBvsi(data, nTest, typeData, textLegend)


%-------------------------- TO BE MODIFY - TEMPORARY --------------------------%
% data(1:15, :) = data(4:18,:);
% data(16:18,:) = [];
%-------------------------- TO BE MODIFY - TEMPORARY --------------------------%

%-- DEFINE
err=1; cpn=2; ens=3;
%-- END DEFINE

%-- Graphic's parameters
sizeLine   = 1;      sizeLineCI = 0.5;
sizeMrkr   = 4;      sizeMrkrCI = 2;
sizeFont   = 8;      sizeFontL  = 10;
typeMrkr   = 'o^svd<*>os';

if typeData == err,   s(1) = subplot(121);   end

for typeTest = 1:nTest
    
    nBlocks = max( data(typeTest*3-2, :) );

    %-- Colors & Markers
    if(nTest == 1)
        TONTE = 0;
    else
        TONTE = 0.2 + 0.5 * (typeTest-1)/(nTest-1);
    end
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(typeTest), '-'];

    hold on
        %-- Curves (incremental)
        g(typeTest) = plot(data(typeTest*3-2, 1:nBlocks),  ...
                           data(typeTest*3-1, 1:nBlocks),  ...
                           graphType, 'Color', graphColor, ...
                           'MarkerSize', sizeMrkr,         ...
                           'MarkerFaceColor', graphColor,  ...
                           'LineWidth', sizeLine);

        %-- Custumized error bar
        for b = 1:nBlocks
            x        = data(typeTest*3-2, b);
            negBound = data(typeTest*3-1, b) - data(typeTest*3, b);
            posBound = data(typeTest*3-1, b) + data(typeTest*3, b);

            z = plot([x,x], [negBound, posBound], graphType,     ...
                   'LineWidth', sizeLineCI, 'Color', graphColor, ...
                   'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
        end
        set(gca, 'FontSize', sizeFont, 'FontName', 'Times New Roman',...
                 'XLim', [1 nBlocks]);
    hold off     
end

%-- Adjust axes if necessary
if typeData == err && nBlocks == 10
    minData = 69;
    maxData = max(max(data));
    set(gca,'YLim',[minData maxData]);
elseif typeData == ens
    maxData = max(max(data));
    set(gca,'YLim',[0 maxData+1]);
end

%-- Zoom in
if typeData == err
    s(2) = subplot(122);
    for typeTest = 1:nTest
        
        nBlocks = max( data(typeTest*3-2, :) );

        %-- Colors & Markers
        if(nTest == 1)
            TONTE = 0;
        else
            TONTE = 0.2 + 0.5 * (typeTest-1)/(nTest-1);
        end
        graphColor  = [TONTE, TONTE, TONTE];
        graphType   = [typeMrkr(typeTest), '-'];

        hold on
            %-- Curves (incremental)
            plot( data(typeTest*3-2, nBlocks-2:nBlocks), ...
                  data(typeTest*3-1, nBlocks-2:nBlocks), ...
                  graphType, 'Color', graphColor, 'LineWidth', sizeLine, ...
                  'MarkerSize', sizeMrkr, 'MarkerFaceColor', graphColor );
                           
            %-- Custumized error bar
            for b = nBlocks-2:nBlocks
                x        = data(typeTest*3-2, b);
                limData(typeTest,b) = data(typeTest*3-1, b);
                negBound = data(typeTest*3-1, b) - data(typeTest*3, b);
                posBound = data(typeTest*3-1, b) + data(typeTest*3, b);

                z = plot([x,x], [negBound, posBound], graphType, ...
                       'LineWidth', sizeLineCI, 'Color', graphColor,...
                       'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
            end
            set(gca,'FontSize', sizeFont-2, 'FontName', 'Times New Roman');
        hold off     
    end

    if nBlocks == 12
        set(s(2), 'Position', [0.63 0.25 0.12 0.35]);
        xtickl = ['11'; '12'];
        set(s(1), 'YLim', [40 85]);
    else
        set(s(2), 'Position', [0.635 0.625 0.12 0.38]);
        xtickl = [' 9'; '10'];
        set(s(1), 'YLim', [73 100]);
    end

    limData(3,:) = [];
    set(s(1), 'YGrid', 'on', 'Position', [0.11 0.165 0.87 0.83]);
    set(s(2), 'XLim', [nBlocks-1.5 nBlocks], 'XTick', [nBlocks-1 nBlocks], ...
        'XTickLabel', xtickl, 'YGrid', 'on', ...
        'YLim', [min(min(limData(:,nBlocks-1:nBlocks))) - 0.6, ...
                 max(max(limData(:,nBlocks-1:nBlocks))) + 0.6] );
end

%-- Y Label definition
if typeData == err
    labely = 'Classification rate (\%)';
elseif typeData == cpn
    labely = 'Compression';
elseif typeData == ens
    labely = 'Ensemble size';
end    

if typeData == err
    xlabel(s(1),'Learning data set $D_t$', 'FontName', 'Times New Roman',...
           'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');
    ylabel(s(1),labely, 'FontName', 'Times New Roman',...
           'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');
else
    xlabel('Learning data set $D_t$', 'FontName', 'Times New Roman',...
           'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');
    ylabel(labely, 'FontName', 'Times New Roman',...
           'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');
end
%-- Legend only for the classification rate
if typeData == err
    leg = legend(g(1:nTest), textLegend(1:nTest,:));
    set(leg, 'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
             'Orientation', textLegend(nTest+2,:), 'Interpreter','Latex', ...
             'Location',    textLegend(nTest+1,:) );
         bp = 1;
end