function  g = util_graphDiversity(data, nTest, typeData)

%-- DEFINE
divFam=1; divSwarm=2; corr=3;
%-- END DEFINE

%-- Graphic's parameters
    sizeLine   = 1.5;    sizeLineCI = 0.5;
    sizeMrkr   = 4;      sizeMrkrCI = 2;
    sizeFont   = 8;      sizeFontL  = 10;
    typeMrkr   = 'o^svd<*>';

if typeData ~= corr

    for typeTest = 1:nTest

        nBlocks = max( data(typeTest*3-2, :) );

        %-- Colors & Markers
        if(nTest == 1)
            TONTE = 0;
        else
            TONTE = 0 + 0.7 * (typeTest-1)/(nTest-1);
        end
        graphColor  = [TONTE, TONTE, TONTE]
        graphType   = [typeMrkr(typeTest), '-'];

        hold on
            %-- Curves (incremental)
            g(typeTest) = plot(data(typeTest*3-2, 1:nBlocks), ...
                               data(typeTest*3-1, 1:nBlocks), ...
                               graphType, 'Color', graphColor, ...
                               'MarkerSize', sizeMrkr, ...
                               'MarkerFaceColor', graphColor, ...
                               'LineWidth', sizeLine);

            %-- Custumized error bar
            for b = 1:nBlocks
                x        = data(typeTest*3-2, b);
                negBound = data(typeTest*3-1, b) - data(typeTest*3, b);
                posBound = data(typeTest*3-1, b) + data(typeTest*3, b);

                z = plot([x,x], [negBound, posBound], graphType, ...
                       'LineWidth', sizeLineCI, 'Color', graphColor,...
                       'MarkerFaceColor', graphColor, 'MarkerSize', sizeMrkrCI);
            end
            set(gca,'FontSize', sizeFont, 'XLim', [1 nBlocks]);

        hold off     
    end
else
    for typeTest = 1:nTest
        %-- Colors & Markers
        if(nTest == 1)
            TONTE = 0;
        else
            TONTE = 0 + 0.5 * (typeTest-1)/(nTest-1);
        end
        graphColor = [TONTE, TONTE, TONTE];
        graphType   = [typeMrkr(typeTest), ''];

        %-- Original data
        x = data(2*typeTest-1,:);    y = data(2*typeTest,:);

        %-- Regression (polynomial or exponential)
        xg = min(x)-0.05:0.0001:max(x)+0.05;
%         p = polyfit(x,y,2)
%         f = polyval(p,xg);
        pExp = polyfit(x, log(y), 1);
        f = exp(pExp(2)).*exp(pExp(1).*xg);
        
        hold on
            g = plot(x, y, graphType, 'Color', graphColor, ...
                     'MarkerSize', sizeMrkr, 'MarkerFaceColor', graphColor, ...
                     'MarkerSize', sizeMrkr);
                 
            h = plot(xg, f, ':', 'Color', graphColor, 'LineWidth', 1);
        hold off
        
        if size(data,2) == 12,   set(gca,'XLim', [0.5 0.95]);
        else                     set(gca,'XLim', [0.35 0.9]);
        end
    end
end

%-- Adjust axes if necessary
% if typeData == ens
%     maxData = max(max(data));
%     set(gca,'YLim',[0 maxData+1]);
% end
% ax = axis;
% ax = [min(data(1,1:nBlocks)), max(data(1,1:nBlocks)), ax(3), ax(4)];
% axis(ax);

%-- Y Label definition
if     typeData == divFam
    labelx = 'Learning data set $D_t$';
    labely = '$\overline{\Delta\theta_{e_1e_2}}$';
elseif typeData == divSwarm
    labelx = 'Learning data set $D_t$';
    labely = '$\overline{\delta_{e_1e_2}}$';
elseif typeData == corr
    labelx = '$\overline{\delta_{e_1e_2}}$';
    labely = '$\overline{\Delta\theta_{e_1e_2}}$';
end    

xlabel(labelx, 'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
               'FontSize', sizeFont, 'Interpreter','Latex');
ylabel(labely, 'FontName', 'Times New Roman', 'FontWeight', 'Bold', ...
               'FontSize', sizeFont, 'Interpreter','Latex');
   
