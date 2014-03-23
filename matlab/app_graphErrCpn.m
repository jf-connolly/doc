%------------------------------------------------------------------------------%
%-- Results: batch -vs- inremental for different settings
%-- Works with util_getName
%------------------------------------------------------------------------------%
fprintf('\n/*---------------- Results - static ----------------*/\n');

%-- DEFINE
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
time           =  3;
nameTest       = 'total';%   {'Ens','Vid','tst'}
nTest          =  4;
legLocation    = 'EastOutside';
legOrientation = 'vertical';
nReplications  =  50;

compute    = 1;
saveResult = 0;
loadResult = 0;
graph      = 1;
show       = 1;
hypTest    = 0;

%-- Graphs
errG = 1;
cpnG = 1;
ensG = 1;

%-- Graph width & ratio
typeDoc  = elsevier; %  elsevier;
ratioErr = 1.25; %  {1.75, 2.5};
ratioOtr = 2;

%---------------------------- Graphics parameters -----------------------------%
%-- Graphs width
switch typeDoc
    case elsevier
        errWidth = 16.64/3;
    case ieee
        errWidth = 8.96;
    case lnsc
        errWidth = 12.2;
end
otherWidth = errWidth/2.1;

%------------------------------------------------------------------------------%
%------------------------------- Initialization -------------------------------%
close all;
widthErr  = res*errWidth;     heightErr = res*errWidth/ratioErr;
widthOtr  = res*otherWidth;   heightOtr = res*otherWidth/ratioOtr;

%-- Classification rate 
fig_1 = figure(1);     clf(fig_1);
posScrn = [0 sizeScrn(2)-heightErr-104 widthErr heightErr];
set(fig_1, 'Position', posScrn, 'Color', [1,1,1]);

a(1) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
            'Position', [0.17 0.19 0.82 0.80]);

fig_2 = figure(2);     clf(fig_2);
posScrn = [0 sizeScrn(2)-heightErr*2-104 widthErr heightErr];
set(fig_2, 'Position', posScrn, 'Color', [1,1,1]);

a(2) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
            'Position', [0.17 0.19 0.82 0.80]);

%--Legend
ratioErr = 10
widthErr  = res*errWidth;     heightErr = res*errWidth/ratioErr;

fig_3 = figure(3);     clf(fig_3);
posScrn = [0 sizeScrn(2)-heightErr*2-104 widthErr*2 heightErr];
set(fig_3, 'Position', posScrn, 'Color', [1,1,1]);

a(3) = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
            'Position', [0 0 1 1]);

%------------------------------------------------------------------------------%
%---------------------- Calcul des points du graphiques -----------------------%
%------------------------------------------------------------------------------%
%-- Order: APSO, DNPSO, MOSPO, PSO, kNN
%-- Test cases - enrollment: 1 5 10
%-- Test cases - update:     1 6 12

errAdd = zeros(2,6,3);
cpnAdd = zeros(2,6,3);

errUpd = zeros(2,6,3);
cpnUpd = zeros(2,6,3);

errAdd(:,:,1) = [ 1.0 0 1.0 10 0 0 ; ...
                  1.7 0 1.7  7 0 0 ];
cpnAdd(:,:,1) = [ 46 0.37 67 54 0.04 1 ; ...
                   9 0.02  8  8 0.01 0 ];
               
errAdd(:,:,2) = [ 1.0 0 6 2.0 0 0 ; ...
                  0.9 0 2 0.6 0 0 ];
cpnAdd(:,:,2) = [ 2.3 0.35 12 13 0.04 1 ; ...
                  0.2 0.02  2  5 0.01 0 ];

errAdd(:,:,3) = [ 0.6 0 5.2 0.8 0 0 ; ...
                  0.7 0 2 0.5 0 0 ];
cpnAdd(:,:,3) = [ 2.1 0.34 8.8 5.6 0.037 1; ...
                  0.2 0.02 0.9 0.7 0.003 0];

errUpd(:,:,1) = [ 24 20 24 24 21 18 ; ...
                   2  1  2  2  1  2 ];
cpnUpd(:,:,1) = [ 3.3 0.17 3.0 3.3 0.04 1; ...
                  0.7 0.01 0.8 0.6 0.01 0];

errUpd(:,:,2) = [ 10 11 10 11 10 10 ; ...
                   1  1  1  2  1  1 ];
cpnUpd(:,:,2) = [ 1.7 0.28 3 4.7 0.04 1 ; ...
                  0.2 0.01 2 0.4 0.01 0 ];

errUpd(:,:,3) = [ 0 0 5 0.8 0 0; ...
                  0 0 2 0.5 0 0 ];
cpnUpd(:,:,3) = [ 1.4 0.30 4 6 0.037 1; ...
                  0.2 0.03 1 3 0.037 0 ];

realTime = [1 5 10; 1 6 12];


%-- Graphic's parameters
sizeLine   = 2;      sizeLineCI = 0.5;
sizeMrkr   = 6;      sizeMrkrCI = 2;
sizeFont   = 10;      sizeFontL  = 10;
typeMrkr   = 'o^svd<*>os';


for tTest = 1:nTest
    
	name = util_getName('total', tTest);

    %-- Colors & Markers
    if(nTest == 1)
        TONTE = 0;
    else
        TONTE = 0.2 + 0.5 * (tTest-1)/(nTest-1);
    end
    graphColor  = [TONTE, TONTE, TONTE];
    graphType   = [typeMrkr(tTest), '-'];

    set(0,'CurrentFigure', 1);
    hold on 
        %-- Point
        g = plot(cpnAdd(1, tTest, time), errAdd(1, tTest, time), ...
                 graphType, 'Color', graphColor, 'MarkerSize', sizeMrkr, ...
                 'MarkerFaceColor', graphColor, 'LineWidth', sizeLine);
             
        %-- IC
        pos = [ cpnAdd(1, tTest, time) - cpnAdd(2, tTest, time), ...
                errAdd(1, tTest, time) - errAdd(2, tTest, time), ...
                cpnAdd(2, tTest, time) * 2, ...
                errAdd(2, tTest, time) * 2 ];
             
        if pos(3) && pos(4)
            h = rectangle('Position', pos, 'EdgeColor', graphColor, ...
                         'LineWidth', sizeLine);
        end  
        
        annotation('arrow', 'Position', [0.2 0.83 0.15 -0.15], 'LineWidth', 2);
        annotation('textbox', 'Position', [0.2 1 0.15 -0.15], ...
                   'String', 'Good solutions', 'EdgeColor', 'none',...
                   'VerticalAlignment', 'Baseline','Interpreter','Latex');
    hold off

    set(0,'CurrentFigure', 2);
    hold on 
        %-- Curves (incremental)
        g = plot(cpnUpd(1, tTest, time), errUpd(1, tTest, time), ...
                 graphType, 'Color', graphColor, 'MarkerSize', sizeMrkr, ...
                 'MarkerFaceColor', graphColor, 'LineWidth', sizeLine);
             
        %-- IC
        pos = [ cpnUpd(1, tTest, time) - cpnUpd(2, tTest, time), ...
                errUpd(1, tTest, time) - errUpd(2, tTest, time), ...
                cpnUpd(2, tTest, time) * 2, ...
                errUpd(2, tTest, time) * 2 ];
             
        if pos(3) && pos(4)
            h = rectangle('Position', pos, 'EdgeColor', graphColor, ...
                         'LineWidth', sizeLine);
        end     
        
        annotation('arrow', 'Position', [0.2 0.83 0.15 -0.15], 'LineWidth', 2);
        annotation('textbox', 'Position', [0.2 1 0.15 -0.15], ...
                   'String', 'Good solutions', 'EdgeColor', 'none',...
                   'VerticalAlignment', 'Baseline','Interpreter','Latex');
    hold off

    set(0,'CurrentFigure', 3);
    hold on
        plot(tTest,1, graphType, 'Color', graphColor, ...
             'MarkerSize', sizeMrkr+2, 'MarkerFaceColor', graphColor); 
        text(tTest+0.2, 1, name.legend)
    hold off
    bp=1;
end

set(0,'CurrentFigure', 1);
ylabel('Error rate (\%)', 'FontName', 'Times New Roman',...
       'FontSize', sizeFont, 'Interpreter','Latex');
xlabel('Compression', 'FontName', 'Times New Roman',...
       'FontSize', sizeFont, 'Interpreter','Latex');
set(gca,'FontSize', sizeFont-2, 'FontName', 'Times New Roman', ...
         'XGrid', 'on');
nameFile = sprintf('export_fig ../figures/cnrc64_addErrCpn_%d -pdf', time);
eval(nameFile);

set(0,'CurrentFigure', 2);
ylabel('Error rate (\%)', 'FontName', 'Times New Roman',...
       'FontSize', sizeFont, 'Interpreter','Latex');
xlabel('Compression', 'FontName', 'Times New Roman',...
       'FontSize', sizeFont, 'Interpreter','Latex');
set(gca,'FontSize', sizeFont-2, 'FontName', 'Times New Roman', ...
         'XGrid', 'on');
nameFile = sprintf('export_fig ../figures/cnrc64_updErrCpn_%d -pdf', time);
eval(nameFile);

set(0,'CurrentFigure', 3);
grid off
box on
set(gca, 'XLim', [0.8 4.75], 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
nameFile = sprintf('export_fig ../figures/cnrc64_updErrCpn_leg -pdf', time);
eval(nameFile);
     

fprintf('/*----------------- End of results -----------------*/\n');
