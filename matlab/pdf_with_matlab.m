%------------------------------------------------------------------------------%
%-- How to do graphics with matlab
%-- Comment: Uses the the export_fig function to export a pdf with
%            embedded fonts.
%------------------------------------------------------------------------------%

%-- DEFINE
elsevier = 1;  ieee = 2; lnsc = 3;  thesis = 4;
%-- END DEFINE

%-- Gets screen resolution and gets the resolution in centimeter
temp     = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res      = get(0,'ScreenPixelsPerInch')/2.54;

%------------------------------------------------------------------------------%
%--------------------------- User-defined parameters --------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
typeDoc  = elsevier;                %- See DEFINE
nameFile = 'pdf_with_Matlab';       %- File name
ratio    = 3;                       %- Width/height ratio
axesPos  = [0.07 0.155 0.92 0.82];  %- Axes posiotion within the figure
sizeLine = 2;                       %- Graphic curves size
sizeMrkr = 6;                       %- Markers size
sizeFont = 10;                      %- Font size

%---------------------------- Data for the graphic ----------------------------%
x    = 1:0.2:2*pi;
ySin = sin(x).^2;
yCos = cos(x).^2;

%--------------------------- Figure initialization ----------------------------%
%-- Graphs paper dimension in centimeters
switch typeDoc
    case elsevier
        widthCm = 16.64;
    case ieee
        widthCm = 8.96;
    case lnsc
        widthCm = 12.2;
end
heightCm = widthCm/ratio;

%-- Paper position, NEEDED with the function saveas, but NOT with export_fig
posPaper = [0 0 widthCm heightCm];

%-- Graphs position (upper left corner) and size on screen
width = res*widthCm;   height = width/ratio;
posScrn = [0-width-10 sizeScrn(2)-height-104 width height];

%-- Figure initialization
close all;
fig = figure(1);     clf(fig);
set(fig, 'Position', posScrn, 'PaperUnits', 'centimeters', ...
         'PaperSize', [posPaper(3) posPaper(4)],...
         'PaperPosition', posPaper, 'Color', [1,1,1]);


a = axes('FontName', 'Times New Roman', 'YGrid', 'on', ...
         'Position', axesPos);

%------------------------------- Graphic itself -------------------------------%
%-- Grayscale colors
colorSin = [0.2 0.2 0.2];
colorCos = [0.6 0.6 0.6];

hold on
    g(1) = plot(x, ySin, 'o-', ...
                'Color', colorSin,     'MarkerFaceColor', colorSin, ...
                'LineWidth', sizeLine, 'MarkerSize',      sizeMrkr );

    g(2) = plot(x,yCos, 's-', ...
                'Color', colorCos,     'MarkerFaceColor', colorCos, ...
                'LineWidth', sizeLine, 'MarkerSize',      sizeMrkr );
hold off

%-- Set current axis (see "help gca").
set(gca, 'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
    'XLim', [min(x) max(x)]);

%-- Labels - using the latex interpreter
xlabel('$x$', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');
ylabel('$y$', 'FontName', 'Times New Roman',...
       'FontWeight', 'Bold', 'FontSize', sizeFont, 'Interpreter','Latex');

%-- Legend - using the latex interpreter
textLegend = ['$\textnormal{sin}^2(x)$'; '$\textnormal{cos}^2(x)$'];

leg = legend(g(1:2), textLegend(1:2,:));

set(leg, 'FontSize', sizeFont, 'FontName', 'Times New Roman', ...
         'Orientation', 'Horizontal', 'Location', 'NorthOutside', ...
         'Interpreter', 'Latex');

%----------------------------- Saving the figure ------------------------------%
%-- The usual wy (WITHOU embeded fonts)
saveas(fig, [nameFile, '_with_saveas'], 'pdf');

%-- With export_fig and embeded fonts.
set(0,'CurrentFigure', 1);  %- In case several figures were generated
nameEvaluated = sprintf('export_fig %s%s -pdf', nameFile, '_with_export');
eval(nameEvaluated);

