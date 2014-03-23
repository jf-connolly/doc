%------------------------------------------------------------------------------%
%-- Scenario : Class enrollment without LTM
%------------------------------------------------------------------------------%

%---------------------------- Parametres variables ----------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
nameDb    = 'mpb';
nPeaks    = 40;       
nShifts   = 4;                %-- Number of blocs

nDims     = 2;
normV     = 0.2*ones(1,nDims);
s         = randn(1,nDims);   %-- sigma
l         = 0.2;              %-- lambda
sevHeight = 0.05;
sevWidth  = 0.05;
new   = 0;
graph = 1;
dbase = 0;

close all

fprintf('\n/*------ Start - moving peaks benchmark ------*/\n');

%--------------------------------- Parameters ---------------------------------%
if new
    %-- Initialization
    rand('twister',sum(100*clock));
    height   = zeros(nPeaks, nShifts);
    width    = zeros(nPeaks, nShifts);
    posPeaks = zeros(nPeaks, nDims, nShifts);
    vCurr    = zeros(nPeaks, nDims);
    vPrev    = zeros(nPeaks, nDims);

    for t = 1:nShifts

        if t == 1
            for p = 1:nPeaks
                r          = rand(1,nDims);         %-- Random vector
                normR      = sqrt( sum( r.^2 ));    %-- To normalize
                vCurr(p,:) = normV/normR .* r;      %-- Peaks velocity

                height  (p,t)   = 1   + rand(1,1);
                width   (p,t)   = 0.3 + rand(1,1)/5;
                posPeaks(p,:,t) = rand(1,nDims);
            end

        else
            for p = 1:nPeaks
                %-- Height and width
                height  (p,t)   = height(p,t-1) + randn(1,1) * sevHeight;
                width   (p,t)   = width (p,t-1) + randn(1,1) * sevWidth;

                %-- Peaks velocity
                vPrev(p,:) = vCurr(p,:);                       %-- Previous upd.
                
                r          = rand(1,nDims)*2 -1;             %-- Random vector
                normR      = sqrt( sum( (r + vPrev(p,:)).^2 ));%-- To normalize
                
                vCurr(p,:) = normV ./ normR .* ((1-l)*r + l*vPrev(p,:));
                    
                %-- Peaks position update
                posPeaks(p,:,t) = posPeaks(p,:,t-1) + vCurr(p,:);
            end
        end
    end

    %-- Min-Max normalization of the peaks positions
    for d = 1:nDims
        posPeaks(:,d,:) = ( posPeaks(:,d,:) - min(min(posPeaks(:,d,:))) ) / ...
               ( max(max(posPeaks(:,d,:))) - min(min(posPeaks(:,d,:))) );
    end
end

%------------------------------ Ploting the peaks -----------------------------%
if graph
    
    resolution = 10;
    x = 0 :1/resolution :1;
    f = zeros(resolution, resolution);
    fpeak = zeros(1, nPeaks);
    normS = zeros(1, nPeaks);

    for t = 1:nShifts

        %-- Find the fitnesses
        for i = 1: size(x,2);
            for j = 1: size(x,2);
                posCurr = [x(i) x(j)];

                for p = 1:nPeaks
                    posPeaks(p,:,t);
                    normS(p) = sum( (posCurr - posPeaks(p,:,t)).^2 );
                    fpeak(p) = height(p,t) * exp( - normS(p) / width(p,t)^2 );
                end
                
                f(i,j) = max(fpeak);   %-- Because maximization
            end
        end
        
        fig(t) = figure(t);%     clf(fig(t));
        
        posX = 1700 + 450 * ( nShifts/2 -1 - mod( nShifts/2-t, nShifts/2 ) );
        posY = 75   + 520 * ( 1 - round(t/nShifts-0.01));
        
        errWidth = 16.64/2;
        posPaper = [0 0 errWidth errWidth];

        set(fig(t), 'Position', [posX posY 400 400], 'Color', [1,1,1], ...
                    'PaperUnits', 'centimeters', 'PaperPosition', posPaper, ...
                    'PaperSize', [posPaper(3) posPaper(4)] );
        
        a(1) = axes('Position', [0 0 1 1]);

        hold on
            [Y X] = meshgrid(x,x);
            [tmp c(t)] = contourf(X,Y,f, 50);
            set(c(t), 'LineStyle', 'none');
            colormap bone;
            box off;
            set(gca, 'XTick',[], 'YTick', []);

%             marg = 0.2;
%             rect = rectangle('Position', [marg, marg, 1-2*marg, 1-2*marg],...
%                       'LineWidth', 2, 'LineStyle', '--');
        hold off
        
        %-- Saving
        nameFile = sprintf('../figures/shift_%d', t);
        saveas(fig(t), nameFile, 'pdf');
    end
end

%----------------------- Database formating and writing -----------------------%
if dbase
    
    for t = 1:nShifts
        nPatterns = nPeaks;
        nFeatures = nDims + 2;
        nCls      = nPeaks;

        dataOut  = zeros(nPatterns, nFeatures);
        tagsOut  = zeros(nPatterns, 1);

        for p = 1:nPeaks
            dataOut(p,:) = [posPeaks(p,:,t), height(p,t), width(p,t)];
            tagsOut(p,1) = 1;
        end

        nameFile = sprintf('../database/%s/%dshifts/shift_%d.db', ...
                           nameDb, nShifts, t);
        util_writeDb(nPatterns,nFeatures,nCls, dataOut, tagsOut, nameFile);
    end

end
        
fprintf('/*------- End - moving peaks benchmark -------*/\n');