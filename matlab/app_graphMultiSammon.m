%------------------------------------------------------------------------------%
%-- mpbGraph: 2D graphic of the multipeak benchmark and PSO swarm.
%-- Works with util_getName
%------------------------------------------------------------------------------%
%-- DEFINE
temp = get(0,'MonitorPosition');
sizeScrn = [-temp(1,1)+1 temp(1,4)];
res = get(0,'ScreenPixelsPerInch')/2.56;
%-- END DEFINE

%------------------------------------------------------------------------------%
%-------------------------- User-defined parameters ---------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
nameDb  = 'cnrc64';
nameSc  = 'upd';
nBlocks = 12;
widthMc = 1500;  
nameLn  = 'Inc';
nameAl  = 'multi';

nRep      = 50;
time      = 1;
iteration = -1;
pbest     = 1;
radius    = 0.1;
numbers   = 0;

%-- What to do
loadFiles  = 1;
sammon     = 0;
saveResult = 1;
loadResult = 0;
%-- What to plot
graphSch   = 0;
graphObj   = 0;

%-- plot parameters
resolution = 70;
nlevels    = 20;
knn        = 4;

swarm     = 1;
ax        = [0 65; 0.3 1; -0.6 0.5; 0 0.9];

fprintf('\n/*---------------------- Start ---------------------*/\n');
nBlocksName = nBlocks;
if strcmp('Bth', nameLn), nBlocks = 1; end
    
%------------------------------------------------------------------------------%
%-------------------------------- Preprocessing -------------------------------%
nObjectives = 2;
if loadFiles,   fprintf('/*-- Load \n');
    %---------------------------- Reading the data ----------------------------%
    %-- Swarm
    nameFile = sprintf('../savedStuff/%s_%s%sHdnc%s_%dBlocks_wd%d.pso', ...
                       nameDb, nameSc, nameLn, nameAl, nBlocksName, widthMc);
    pso = util_readPsoMulti(nameFile, nBlocks, nRep);
    szSwarm = pso.sizeSwarm;
    
    %-- Archive
    nameFile = sprintf('../savedStuff/%s_%s%sHdnc%s_%dBlocks_wd%d.arc', ...
                       nameDb, nameSc, nameLn, nameAl, nBlocksName, widthMc);
    archive  = util_readArchive(nameFile, nBlocks, nRep, pso.nIterations);
    szArchive = archive.size;

    %-- Result
    nameFile = sprintf('../savedStuff/%s_%s%sEoMA+t%s_%dBlocks_wd%d.result', ...
                       nameDb, nameSc, nameLn, nameAl, nBlocksName, widthMc);
    result   = util_readResults(nameFile, nBlocks, nRep);
    
    %---------------------- Reduce to transferable size -----------------------%
    %-- Find the replication closest to performance ion the objective space
    avCrate = mean(result.clsRate(nBlocks, :));
    avNcpn  = mean(result.normCpn(nBlocks, :));
    
    distShort = 9999;
    dist = zeros(1,nRep);
    for r = 1:nRep
        cRate = result.clsRate(nBlocks,r);
        nCpn  = result.normCpn(nBlocks,r);
        
        dist(r) = pdist([cRate, nCpn; avCrate, avNcpn]);
    end
    
    [distShort , cRep] = min(dist);
    
%     cRep = 1
        
    %-- Trimming down to save less
    reducedPso = pso;   reducedArc = archive;

    reducedPso.sPrev = [];
    reducedPso.cReplication = [];
    reducedArc.cReplication = [];

    for sr = nRep:-1:1
        if sr ~= cRep
            reducedPso.nIterations(:,sr) =[];
            reducedPso.s(:,:,:,:,sr)     =[];  reducedPso.p(:,:,:,:,sr,:)   =[];
            reducedPso.gbest(:,:,sr,:)   =[];  reducedPso.lbest(:,:,:,sr,:) =[];
            reducedPso.sPm(:,:,:,sr)     =[];  reducedPso.pPm(:,:,:,sr,:)   =[];
            reducedPso.sSz(:,:,:,sr)     =[];
            reducedPso.pSz(:,:,:,sr,:)   =[];
            reducedPso.sCn(:,:,:,sr)     =[];  reducedPso.pCn(:,:,:,sr,:)   =[];

            reducedArc.nIterations(:,sr)  = [];  reducedArc.s(:,:,:,:,sr) = [];
            reducedArc.nFilled(:,:,sr)    = [];  reducedArc.sPm(:,:,:,sr) = [];
            reducedArc.filled(:,:,:,sr)   = [];  reducedArc.sSz(:,:,:,sr) = [];
            reducedArc.nMembers(:,:,:,sr) = [];  
        end
    end
end

%------------------------------------------------------------------------------%
%------------------------------- Compute sammon -------------------------------%
if sammon,    fprintf('/*-- Sammon - time:');
    
    szSwarm   = reducedPso.sizeSwarm;
    maxIt     = max(reducedPso.nIterations);
    sammonMap = zeros(szSwarm,2, maxIt, nBlocks);
    fitnesses = zeros(szSwarm,2, maxIt, nBlocks);
    
    %---------------- Initial mapping - for the initial matrix ----------------%
    positions       = reducedPso.s(:,:,1,1);
    distancesVector = pdist(positions);
    distances       = squareform(distancesVector);

    maxIt   = max(reducedPso.nIterations);
    
    opts    = statset('Display','off', 'MaxIter', 1000, ...
                   'TolFun', 0.000001, 'TolX', 0.000001);
    tempMap = mdscale( distances, 2, 'Criterion', 'sammon', ...
                         'Start', 'random', 'Options', opts );
                               
    sammonMap(:,:,1,1) = tempMap;
    fitnesses(:,:,1,1) = [1-reducedPso.sPm(:,1,1), reducedPso.sSz(:,1,1)];
    
    previousPos = positions;
    previousSam = tempMap;
    
    %--------------------------- Rest of the mapping --------------------------%
    addr = 0;
    for t = 1:nBlocks
        fprintf(' %d', t);
        nIt = reducedPso.nIterations(t);
        i = 1;
        while i <= nIt
            %-- skip the first iteration of the first block
            if t == 1 && i == 1, i = 2; end

            %-- Position vector and fixe equality for the mdscale function
            currentPos = reducedPso.s(:,:,i,t);
            
            nDims  = reducedPso.nDimensions;
            tested = ones(1, nDims);
            while min(tested)
                [tested index] = max( eq(previousPos, currentPos) );
                
                if min(tested)
                    currentPos(index(1),:) = currentPos(index(1),:) + ...
                                             rand(1, nDims)/2-0.25;
                end
            end
            positions  = [ previousPos; currentPos ];
            
            %-- Initial position matrix in the Sammon space
            distancesVector = pdist(positions);
            distances       = squareform(distancesVector);
            
            matrixInit = [previousSam ; previousSam+ rand(szSwarm,2)-0.5];

            opts = statset('Display','off', 'MaxIter', 2000, ...
                           'TolFun', 0.000001, 'TolX', 0.000001);
                       
            tempMap = mdscale( distances, 2, 'Criterion', 'sammon', ...
                               'Start', matrixInit, 'Options', opts );
            
            %-- Keep later half
            sammonMap(:,:,i,t) = tempMap(szSwarm+1:szSwarm*2,:);
            fitnesses(:,:,i,t) = [ 1-reducedPso.sPm(:,i,t), ...
                                     reducedPso.sSz(:,i,t) ];
                                       
            %-- Reassign previous positions
            previousPos = positions(szSwarm+1:2*szSwarm,:);
            previousSam = tempMap  (szSwarm+1:szSwarm*2,:);
            
            i = i + 1;
        end
    end
    fprintf('\n');
end

%------------------------------------------------------------------------------%
%------------------------------ Saving & loading ------------------------------%
if saveResult,   fprintf('/*-- Save reduced data\n');

    nameFile = sprintf('../savedStuff/%s_obj%s%s_wd%d.mat', ...
					   nameDb, nameSc, nameAl, widthMc);
    save(nameFile, 'reducedPso', 'reducedArc', 'positions', 'fitnesses', ...
	               'sammonMap', 'szSwarm', 'szArchive', '-mat');
end

if loadResult
    nameFile = sprintf('../savedStuff/%s_obj%s%s_wd%d.mat', ...
					   nameDb, nameSc, nameAl, widthMc);
    load(nameFile, 'reducedPso', 'reducedArc', 'positions', 'fitnesses', ...
	               'sammonMap', 'szSwarm', 'szArchive', '-mat');
end

%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%----------------------------- Graphic parameters -----------------------------%
if graphSch || graphObj
	%-- General
	graphWidth = 8;
	ratio      = 1/1;

	%-- Marker & font
	mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
	sizeMkr = 4;   sizeFont = 10;   widthLine = 0.5;
	clrMin  = 0;   clrMax   = 1;    pow       = 1;

	%-- Figure initialization
	close all;
    width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    posPaper = [0 0 graphWidth graphWidth/ratio];
end
%------------------------------------------------------------------------------%
%-------------------------------- Search space --------------------------------%
if graphSch,   fprintf('/*-- Plot search space\n');
   
    t = time;
    
    %------------------------------- Format data ------------------------------%
    %-- All positions for a time t
    nIt    = reducedPso.nIterations(t);
    allMap = zeros(szSwarm*nIt,2);
    allFit = zeros(szSwarm*nIt,2);
    allPos = zeros(szSwarm*nIt,4);

    addr = 0;
    for i = 1:nIt
        allMap(addr+1:addr+szSwarm,:) = sammonMap(:,:,i,t);
        allFit(addr+1:addr+szSwarm,:) = fitnesses(:,:,i,t);
        allPos(addr+1:addr+szSwarm,:) = reducedPso.s(:,:,i,t);
        addr = addr +szSwarm;
    end
    
    %-- Filter data outside search space
    nPos = size(allMap,1);
    for n = nPos:-1:1
        if min(allPos(n,:)) < 0 || max(allPos(n,:)) > 1
            allMap(n,:) = [];
            allFit(n,:) = [];
            allPos(n,:) = [];
        end
    end

%     knn = size(allFit,1);

    %-- Normalization network sizes
    minSz = min(allFit(:,2));
%     maxSz = max(allFit(:,2));
    maxSz = 584;
    allFit(:,2) = allFit(:,2)/maxSz;
    
    %---------------------------- Grid & fitnesses ----------------------------%
    %-- Meshgrid & fitnesses
    mins = round(min(allMap) * 100) / 100;
    minX = mins(1);   minY = mins(2);

    maxs = round(max(allMap) * 100) / 100;
    maxX = maxs(1);   maxY = maxs(2);

    xObj = minX : (maxX-minX)/resolution : maxX+(maxX-minX)/resolution;
    yObj = minY : (maxY-minY)/resolution : maxY;

    [Y X] = meshgrid(yObj, xObj);

    nPos = size(allMap,1);
    f    = zeros(resolution, resolution, 2);
    for i = 1:resolution+2
        for j = 1:resolution+1
            testedpos      = repmat([xObj(i) yObj(j)] , nPos, 1);
            distances      = sum((testedpos-allMap).^2,2).^0.5;

            %-- The kNN appraoch
            [sorted index] = sort(distances);
            totDist = sum(distances(index(1:knn)));
            weights = flipud( distances(index(1:knn)) / totDist );
            
            f(i,j,:) = sum(allFit(index(1:knn),:) .* repmat(weights,1,2),1);
            
%             ftemp = f(i,j,:);
%             
%             [mind nearest] = min(distances);
%             f(i,j,:)       = allFit(nearest,:);
        end
        
        if i == resolution+2
            f(i, 1, :) = [0, 0];
            f(i, 2, :) = [1, 1];
        end

    end
    
    %------------------------------ Plot the data -----------------------------%
    for o = 1:nObjectives

        %-- Introduction
		posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
        fig(o)  = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
		set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
					'PaperSize', [posPaper(3) posPaper(4)],...
					'PaperPosition', posPaper, 'Color', [1,1,1]);
		ax = axes('Position', [0.01 0.01 0.98 0.98]);
        cla %- Clear the current axis
        
        hold on

            %--------------------- Explored space & swarm ---------------------%
            %-- Landscape
            clr = [0.1    0.1    0.25];
            
            [tmp ct] = contourf(X, Y, f(:,:,o), nlevels);
            set(ct, 'LineStyle', '-', 'LineWidth', 0.2, 'LineColor', clr);
            colormap bone;
            
            %-- Particles
            cmap    = colormap(bone(101));
            clr     = [0.7179    0.8194    0.8194];
            sizePos = size(allMap,1);

            for n = 1:sizePos
                clr     = cmap( ceil(allFit(n,o)^pow *100) +1, : );
                clrEdge = [0.2100    0.2100    0.2945];
                plot(allMap(n,1), allMap(n,2), 'o', ...
                     'Color', clrEdge, 'MarkerFaceColor', clr, ...
                     'MarkerSize', sizeMkr-1, 'LineWidth', 0.5);
            end

%             %-- Local best other than the gbest
%             gb = reducedPso.gbest(nIt,t,1,o);
%             for j =1:nSS(o)
%                 lb = lbests(j,o);
%                 if lb ~= gb
%                     testedPos = repmat(reducedPso.p(lb, :, nIt, t, 1, o), ...
%                                        sizePos,1);
%                     [val index] = max( eq(testedPos, positions) );
%                 
%                     currentP = sammonMap(index(1),:);
%                     if ~index, fprintf('/*-- Error plotting\n'); end
% 
%                     if subswarms(lb,o),  marker = mkrSubswarms(subswarms(lb,o));
%                     else                 marker = mkrFree;
%                     end
%                     
%                     plot(currentP(1), currentP(2), marker, ...
%                          'Color', [0 0.7 0], 'MarkerFaceColor', [0 1 0], ...
%                          'MarkerSize', sizeMkr, 'LineWidth', 0.5);
%                 end
%             end

            %----------------------------- Archive ----------------------------%
            clr     = [1 1 0];
            clrEdge = [0.6 0.6 0];
            for n = 1:szArchive
                if reducedArc.filled(n, nIt, t)
                    
                    %-- Find the positions by scanning
                    testedPos = repmat(reducedArc.s(n, :, nIt, t), szSwarm, 1);
                    tp = 1;
                    while tp <= t
                        for i = 1:nIt
                            currentPos     = reducedPso.s(:,:,i,tp);
                            [tested index] = max( eq(testedPos, currentPos) );
                            
                            if min(tested)
                               arcS = sammonMap(index(1),:,i,tp);
                               tp   = nBlocks+1;
                               break;
                            end
                        end
                        tp = tp +1;
                    end

                    plot(arcS(1), arcS(2), 's', 'LineWidth', 0.5, ...
                         'Color', clrEdge, 'MarkerFaceColor', clr, ...
                         'MarkerSize', sizeMkr+1);
                end
            end

            %------------------------------ Gbest -----------------------------%
            %-- Filter for outside search space
            for n = 1:szSwarm
                if min(reducedPso.s(n,:,nIt,t)) < 0 ||...
                   max(reducedPso.s(n,:,nIt,t) > 1)
                    fitnesses(n,:,nIt,t) = [1 9999];
                end
            end

            [minGb indexGb] = min(fitnesses(:, o, nIt, t));
            gbest           = sammonMap(indexGb,:,nIt,t);

%             if o == 1
%                 [minGb indexGb] = min(fitnesses(:, 1, nIt, t));
%                 gbest = sammonMap(indexGb,:,nIt,t);
%             else
%                 minSz  = min(fitnesses(:, 2, nIt, t));
%                 sizeFt = size(fitnesses(:, 2, nIt, t),1);
%                 keep   = eq(fitnesses(:, 2, nIt, t), minSz*ones(sizeFt,1));
%                 tempFt = fitnesses(:, 1, nIt, t) + 1-keep;
%                 [minGb indexGb] = min(fitnesses(:, 1, nIt, t));
%             end
            plot(gbest(1), gbest(2), 'd', 'LineWidth', 1.5, ...
                 'Color', [0.7 0 0], 'MarkerFaceColor', [1 0 0], ...
                 'MarkerSize', sizeMkr+2);

        hold off

        box off
        set(gca,'XTick',[], 'YTick',[], 'XColor', [1 1 1], 'YColor', [1 1 1],...
            'XLim', [minX maxX])
%         title('\textbf{Sammon''s mapping}', 'FontName', 'Times New Roman',...
%        'FontWeight', 'Bold', 'FontSize', sizeFont+2, 'Interpreter','Latex');

        %-- Saving the figure
        if pbest,    namePb  = 'p'  ;   else namePb  = ''   ;   end
        if o == 1,   nameObj = 'Clr';   else nameObj = 'Cpn';   end
        nameFile = sprintf('../figures/%s%s%s_%s%s_t%d_wd%d', ...
                           nameDb, nameSc, nameObj, nameAl, namePb, t, widthMc);
                       
%         saveas(fig(o), nameFile, 'pdf');
        print('-dpng', '-r300', nameFile);


    end
end    
%------------------------------------------------------------------------------%
%------------------------------- Objective space ------------------------------%
if graphObj,   fprintf('/*-- Plot objectives space\n');

    t = time;

    o       = nObjectives +1;
    posScrn = [1300 sizeScrn(2)-104-550*(o-1) width height];
    fig(o)  = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
    ax = axes('Position', [0 0 1 1]);

    it = reducedPso.nIterations(t);
    
	hold on
		%---------------------- Explored objective space ----------------------%
		f = [];   nIt = reducedPso.nIterations(t, 1);
		
		x = (1-reducedPso.sPm(:, 1:nIt, t, 1))*100;    f(:,1) = x(:);
		y =    reducedPso.sSz(:, 1:nIt, t, 1);         f(:,2) = y(:);
		
		for n = size(f,1):-1:1
			if f(n,1) == 0,   f(n,:) = [];   end
		end
		
        yMin = min(f(:,1))-0.01;
        yMax = max(f(:,1))+0.01;

        set(gca, 'Position', [0.15 0.13 0.82 0.85]);%, ...
%                  'XLim', [0 1], 'YLim', [yMin-0.01 1.01]);

		clr = [0.7179    0.8194    0.8194];
	%     clr = [0.2222    0.2222    0.3108];
		plot(f(:,2), f(:,1), 'ko', 'Color', clr, 'markerfacecolor', clr, ...
			'MarkerSize',sizeMkr);
		
		%----------------------------- Boundaries -----------------------------%
		x = [0, reducedArc.boundaries'];
		for m = 1:size(x,2)
			plot([x(m) x(m)], [0 100], 'k--','LineWidth', widthLine);
        end
        
% 		%-- Particles
% 		pl = plot(z(:,2), z(:,1), 'o', 'Color', [0 0 0], ...
% 			  'markerfacecolor', [1 1 1], 'MarkerSize',sizeMkr,...
% 			  'LineWidth', widthLine);

		%------------------------------- Archive ------------------------------%
        %-- Archive
		clr = [0.2222    0.2222    0.3108];
		zArc = [ (1-reducedArc.sPm(:,nIt, t))*100, reducedArc.sSz(:,nIt, t) ];
		for n = 1:reducedArc.size
			if reducedArc.filled(n, reducedPso.nIterations(t), t)
                
%                 if zArc(n,1) > avCrate
%                     pl = plot(zArc(n,2), zArc(n,1), 'p', 'Color', [0 0 0], ...
%                       'markerfacecolor', clr, 'MarkerSize',sizeMkr+2,...
%                       'LineWidth', widthLine);
%                 else
                    pl = plot(zArc(n,2), zArc(n,1), 'o', 'Color', [0 0 0], ...
                      'markerfacecolor', clr, 'MarkerSize',sizeMkr,...
                      'LineWidth', widthLine);
%                 end
			end
        end
		
        %-------------------------------- Gbest -------------------------------%
%         %-- Filter for outside search space
%         for n = 1:szSwarm
%             if min(reducedPso.s(n,:,nIt,t)) < 0 ||...
%                max(reducedPso.s(n,:,nIt,t) > 1)
%                 fitnesses(n,:,nIt,t) = [1 9999];
%             end
%         end
% 
%         [minGb indexGb] = min(fitnesses(:, 1, nIt, t));
%         gbest           = fitnesses(indexGb, :, nIt, t);
%         
%         plot(gbest(2), gbest(1)*100, 'd', 'LineWidth', 1.5, ...
%              'Color', [0.7 0 0], 'MarkerFaceColor', [1 0 0], ...
%              'MarkerSize', sizeMkr+2);
% 
%         minSz  = min(fitnesses(:, 2, nIt, t));
%         sizeFt = size(fitnesses(:, 2, nIt, t),1);
%         keep   = eq(fitnesses(:, 2, nIt, t), minSz*ones(sizeFt,1));
%         tempFt = fitnesses(:, 1, nIt, t) + 1-keep;
%         [minGb indexGb] = min(tempFt);
%         gbest           = fitnesses(indexGb, :, nIt, t);
%         
%         plot(gbest(2), gbest(1)*100, 'd', 'LineWidth', 1.5, ...
%              'Color', [0.7 0 0], 'MarkerFaceColor', [1 0 0], ...
%              'MarkerSize', sizeMkr+2);

        %------------------------------- Numbers ------------------------------%
		if numbers
			for n = 1:szSwarm
				text( z(n,2)+0.01, z(n,1)-0.02, int2str(n), ...
					  'FontName', 'Times New Roman', 'Interpreter', 'latex', ...
					  'FontSize', sizeFont-3, 'FontWeight','Bold' );
			end
			for n = szSwarm+1:szSwarm*2
				np = n-szSwarm;
				text( z(n,2)+0.01, z(n,1)-0.02, int2str(np), ...
					  'FontName', 'Times New Roman', 'Interpreter', 'latex', ...
					  'FontSize', sizeFont-3, 'FontWeight','Bold' );
			end
        end
        
	hold off

    %---------------------------- Labels and stuff ----------------------------%
	xlabel('Network size ($F_2$ layer)', 'FontSize', sizeFont, ...
           'Interpreter','latex');
	ylabel('Error rate (\%)', 'FontSize', sizeFont, 'Interpreter','latex');
       
    maxY = max(f(:,1));
    maxX = max(max(max(reducedArc.sSz))) + 20;
    
%     set(gca, 'YLim', [0 maxY], 'XLim', [0 maxX]);
    set(gca, 'Position', [0.15 0.13 0.81 0.85], ...
             'YLim', [0 maxY], 'XLim', [0 maxX]);

	%-- Saving the figure
	if pbest,    namePb  = 'p'  ;   else namePb  = ''   ;   end
	if o == 1,   nameObj = 'Clr';   else nameObj = 'Cpn';   end

	nameFile = sprintf('../figures/%s%sPar_%s%s_t%d_wd%d', ...
                       nameDb, nameSc, nameAl, namePb, t, widthMc);
	saveas(fig(nObjectives+1), nameFile, 'pdf');
end
fprintf('/*----------------------- End ----------------------*/\n');
