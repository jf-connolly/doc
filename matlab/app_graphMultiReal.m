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
widthMc = 100;  
nameLn  = 'Inc';
nameAl  = 'multi';
nameHp  = 'Hdnc';
nRep    = 50;

time      = 12;
iteration = -1;
numbers   = 0;

%-- What to do
loadFiles  = 1;
saveResult = 0;
loadResult = 0;
%-- What to plot
graphObj   = 0;

fprintf('\n/*---------------------- Start ---------------------*/\n');
nBlocksName = nBlocks;
if strcmp('Bth', nameLn), nBlocks = 1; end
    
%------------------------------------------------------------------------------%
%-------------------------------- Preprocessing -------------------------------%
nObjectives = 2;
if loadFiles,   fprintf('/*-- Load \n');
    %---------------------------- Reading the data ----------------------------%
    %-- Swarm
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.pso', ...
                       nameDb, nameSc, nameLn, nameHp, nameAl, nBlocksName, widthMc)
    pso = util_readPsoMulti(nameFile, nBlocks, nRep);
    szSwarm = pso.sizeSwarm;
    
    %-- Archive
    nameFile = sprintf('../savedStuff/%s_%s%s%s%s_%dBlocks_wd%d.arc', ...
                       nameDb, nameSc, nameLn, nameHp, nameAl, nBlocksName, widthMc)
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
    
%     cRep = 9
        
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
%------------------------------ Saving & loading ------------------------------%
if saveResult,   fprintf('/*-- Save reduced data\n');

    nameFile = sprintf('../savedStuff/%s_obj%s%s_wd%d.mat', ...
					   nameDb, nameSc, nameAl, widthMc);
    save(nameFile, 'reducedPso', 'reducedArc', 'szSwarm', 'szArchive', '-mat');
end

if loadResult
    nameFile = sprintf('../savedStuff/%s_obj%s%s_wd%d.mat', ...
					   nameDb, nameSc, nameAl, widthMc);
    load(nameFile, 'reducedPso', 'reducedArc', 'szSwarm', 'szArchive', '-mat');
end

%------------------------------------------------------------------------------%
%------------------------------- Objective space ------------------------------%
if graphObj,    fprintf('/*-- Plot objectives space\n');
    %--------------------------- Graphic parameters ---------------------------%
	%-- General
	graphWidth = 16.65/3;
	ratio      = 1/1;

	%-- Marker & font
	mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
	sizeMkr = 7;   sizeFont = 10;   widthLine = 0.5;
	clrMin  = 0;   clrMax   = 1;    pow       = 1;

	%-- Figure initialization
	close all;
    width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    posPaper = [0 0 graphWidth graphWidth/ratio];

    posScrn = [1300 sizeScrn(2)-104-550*(1) width height];
    fig(1)  = figure(1);   clf(fig(1));   set(fig(1),'Color',[1,1,1]);
    set(fig(1), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
    ax = axes('Position', [0 0 1 1]);

    %----------------------------- Objective space ----------------------------%
    t   = time;
    nIt = reducedPso.nIterations(t);
    
	hold on
        %--------------------------- Explored space ---------------------------%
        clr = flipud( bone(nIt) );

        maxX = 0;
        for i = 1:nIt
            for n = 1:szSwarm
                z = [ reducedPso.sSz(:,i,t) (1-reducedPso.sPm(:,i,t))*100 ];
            end

            for n = size(z,1):-1:1
               if z(n,1) == 0
                   z(n,:) = [];
               end
            end
            plot(z(:,1), z(:,2), 'o', 'Color', clr(i,:), ...
                      'markerfacecolor', clr(i,:), 'MarkerSize',sizeMkr-4,...
                      'LineWidth', widthLine);
                  
            maxX = max([z(:,1); maxX]);
        end
        
		%----------------------------- Boundaries -----------------------------%
		x = [0, reducedArc.boundaries'];
		for m = 1:size(x,2)
			plot([x(m) x(m)], [0 100], 'k--','LineWidth', widthLine);
        end
        
		%------------------------------- Archive ------------------------------%
        %-- Archive
		clr = [0.2222    0.2222    0.3108];
		zArc = [ reducedArc.sSz(:,nIt, t), (1-reducedArc.sPm(:,nIt, t))*100 ];
		for n = 1:reducedArc.size
			if reducedArc.filled(n, reducedPso.nIterations(t), t)
                
                pl = plot(zArc(n,1), zArc(n,2), 's', 'Color', [0 0 0], ...
                  'markerfacecolor', [1 1 1], 'MarkerSize',sizeMkr,...
                  'LineWidth', widthLine);
			end
        end
		
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
       
%     maxY = max(f(:,1));
%     maxX = max(f(:,2));
%     maxX = max(max(max(reducedArc.sSz))) + 20;
    
%     set(gca, 'YLim', [0 maxY], 'XLim', [0 maxX]);
    set(gca, 'Position', [0.2 0.162 0.78 0.81], ...
             'YLim', [0 100], 'XLim', [0 474], ...
             'FontName', 'Times New Roman', ...
             'FontSize', 9);

	%-- Saving the figure
	nameFile = sprintf('../figures/%s%sPar_%s_t%d_wd%d', ...
                       nameDb, nameSc, nameAl, t, widthMc);
	saveas(fig(1), nameFile, 'pdf');
else
    %---------------------------------- Graph ---------------------------------%
    t = time;
	%------------------------ Explored objective space ------------------------%
	nIt = reducedPso.nIterations(t, 1);
	
	sz = szSwarm*nIt;
	sz = szSwarm;
	f = [ [ reshape((1-reducedPso.sPm(:, nIt, t))*100, sz, 1), ...
			reshape((  reducedPso.sSz(:, nIt, t)),     sz, 1) ] ];

	for n = size(f,1):-1:1
		if   f(n,2) > 9000 || f(n,2) == 0,   f(n,:) = [];   end
	end

	nSz     = 15;
	nAc     = 15;
	widthSz = max(f(:,2))/nSz;
	widthAc = 100        /nAc;
	graph = zeros(nSz,nAc);

	for n = 1:size(f,1)
		coord = [ceil(f(n,1)/widthAc) ceil((f(n,2)-0.01)/widthSz)];
		graph(coord(1), coord(2)) = graph(coord(1), coord(2))+1;
	end

	flipud(graph)

	%--------------------------------- Archive --------------------------------%
	%-- Archive
	zArc = [ (1-reducedArc.sPm(:, nIt, t))*100 ...
			    reducedArc.sSz(:, nIt, t)      ];
	zArc(zArc(:,2) == 0,:) = []
	
	%------------------------------- Performance ------------------------------%
	
	for r = 1:50
		zArc = [ (1-archive.sPm(:, nIt, t, r))*100 ...
				    archive.sSz(:, nIt, t, r)      ];
		zArc(zArc(:,2) == 0,:) = [];
		dataMeanArc(r,:) = mean(zArc);
		dataStdArc(r,:) = std(zArc);

		zPso = [ (1-pso.pPm(:, nIt, t, r))*100 ...
				    pso.pSz(:, nIt, t, r)      ];
		zPso(zPso(:,2) == 0,:) = [];
		dataMeanPso(r,:) = mean(zPso);
		dataStdPso(r,:) = std(zPso);
	end

	meanMeanArc = mean(dataMeanArc)
	stdMeanArc  = std(dataMeanArc)
	meanStdArc  = mean(dataStdArc)
	stdStdArc   = std(dataStdArc)

	meanMeanPso = mean(dataMeanPso)
	stdMeanPso  = std(dataMeanPso)
	meanStdPso  = mean(dataStdPso)
	stdStdPso   = std(dataStdPso)
end
fprintf('/*----------------------- End ----------------------*/\n');
