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
widthMc = 0;  
nameLn  = 'Inc';
nameAl  = 'dnpso';

nRep      = 50;
time      = 12;
iteration = -1;
pbest     = 1;
radius    = 0.1;
numbers   = 0;

%-- What to do
loadFiles  = 1;
saveResult = 0;
loadResult = 0;
%-- What to plot
graphObj   = 0;

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
    pso = util_readMopso(nameFile, nBlocks, nRep);
    szSwarm = pso.sizeSwarm;
    
    %-- Archive
    nameFile = sprintf('../savedStuff/%s_%s%sHdnc%s_%dBlocks_wd%d.arc', ...
                       nameDb, nameSc, nameLn, nameAl, nBlocksName, widthMc);
    archive  = util_readArchiveMopso(nameFile, nBlocks, nRep, pso.nIterations);
    szArchive = archive.size;

%     %-- Result
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
    
        
    %-- Trimming down to save less
    reducedPso = pso;   reducedArc = archive;

    reducedPso.sPrev = [];
    reducedPso.cReplication = [];
    reducedArc.cReplication = [];

    for sr = nRep:-1:1
        if sr ~= cRep
            reducedPso.nIterations(:,sr) =[];
            reducedPso.s(:,:,:,:,sr)     =[];  reducedPso.p(:,:,:,:,sr,:)   =[];
            reducedPso.gbest(:,:,sr,:)   =[];  
            reducedPso.sPm(:,:,:,sr)     =[];  reducedPso.pPm(:,:,:,sr,:)   =[];
            reducedPso.sSz(:,:,:,sr)     =[];
            reducedPso.pSz(:,:,:,sr,:)   =[];

            reducedArc.nIterations(:,sr) =[]; reducedArc.nFilled(:,sr)      =[];
            reducedArc.grid(:,:,:,:,sr)  =[]; reducedArc.coords(:,:,:,:,sr) =[];
            reducedArc.s(:,:,:,:,sr)     =[];
            reducedArc.sPm(:,:,:,sr)     =[]; reducedArc.sSz(:,:,:,sr)      =[];
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
if graphObj, fprintf('/*-- Plot objectives space\n');
    %--------------------------- Graphic parameters ---------------------------%
	%-- General
	graphWidth = 8;
	ratio      = 1/1;

	%-- Marker & font
	mkrSubswarms = 'os^dv><hpx+os^dv><hpx+';  mkrFree = '*';
	sizeMkr = 7;   sizeFont = 10;   widthLine = 0.5;
	clrMin  = 0;   clrMax   = 1;    pow       = 1;

	%-- Figure initialization
	close all;
    width    = res*graphWidth+2;   height = res*graphWidth/ratio+2;
    posPaper = [0 0 graphWidth graphWidth/ratio];

    %---------------------------------- Graph ---------------------------------%
    t = time;

    o       = nObjectives +1;
    posScrn = [1300 sizeScrn(2)-104-50*(o-1) width height];
    fig(o)  = figure(o);   clf(fig(o));   set(fig(o),'Color',[1,1,1]);
    set(fig(o), 'Position', posScrn, 'PaperUnits', 'centimeters', ...
                'PaperSize', [posPaper(3) posPaper(4)],...
                'PaperPosition', posPaper, 'Color', [1,1,1]);
    ax = axes('Position', [0 0 1 1]);

	hold on
		%---------------------- Explored objective space ----------------------%
		nIt = reducedPso.nIterations(t, 1);
		
        sz = szSwarm*nIt;
		f = [ [ reshape((1-reducedPso.sPm(:, 1:nIt, t, r))*100, sz, 1), ...
                reshape((  reducedPso.sSz(:, 1:nIt, t, r)),     sz, 1) ] ];
		
		for n = size(f,1):-1:1
            if f(n,2) > 9000,   f(n,:) = [];   end
		end
		
        yMin = min(f(:,1))-0.01;
        yMax = max(f(:,1))+0.01;

        set(gca, 'Position', [0.15 0.13 0.82 0.85]);%, ...
%                  'XLim', [0 1], 'YLim', [yMin-0.01 1.01]);

		clr = [0.7179    0.8194    0.8194];
	%     clr = [0.2222    0.2222    0.3108];
% 		plot(f(:,1), f(:,2), 'ko', 'Color', clr, 'markerfacecolor', clr, ...
% 			'MarkerSize',sizeMkr);
		
        
        %--------------------------- Explored space ---------------------------%
        clr = bone(nIt+5);

        for i = 1:nIt
            for n = 1:szSwarm
                z = [ pso.sSz(:,i,t,r) (1-pso.sPm(:,i,t,r))*100 ];
            end

            for n = size(z,1):-1:1
               if z(n,1) > 9000
                   z(n,:) = [];
               end
            end
            plot(z(:,1), z(:,2), 'o', 'Color', clr(i,:), ...
                      'markerfacecolor', clr(i+5,:), 'MarkerSize',sizeMkr-4,...
                      'LineWidth', widthLine);
        end

% 		%----------------------------- Boundaries -----------------------------%
% 		x = [0, reducedArc.boundaries'];
% 		for m = 1:size(x,2)
% 			plot([x(m) x(m)], [0 100], 'k--','LineWidth', widthLine);
%         end
        
		%-- Particles
%		pl = plot(z(:,1), z(:,2), 'o', 'Color', [0 0 0], ...
%			  'markerfacecolor', [1 1 1], 'MarkerSize',sizeMkr,...
%			  'LineWidth', widthLine);

		%------------------------------- Archive ------------------------------%
        %-- Archive
		clr = [0.2222    0.2222    0.3108];
        nFilled = reducedArc.nFilled(nIt,t);
		zArc = [    reducedArc.sSz(1:nFilled, nIt, t) ...
                 (1-reducedArc.sPm(1:nFilled, nIt, t))*100 ];
        pl = plot(zArc(:,1), zArc(:,2), 's', 'Color', [0 0 0], ...
          'markerfacecolor', [1 1 1], 'MarkerSize',sizeMkr,...
          'LineWidth', widthLine);
		
	hold off

    %---------------------------- Labels and stuff ----------------------------%
    xlabel('Network size ($F_2$ layer)', 'FontName', 'Times New Roman',...
     'FontSize', sizeFont, 'Interpreter','latex');
    ylabel('Error rate (\%)', 'FontName', 'Times New Roman',...
     'FontSize', sizeFont, 'Interpreter','latex');

    maxY = max(f(:,1)) + 1;
    maxX = max(f(:,2)) + 1;
    
%     set(gca, 'YLim', [0 maxY], 'XLim', [0 maxX]);
    set(gca, 'Position', [0.16 0.13 0.81 0.85], ...
             'YLim', [0 maxY], 'XLim', [0 maxX]);

	%-- Saving the figure
	if pbest,    namePb  = 'p'  ;   else namePb  = ''   ;   end
	if o == 1,   nameObj = 'Clr';   else nameObj = 'Cpn';   end

	nameFile = sprintf('../figures/%s%sPar_%s%s_t%d_wd%d', ...
                       nameDb, nameSc, nameAl, namePb, t, widthMc);
    nameFig  = sprintf('-f3');
    print(nameFig, '-dpdf', '-r300', nameFile);
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
end
fprintf('/*----------------------- End ----------------------*/\n');
