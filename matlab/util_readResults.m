function  [result] = util_readResults(nameFile, nBlocks, nReps)
f = fopen(nameFile,'r');

    %-- Read the info number of classes for intialization and rewind
    temp   = fread(f, 1, 'int32');
    sizeDB = fread(f, 1, 'int32');
    sizeDB = 15;
    nCls   = fread(f, 1, 'int32');
    nPats  = fread(f, 1, 'int32');
    nPats  = fread(f, 1, 'int32');
    sizeSw = fread(f, 1, 'int32');
    
    fseek(f,0,-1);

    result = struct( ...
                'nBlocks',          nBlocks,                                 ...
                'nReplications',    nReps,                                   ...
                'nClasses',         nCls,                                    ...
                'nPatternsTest',    zeros(                nBlocks, nReps),   ...
                'nPatternsLearned', zeros(                nBlocks, nReps),   ...
                'sizeSwarm',        sizeSw,                                  ...
                'sizeEns',          zeros(                nBlocks, nReps),   ...
                'members',          zeros(sizeSw,         nBlocks, nReps)-1, ...
                'nClassOk',         zeros(                nBlocks, nReps),   ...
                'nClassWg',         zeros(                nBlocks, nReps),   ...
                'normCpn',          zeros(                nBlocks, nReps),   ...
                'predictions',      zeros(sizeDB, sizeSw, nBlocks, nReps)-1, ...
                'predisemble',      zeros(sizeDB,         nBlocks, nReps)-1, ...
                'trueClasses',      zeros(sizeDB,         nBlocks, nReps)-1, ...
                'sizeClasses',      zeros(nCls,   sizeSw, nBlocks, nReps)-1, ...
                'membersTime',      zeros(        sizeSw, nBlocks, nReps)-1, ...
                'cRep',             zeros(                nBlocks, nReps)    ...
             );
	flagRestart = 1;
	r = 1;
    while r <= nReps
		t = 1;
        while t <= nBlocks

%if r == 13 && t == 16 && flagRestart
%	r = 13;
%	t = 1;
%	flagRestart = 0;
%end
% fprintf('r: %d/%d, t: %d\n', r, nReps, t);

            %-- General informations
            result.cRep(t,r)             = fread(f, 1, 'int32');
            sizeDb                       = fread(f, 1, 'int32');
            result.nClasses              = fread(f, 1, 'int32');
            result.nPatternsTest(t,r)    = fread(f, 1, 'int32');
            result.nPatternsLearned(t,r) = fread(f, 1, 'int32');
            result.sizeSwarm             = fread(f, 1, 'int32');

            %-- Ensemble info
            result.sizeEns(t,r)               = fread(f, 1,       'int32' );
            sizeEns                           = result.sizeEns(t,r);
			if sizeEns
	            result.members(1:sizeEns,t,r) = fread(f, sizeEns, 'int32' ) +1;
            end
            
 %           fprintf('r: %d, t: %d\n  cRep: %d, Size db: %d, nClasses: %d\n', ... 
%		    		r, t, result.cRep(t,r), sizeDb, result.nClasses);
%            fprintf('  nPatTest: %d, nPatLearn: %d, szSwarm %d, sizeEns: %d\n',...
%					result.nPatternsTest(t,r), result.nPatternsLearned(t,r), ...
%					result.sizeSwarm, result.sizeEns(t,r));

            %-- Individual classification results
            result.nClassOk(t,r) = fread(f,1, 'single');
            result.nClassWg(t,r) = fread(f,1, 'single');
            result.clsRate (t,r) = fread(f,1, 'single');
            result.normCpn (t,r) = fread(f,1, 'single');

%			fprintf('  nClassOk: %1.2f, nClassWg : %1.2f, clsRate: %1.2f, normCpn, %1.2f\n',...
%			         result.nClassOk(t,r), result.nClassWg(t,r), ...
%					 result.clsRate (t,r), result.normCpn (t,r));

            %-- General classification results & true class labels
            nPats = result.nPatternsTest(t,r);
			if sizeEns
	            for n = 1:sizeEns
    	            result.predictions(1:nPats,n,t,r) = fread(f,nPats,'int32');
        	    end
            end

			if nPats
	            result.predisemble( 1:nPats,t,r )     = fread(f,nPats,'int32');
    	        result.trueClasses( 1:nPats,t,r )     = fread(f,nPats,'int32');
            end

            %-- Class sizes
%	    fprintf('  Size network:');
            for n = 1:sizeEns
                result.sizeClasses(:,n,t,r) = fread(f, nCls, 'int32');
%		fprintf(' %d', sum(result.sizeClasses(:,n,t,r)));
            end
 %           fprintf('\n');

            %-- Convergenc time
  %          fprintf('  Conv. time:');
	    if sizeEns
            	result.membersTime(1:sizeEns,t,r) = fread(f, sizeEns, 'int32');
 		for m=1:sizeEns
%			fprintf(' %d', result.membersTime(m,t,r));
		end
            end
 %           fprintf('\n');
t = t+1;
    	end
r = r+1;
    end
fclose(f);
