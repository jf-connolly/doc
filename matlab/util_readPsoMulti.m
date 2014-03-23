function  [pso] = util_readPsoMulti(nameFile, nBlocks, nRep)

f = fopen(nameFile,'r');

    %-- General informations
    temp      = fread(f, 1, 'int32');
    sizeSwarm = fread(f, 1, 'int32');
    nDim      = fread(f, 1, 'int32');
    nObj      = fread(f, 1, 'int32');
    nItMax    = 10;

    %-- Matrixes initialization
    pso = struct( ...
     'nBlocks',       nBlocks, ...
     'nRep',          nRep, ...
     'cReplication' , zeros(                 nItMax, nBlocks, nRep),...
     'nIterations',   zeros(                         nBlocks, nRep),...
     'sizeSwarm',     sizeSwarm, ...
     'nDimensions',   nDim, ...
     'nObjectives',   nObj, ...
     'gbest',         zeros(nItMax,                  nBlocks, nRep, nObj),...
     'lbest',         zeros(sizeSwarm,       nItMax, nBlocks, nRep, nObj),...
     's',             zeros(sizeSwarm, nDim, nItMax, nBlocks, nRep),...
     'sPrev',         zeros(sizeSwarm, nDim, nItMax, nBlocks, nRep),...
     'p',             zeros(sizeSwarm, nDim, nItMax, nBlocks, nRep, nObj),...
     'sPm',           zeros(sizeSwarm,       nItMax, nBlocks, nRep),      ...
     'sSz',           zeros(sizeSwarm,       nItMax, nBlocks, nRep),      ...
     'sCn',           zeros(sizeSwarm,       nItMax, nBlocks, nRep),      ...
     'pPm',           zeros(sizeSwarm,       nItMax, nBlocks, nRep, nObj),...
     'pSz',           zeros(sizeSwarm,       nItMax, nBlocks, nRep, nObj), ...
     'pCn',           zeros(sizeSwarm,       nItMax, nBlocks, nRep, nObj) ...
    );

    gbestTemp = zeros(1,nObj);
    
    %-- Reset and read
    fseek(f, 0, -1);
    flag = 1;
    
    for r = 1:nRep
        for t = 1:nBlocks
[r t];
            i=1;
            while flag
                pos = ftell(f);

                %-- General informations
                cRepTemp  = fread(f, 1, 'int32');
                sizeTemp  = fread(f, 1, 'int32');
                dimTemp   = fread(f, 1, 'int32');
                objTemp   = fread(f, 1, 'int32');
                nItTemp   = fread(f, 1, 'int32');

                gbestTemp = fread(f, nObj, 'int32')+1;
                
                if pso.nIterations(t,r) < nItTemp
                    pso.cReplication(i,t,r) = cRepTemp;
                    pso.sizeSwarm           = sizeTemp;
                    pso.nDimensions         = dimTemp;
                    pso.nObjectives         = objTemp;
                    pso.nIterations(t,r)    = nItTemp;
                    pso.gbest(i,t,r,:)      = gbestTemp;
                else
                    fseek(f, pos, -1);
                    break
                end

                %-- Position, previous position (for velocity), & best position
                for n = 1:pso.sizeSwarm
                    
                    pso.s    (n,:,i,t,r) = fread(f, nDim, 'single');
                    pso.sPrev(n,:,i,t,r) = fread(f, nDim, 'single');
                    
                    for o = 1:nObj
                        pso.p(n,:,i,t,r,o) = fread(f, nDim, 'single');
                    end
                    
                    %-- Local Best (used to identify the subswarm)
                    for o = 1:nObj
                        pso.lbest(n,i,t,r,o) = fread(f, 1, 'int32') +1;
                    end

                    %-- fitness (current and best position)
                    pso.sPm(n,i,t,r) = fread(f, 1, 'single');
                    pso.sSz(n,i,t,r) = fread(f, 1, 'single');
                    pso.sCn(n,i,t,r) = fread(f, 1, 'single');
                    for o = 1:nObj
                        pso.pPm(n,i,t,r,o) = fread(f, 1, 'single');
                        pso.pSz(n,i,t,r,o) = fread(f, 1, 'single');
                        pso.pCn(n,i,t,r,o) = fread(f, 1, 'single');
                    end
                end
                i = i+1;
            end
        end
    end

fclose(f);
