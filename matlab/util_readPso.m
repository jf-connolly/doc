function  [pso] = util_readPso(nameFile, nBlocks, nReplications)

f = fopen(nameFile,'r');

    %-- General informations
    temp      = fread(f, 1, 'int32');
    sizeSwarm = fread(f, 1, 'int32');
    nDim      = fread(f, 1, 'int32');
    nItMax    = 10;

    %-- Matrixes initialization
    pso = struct( ...
     'nBlocks',       nBlocks, ...
     'nReplications', nReplications, ...
     'cReplication' , zeros(                 nItMax, nBlocks, nReplications),...
     'nIterations',   zeros(                         nBlocks, nReplications),...
     'sizeSwarm',     sizeSwarm, ...
     'nDimensions',   nDim, ...
     'gbest',         zeros(nItMax,                  nBlocks, nReplications),...
     'lbest',         zeros(sizeSwarm,       nItMax, nBlocks, nReplications),...
     's',             zeros(sizeSwarm, nDim, nItMax, nBlocks, nReplications),...
     'sPrev',         zeros(sizeSwarm, nDim, nItMax, nBlocks, nReplications),...
     'p',             zeros(sizeSwarm, nDim, nItMax, nBlocks, nReplications),...
     'sPm',           zeros(sizeSwarm,       nItMax, nBlocks, nReplications),...
     'sSz',           zeros(sizeSwarm,       nItMax, nBlocks, nReplications),...
     'pPm',           zeros(sizeSwarm,       nItMax, nBlocks, nReplications),...
     'pSz',           zeros(sizeSwarm,       nItMax, nBlocks, nReplications) ...
    );

    %-- Reset and read
    fseek(f, 0, -1);
    flag = 1;
    
    for r = 1:pso.nReplications
        for t = 1:pso.nBlocks

            i=1;
            while flag
                pos = ftell(f);

                %-- General informations
                cRepTemp  = fread(f, 1, 'int32');
                sizeTemp  = fread(f, 1, 'int32');
                dimTemp   = fread(f, 1, 'int32');
                gbestTemp = fread(f, 1, 'int32')+1;
                nItTemp   = fread(f, 1, 'int32');

                if pso.nIterations(t,r) < nItTemp
                    pso.cReplication(i,t,r) = cRepTemp;
                    pso.sizeSwarm           = sizeTemp;
                    pso.nDimensions         = dimTemp;
                    pso.nIterations(t,r)    = nItTemp;
                    pso.gbest(i,t,r)        = gbestTemp;
                else
                    fseek(f, pos, -1);
                    break
                end

                %-- Position, previous position (for velocity), & best position
                for n = 1:pso.sizeSwarm
                    for d = 1:pso.nDimensions
                        pso.s    (n,d,i,t,r) = fread(f, 1, 'single');
                        pso.sPrev(n,d,i,t,r) = fread(f, 1, 'single');
                        pso.p    (n,d,i,t,r) = fread(f, 1, 'single');
                    end

                    %-- Local Best (used to identify the subswarm)
                    pso.lbest(n,i,t,r) = fread(f, 1, 'int32') +1;

                    %-- fitness (current and best position)
                    pso.sPm(n,i,t,r) = fread(f, 1, 'single');
                    pso.sSz(n,i,t,r) = fread(f, 1, 'single');
                    pso.pPm(n,i,t,r) = fread(f, 1, 'single');
                    pso.pSz(n,i,t,r) = fread(f, 1, 'single');
                end
                i = i+1;
            end
        end
    end

fclose(f);