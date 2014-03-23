function  [archive] = util_readArchive(nameFile,nBlocks,nReplications,nIterations)

nItMax = max(max(nIterations));

f = fopen(nameFile,'r');

    %-- General informations
    temp   = fread(f, 1, 'int32');
    asize  = fread(f, 1, 'int32');
    temp   = fread(f, 1, 'int32');
    nMems  = fread(f, 1, 'int32');
    nDims  = fread(f, 1, 'int32');
    
    %-- Matrixes initialization
    archive = struct( ...
     'nBlocks',       nBlocks,                                            ...
     'nReplications', nReplications,                                      ...
     'nIterations',   nIterations,                                        ...
     'cReplication' , zeros(             nItMax, nBlocks, nReplications), ...
     'size',          asize,                                              ...
     'nMemetics',     nMems,                                              ...
     'nDimensions',   nDims,                                              ...
     'nFilled',       zeros(             nItMax, nBlocks, nReplications), ...
     'grid',          zeros(nMems, 2,    nItMax, nBlocks, nReplications), ...
     'coords',        zeros(2, asize,     nItMax, nBlocks, nReplications), ...
     's',             zeros(asize, nDims, nItMax, nBlocks, nReplications), ...
     'sPm',           zeros(asize,        nItMax, nBlocks, nReplications), ...
     'sSz',           zeros(asize,        nItMax, nBlocks, nReplications)  ...
    );

    %-- Reset and read
    fseek(f, 0, -1);
    
    for r = 1:nReplications
        for t = 1:nBlocks
            for i=1:nIterations(t,r);
                
                %-- General informations
                archive.cReplication(i,t,r) = fread(f, 1, 'int32');
                archive.size                = fread(f, 1, 'int32');
                archive.nFilled(i,t,r)      = fread(f, 1, 'int32');
                archive.nMemetics           = fread(f, 1, 'int32');
                archive.nDimensions         = fread(f, 1, 'int32');

                %-- Mopso grid
                archive.grid(1:nMems,1,i,t,r) = fread(f, nMems, 'single');
                archive.grid(1:nMems,2,i,t,r) = fread(f, nMems, 'single');
                
                for n = 1:archive.nFilled(i,t,r)
                    archive.coords(:,n,i,t,r) = fread(f, 2, 'int32');
                end
                
                %-- Position and fitnesses
                for n = 1:archive.nFilled(i,t,r)
                    archive.s(n,:,i,t,r) = fread(f, nDims, 'single');
                    archive.sPm(n,i,t,r) = fread(f, 1,     'single');
                    archive.sSz(n,i,t,r) = fread(f, 1,     'single');
                end
            end
        end
    end

fclose(f);
