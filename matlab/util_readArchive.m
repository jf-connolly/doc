function  [archive] = util_readArchive(nameFile,nBlocks,nReplications,nIterations)

nItMax = max(max(nIterations));

f = fopen(nameFile,'r');

    %-- General informations
    temp         = fread(f, 1, 'int32');
    asize        = fread(f, 1, 'int32');
    sizeMemetic  = fread(f, 1, 'int32');
    widthMemetic = fread(f, 1, 'int32');
    temp         = fread(f, 1, 'int32');
    nDim         = fread(f, 1, 'int32');
    
    %-- Matrixes initialization
    archive = struct( ...
     'nBlocks',       nBlocks,                                           ...
     'nReplications', nReplications,                                     ...
     'nIterations',   nIterations,                                       ...
     'cReplication' , zeros(            nItMax, nBlocks, nReplications), ...
     'size',          asize,                                              ...
     'sizeMemetic',   sizeMemetic,                                       ...
     'widthMemetic',  widthMemetic,                                      ...
     'nMemetics',     zeros(            nItMax, nBlocks, nReplications), ...
     'nDimensions',   nDim,                                              ...
     'boundaries',    zeros(asize, 1),                                    ...
     'nFilled',       zeros(            nItMax, nBlocks, nReplications), ...
     'filled',        zeros(asize,       nItMax, nBlocks, nReplications), ...
     'nMembers',      zeros(asize,       nItMax, nBlocks, nReplications), ...
     's',             zeros(asize, nDim, nItMax, nBlocks, nReplications), ...
     'sPm',           zeros(asize,       nItMax, nBlocks, nReplications), ...
     'sSz',           zeros(asize,       nItMax, nBlocks, nReplications)  ...
    );

    %-- Reset and read
    fseek(f, 0, -1);
    
    for r = 1:nReplications
        for t = 1:nBlocks
            for i=1:nIterations(t,r);
                %-- General informations
                archive.cReplication(i,t,r) = fread(f, 1, 'int32');
                asize                       = fread(f, 1, 'int32');
                sizeMemetic                 = fread(f, 1, 'int32');
                widthMemetic                = fread(f, 1, 'int32');
                nMem                        = fread(f, 1, 'int32');
                archive.nMemetics(i,t,r)    = nMem;
                nDim                        = fread(f, 1, 'int32');

                %-- Memetic boundaries
                archive.boundaries             = fread(f, nMem, 'int32');
                
                %-- How the archive is filled
                archive.nFilled(i,t,r)         = fread(f, 1,    'int32');
                archive.filled(:,i,t,r)        = fread(f, asize, 'int32');
                archive.nMembers(1:nMem,i,t,r) = fread(f, nMem, 'int32');

                %-- Position and fitnesses
                for m = 1:asize
                    if archive.filled(m,i,t,r)
                        archive.s(m,:,i,t,r) = fread(f, nDim, 'single');
                        archive.sPm(m,i,t,r) = fread(f, 1,    'single');
                        archive.sSz(m,i,t,r) = fread(f, 1,    'single');
                    end
                end
                
            end
        end
    end

fclose(f);
