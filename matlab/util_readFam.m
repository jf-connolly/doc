function  artmap = util_readFam(nameFile, nBlocks, nReplications, ...
                                sizeEnsemble, currRep, currBlk, currMbr)

% fprintf('/*-- block %d, rep %d/%d, member %d/%d', currBlk, currRep, ...
%         nReplications, currMbr, sizeEnsemble(currBlk, currRep));

f = fopen(nameFile,'r');

    for r = 1:nReplications
        for t = 1:nBlocks

            savedRep = fread(f, 1, 'int32')+1;
            savedBlk = fread(f, 1, 'int32')+1;
            fseek(f, -8, 0);

            if nBlocks > 1
                sizeEns = sizeEnsemble(savedBlk, savedRep);
            else
                sizeEns  = sizeEnsemble(savedRep);
                savedBlk = 1;
            end
            for e = 1:sizeEns
                %-- General informations
                savedRep  = fread(f, 1, 'int32')+1;
                savedBlk  = fread(f, 1, 'int32')+1;
                nClasses  = fread(f, 1, 'int32');
                nFeatures = fread(f, 1, 'int32');
                sizeNn    = fread(f, 1, 'int32');

                %-- Adjust for batch learning
                if nBlocks == 1,   savedBlk = 1;   end

                %-- Pass over all the other neural networks ...
                if r == currRep && savedBlk == currBlk && e == currMbr
                    %-- Matrixes initialization (with the general informations)
                    artmap = struct( ...
                     'nClasses',   nClasses,                   ...
                     'nFeatures',  nFeatures,                  ...
                     'sizeNn',     sizeNn,                     ...
                     'sizeClasses',zeros(1,nClasses),          ...
                     'h',          zeros(1,4),                 ...
                     'W',          zeros(sizeNn, nFeatures*2), ...
                     'normW',      zeros(sizeNn, 1),           ...
                     'normWalpha', zeros(sizeNn, 1),           ...
                     'Wab',        zeros(sizeNn, nClasses)     ...
                    );

                    %-- Class sizes
                    artmap.sizeClasses(1,:) = fread(f, nClasses, 'int32');

                    %-- Hyperparameters
                    artmap.h(:) = fread(f, 4, 'single');

                    %-- Wij, node F2 by node F2
                    for j = 1:artmap.sizeNn
                        artmap.W(j,:) = fread(f, 2*artmap.nFeatures, 'single');
                    end

                    %-- |W| for each node
                    artmap.normW(:) = fread(f, artmap.sizeNn, 'single');

                    %-- 1 / (|W| + alpha) for each node
                    artmap.normWalpha(:) = fread(f, artmap.sizeNn, 'single');

                    %-- Wab, node F2 by node F2
                    for j = 1:artmap.sizeNn
                        artmap.Wab(j,:) = fread(f, nClasses,'int32');
                    end

                    fclose(f);
                    return;
                else
                    %-- Offset computed on the basis of the genral informations
                    sizeData = 4;
                    sizeF1   = nFeatures*2;
                    sizeF2   = sizeNn;
                    offset   = sizeData *( 4        ... %- Hyperparameters
                                + nClasses          ... %- Class sizes
                                + sizeF2*(sizeF1+2) ... %- W, normW, & normAlpha
                                + sizeF2*nClasses   );  %- Wab

                    fseek(f, offset, 0);
                end
                
                %-- Next member
                e = e+1;
            end
        end
    end
    
fprintf('!!! Returned nothing !!!\n');
fclose(f);
