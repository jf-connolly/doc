%------------------------------------------------------------------------------%
%-- Incremental learning scenario - Class update
%------------------------------------------------------------------------------%
fprintf('\n/*------ Generating data for learning ------*/\n');
%---------------------------------- Settings ----------------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
%
%-- General
nBlocks     = 26;
pourcentage = 10;           %-- Percentage of the data for the first bloc
nameDb      = 'mobof32';
nRepExt     = 5;
nRep        = 10;

NEW         = 0;
BATCH       = 1;

%-- External db parameters
ratioEt = 1/6;              %-- Fraction of Db going to E_b
normEt  = 30;               %-- Maximum number of patterns  / class
if BATCH  ratioEt = 0;  end %-- No external data base for batch mode learning

%-- Names
if BATCH     nameLn = 'B';    else  nameLn = 'D';  end
if ratioEt,  nameEx = 'Ext';  else  nameEx = '';   end

%-- Number of classes per block
nameFile = sprintf('../database/%s/learn.db', nameDb);
dbase = util_readHeader(nameFile);

nClsPblock = zeros(1,nBlocks);
for b = 1:nBlocks
    nClsPblock(b) = ceil(dbase.nClasses/nBlocks - (b-1)/nBlocks);
end
nClsPblockSum = [0 cumsum(nClsPblock)];
if sum(nClsPblock) ~= dbase.nClasses,
    fprintf('/*-- Error - Nb. of class per block');
end

%-- Random class order definition for different replications
[order, success] = util_order(NEW, nRepExt, dbase.nClasses, nameDb);
if ~success, break; end;

%------------------------------------------------------------------------------%
%-------------------------------- Learning data -------------------------------%
%-- Reading the db
fprintf('/*--- Class update - learning data sets ----*/\n');
fprintf('/*-- Data reading\n');
nameFile = sprintf('../database/%s/learn.db', nameDb);

dbase = util_readDb(nameFile);
nPatterns = size(dbase.data,1);
nFeatures = size(dbase.data,2);
nCls      = max(dbase.tags) + 1;

%-- Test for the number of blocks
if nBlocks > nCls+1,   fprintf('/*-- To much blocks...\n');   return;   end

%-- Nb of patterns per classes
sizeClsLoad = zeros(1, nCls);
for i = 1:nPatterns
    sizeClsLoad(dbase.tags(i)+1) = sizeClsLoad(dbase.tags(i)+1) + 1;
end
addrCls = [1, (cumsum(sizeClsLoad)+1)];

%----------------------------- Classes separation -----------------------------%
data  = zeros(max(sizeClsLoad), nFeatures, nCls);
tags  = zeros(max(sizeClsLoad), nCls);
sizeL = zeros(1, nCls);

%-- Separation (tags = tags+1)
for k = 1:nCls
    dataL(1:sizeClsLoad(k), :, k) = dbase.data(addrCls(k):addrCls(k+1)-1, :);
    tagsL(1:sizeClsLoad(k), k)    = ones(sizeClsLoad(k),1) * (k);
    sizeL(1,k) = sizeClsLoad(k);
end

%-- E initialization
dataNw = zeros(1, nFeatures, nCls);         tagsNw = zeros(1, nCls);
dataEt = zeros(normEt, nFeatures, nCls);    tagsEt = zeros(normEt, nCls);

%-- Verification if the folder exists
nameFile = sprintf('upd%dBlocks%s', nBlocks, nameEx);
util_CheckAndCreateFolder(nameDb, nameFile);

%------------------------------------------------------------------------------%
%-------------------------- Learning blocs creation ---------------------------%
for rExt = 1:nRepExt
fprintf('/*-- External replication %d/%d - blocks: ', rExt, nRepExt);

%-- Blocks size 
size1st = round(sizeL*pourcentage/100);  
sizeOther = sizeL - size1st;

for c = 1:nCls
    dataOtr(1:sizeOther(c),:,c) = dataL(size1st(c)+1:sizeL(c),:,c);
    tagsOtr(1:sizeOther(c),c)   = tagsL(size1st(c)+1:sizeL(c),c);
end
addrBlk = util_defineFolds( nBlocks, 1, nCls, sizeOther );

%---------------------------- Learning block b (Db) ---------------------------%
%-- For batch learning
if BATCH,	tStart = nBlocks;
else	    tStart = 1;
end

for b = tStart:nBlocks
    fprintf('%d ', b);
    
    %-- Reinitialize if not batch
    if b == 1 || ~BATCH
        sizeDb = zeros(1, nCls);
        dataDb = zeros(1, nFeatures, nCls);   tagsDb = zeros(1, nCls);
    end

    %-- D^b definition - 1rst block ...
    if b == 1 
        sizeDb = size1st;  
        for c = 1:nCls
            dataDb(1:sizeDb(c),:,c,1) = dataL(1:sizeDb(c),:,c);
            tagsDb(1:sizeDb(c),c,1)   = tagsL(1:sizeDb(c),c);
        end
    %-- ... other blocks
    else
        %-- For each classes in D^b
        for cls = 1:nClsPblock(b-1)
            cCls = order( rExt, nClsPblockSum(b-1) + cls);

            sizeCls      = sizeDb(cCls);
            sizeDb(cCls) = sizeDb(cCls) + sizeOther(cCls);
            dataDb(sizeCls+1 : sizeDb(cCls), :, cCls) = ...
                                    dataL(size1st(cCls)+1:sizeL(cCls), :, cCls);
            tagsDb(sizeCls+1 : sizeDb(cCls),    cCls) = ...
                                    tagsL(size1st(cCls)+1:sizeL(cCls),    cCls);
        end
    end
    
    %------------------- Spliting data for E_t and learning -------------------%
    %-- Number of pattern for E_t and for learning
    sizeNw = round(sizeDb*ratioEt);
    sizeLn = sizeDb - sizeNw;
    
    %-- Data for E and for learning
    dataNw = zeros(1, nFeatures, nCls); tagsNw = zeros(1, nCls);
    dataLn = zeros(1, nFeatures, nCls); tagsLn = zeros(1, nCls);
    for k = 1:nCls
        dataNw(1:sizeNw(k), :, k) = dataDb(1:sizeNw(k), :, k);
        tagsNw(1:sizeNw(k), k)    = tagsDb(1:sizeNw(k), k);

        dataLn(1:sizeLn(k),:,k) = dataDb(sizeNw(k)+1:sizeDb(k),:,k);
        tagsLn(1:sizeLn(k),k)   = tagsDb(sizeNw(k)+1:sizeDb(k),k);
    end
    
    %------------------- E_b creation & writing it in a file ------------------%
    %-- Filling / updating Eb
    [dataRd, tagsRd, dataEt, tagsEt] = ...
                          util_externalDbUpdate(dataEt, tagsEt, dataNw, tagsNw);                           
    sizeEt = util_sizeofDb(dataEt);

    %-- Returning the unused data to D_t
    sizeRd = util_sizeofDb(dataRd);
    for k = 1:nCls
        dataLn(sizeLn(k)+1:sizeLn(k)+ sizeRd(k),:,k) = dataRd(1:sizeRd(k),:,k);
        tagsLn(sizeLn(k)+1:sizeLn(k)+ sizeRd(k),k)   = tagsRd(1:sizeRd(k),k);
        sizeLn(k) = sizeLn(k) + sizeRd(k);
    end
        
    %-- Output file - External data base
    dataOut = [];   tagsOut = [];
    for k = 1:nCls
        dataOut = [dataOut ; dataEt(1:sizeEt(k), :, k)];
        tagsOut = [tagsOut ; tagsEt(1:sizeEt(k), k)];
    end
    tagsOut = tagsOut - 1;

    nPatterns = size(dataOut, 1);
    nameFile  = sprintf('../database/%s/upd%dBlocks%s/E%d_%d.db', ...
                      nameDb, nBlocks, nameEx, b, rExt);
    util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);

    %-------------------------------------------------------------------------%
    %------- Training (Db^t), validation (Db^v), & performance (Db^f) --------%
    %-- Folds definition for learning data
    [addrFold] = util_defineFolds(1, nRep, nCls, sizeLn);

    for r = 1:nRep
        %-- When not to train (folds dedicated to validation and performance)
        foldVal = zeros(1,nRep);
        foldPfm = zeros(1,nRep);
        foldTrn = zeros(1,nRep);

         foldVal(r) = 1;                     %- Validation fold
%        foldPfm( r )               = 1;       %- Performance evaluation folds 
        foldPfm( mod(r,  nRep)+1 ) = 1;       %- Performance evaluation folds 
%         foldPfm( mod(r+1,nRep)+1 ) = 1;         %-- Performance evaluation folds 
        foldTrn = not(or(foldVal, foldPfm));  %- Training folds 

        %-------------------------- Training (Db^t) ---------------------------%
        dataOut = []; tagsOut = [];
        for f = 1:nRep
            if foldTrn(f)
                for c = 1:nCls
                    dataOut = [dataOut ; dataLn(addrFold(c,f,1) : ...
                                                addrFold(c,f+1,1)-1, :, c)];
                    tagsOut = [tagsOut ; tagsLn(addrFold(c,f,1) : ...
                                                addrFold(c,f+1,1)-1, c)];
                end
            end
        end
        
        %-- Fichier de sortie - Entra�nement
        tagsOut = tagsOut-1;
        nPatterns = size(dataOut, 1);
        nameFile  = sprintf('../database/%s/upd%dBlocks%s/%s%d^t%d.db',...
                          nameDb, nBlocks, nameEx, nameLn, b, (rExt-1)*nRep + r);
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
        
        %------------------------- Validation (Db^v) --------------------------%
        dataOut = []; tagsOut = [];
%         sizeEv = zeros(1,nCls);    %-- External db
        sizeEv = floor(sizeEt/2);    %-- External db

        %-- LTM
        for c = 1:nCls
            dataOut = [dataOut ; dataEt(1:sizeEv(c),:,c)];
            tagsOut = [tagsOut ; tagsEt(1:sizeEv(c),c)  ];
        end
        
        %-- D_t
        for f =1:nRep, if foldVal(f)
            for c = 1:nCls
                dataOut = [dataOut ; ...
                     dataLn(addrFold(c,f,1) : addrFold(c,f+1,1)-1, :, c)];
                tagsOut = [tagsOut ; ...
                     tagsLn(addrFold(c,f,1) : addrFold(c,f+1,1)-1, c)];
            end
        end, end

        %-- Fichier de sortie - Validation
        tagsOut = tagsOut-1;
        nPatterns = size(dataOut, 1);
        nameFile  = sprintf('../database/%s/upd%dBlocks%s/%s%d^v%d.db',...
                          nameDb, nBlocks, nameEx, nameLn, b, (rExt-1)*nRep + r);
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);

        %--------------------- Fitness evaluation (Db^f) ----------------------%
        dataOut = []; tagsOut = [];
  
        %-- LTM
        for c = 1:nCls
            dataOut = [dataOut ; dataEt(sizeEv(c)+1:sizeEt(c),:,c)];
            tagsOut = [tagsOut ; tagsEt(sizeEv(c)+1:sizeEt(c),c)  ];
        end
        
        %-- D_t
        for f =1:nRep, if foldPfm(f)
            for c = 1:nCls
                dataOut = [dataOut ; ...
                     dataLn(addrFold(c,f,1) : addrFold(c,f+1,1)-1, :, c)];
                tagsOut = [tagsOut ; ...
                     tagsLn(addrFold(c,f,1) : addrFold(c,f+1,1)-1, c)];
            end
        end, end

        %-- Fichier de sortie - Fitness
        tagsOut = tagsOut-1;
        nPatterns = size(dataOut, 1);
        nameFile= sprintf('../database/%s/upd%dBlocks%s/%s%d^f%d.db',...
                        nameDb, nBlocks, nameEx, nameLn, b, (rExt-1)*nRep + r);
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
    end
end
fprintf('\n');
end
fprintf('/*-------------- end - Learn ---------------*/\n');
