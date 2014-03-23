%------------------------------------------------------------------------------%
%-- Scenario : Class enrollment with LTM
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
nBlocks = 24;
nameDb  = 'mobof32';
nRep    = 10;
nRepExt = 5;

NEW   = 0;        %-- New pattern presentation order?
BATCH = 1;        %-- Batch learning data sets?

%-- External db parameters
ratioEt = 1/5;      %-- Fraction of Dt going to E_t
normEt  = 30;       %-- Maximum number of patterns  / class

%-------------------------- Initialization sequence ---------------------------%
%-- Names and ratio (for batch learning)
if BATCH,    nameLn = 'B';    ratioEt = 0;  else  nameLn = 'D';  end
if ratioEt,  nameEx = 'Ext';                else  nameEx = '';   end

%-- Number of classes per block
nameFile = sprintf('../database/%s/learn.db', nameDb);
dbase = util_readHeader(nameFile);

%-- Test for the number of blocks
if nBlocks > dbase.nClasses-1, fprintf('/*-- To much blocks...\n'); return; end

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

%-- Reading the db
fprintf('/*-- Class addition - learning data sets ---*/\n');
fprintf('/*-- Data reading\n');
nameFile = sprintf('../database/%s/learn.db', nameDb);

dbase = util_readDb(nameFile);
nPatterns = size(dbase.data,1);
nFeatures = size(dbase.data,2);
nCls      = max(dbase.tags) + 1;

%-- Nb of patterns per classes
sizeClsLoad = zeros(1, nCls);
for i = 1:nPatterns
    sizeClsLoad(dbase.tags(i)+1) = sizeClsLoad(dbase.tags(i)+1) + 1;
end
addrCls = [1, (cumsum(sizeClsLoad)+1)];

%----------------------------- Classes separation -----------------------------%
data    = zeros(max(sizeClsLoad), nFeatures, nCls);
tags    = zeros(max(sizeClsLoad), nCls);
sizeCls = zeros(1, nCls);

%-- Separation (tags = tags+1)
for k = 1:nCls
    data(1:sizeClsLoad(k), :, k) = dbase.data(addrCls(k):addrCls(k+1)-1, :);
    tags(1:sizeClsLoad(k), k)    = ones(sizeClsLoad(k),1) * (k);
    sizeCls(1,k) = sizeClsLoad(k);
end

%-- Verification if the folder exists
nameFile = sprintf('add%dBlocks%s', nBlocks, nameEx);
util_CheckAndCreateFolder(nameDb, nameFile);

%------------------------------------------------------------------------------%
%-------------------------- Learning blocs creation ---------------------------%
for rExt = 1:nRepExt
fprintf('/*-- External replication %d/%d - blocks: ', rExt, nRepExt);

%-- E initialization
dataNw = zeros(1, nFeatures, nCls);         tagsNw = zeros(1, nCls);
dataEt = zeros(normEt, nFeatures, nCls);    tagsEt = zeros(normEt, nCls);

for t = 1:nBlocks
    fprintf('%d ', t);

    %------------------------- Learning block t (D_t) -------------------------%
    if t == 1 || ~BATCH
        sizeDt = zeros(1, nCls);
        dataDt = zeros(1, nFeatures, nCls);   tagsDt = zeros(1, nCls);
    end

    for cls = 1:nClsPblock(t)
        cCls = order( rExt, nClsPblockSum(t) + cls);

        sizeDt(cCls)                    = sizeCls(cCls);
        dataDt(1:sizeDt(cCls), :, cCls) = data(1:sizeDt(cCls), :, cCls);
        tagsDt(1:sizeDt(cCls),    cCls) = tags(1:sizeDt(cCls),    cCls);
    end
    
    if ~BATCH
        %----------------- Spliting data for E_t and learning -----------------%
        %-- Number of pattern for E_t and for learning
        sizeNw = round(sizeDt*ratioEt);
        sizeLn = sizeDt - sizeNw;

        %-- Data for E and for learning
        dataNw = zeros(1, nFeatures, nCls); tagsNw = zeros(1, nCls);
        dataLn = zeros(1, nFeatures, nCls); tagsLn = zeros(1, nCls);
        for k = 1:nCls
            dataNw(1:sizeNw(k), :, k) = dataDt(1:sizeNw(k), :, k);
            tagsNw(1:sizeNw(k), k)    = tagsDt(1:sizeNw(k), k);

            dataLn(1:sizeLn(k),:,k) = dataDt(sizeNw(k)+1:sizeDt(k),:,k);
            tagsLn(1:sizeLn(k),k)   = tagsDt(sizeNw(k)+1:sizeDt(k),k);
        end

        %----------------- E_b creation & writing it in a file ----------------%
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
        nameFile  = sprintf('../database/%s/add%dBlocks%s/E%d_%d.db', ...
                          nameDb, nBlocks, nameEx, t, rExt);
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
    else
        %------------------------- For batch learning -------------------------%
        dataLn = dataDt;
        tagsLn = tagsDt;
        sizeLn = sizeDt;
		sizeEt = zeros(1,nCls);
    end
    

    %--------------------------------------------------------------------------%
    %------- Training (Dt^t), validation (Dt^v), & performance (Dt^p) ---------%
    %-- Folds definition for learning data
    [addrFold] = util_defineFolds(1, nRep, nCls, sizeLn);

    for r = 1:nRep
        %-- When not to train (folds dedicated to validation and performance)
        foldVal = zeros(1,nRep);
        foldPfm = zeros(1,nRep);
        foldTrn = zeros(1,nRep);

        foldVal(r) = 1;                       %-- Validation fold
        foldPfm( mod(r,  nRep)+1 ) = 1;         %-- Performance evaluation folds 
        foldPfm( mod(r+1,nRep)+1 ) = 1;         %-- Performance evaluation folds 
        foldTrn = not(or(foldVal, foldPfm));  %-- Training folds 

        %-------------------------- Training (Dt^t) ---------------------------%
        dataOut = [];   tagsOut = [];
        for kTr = 1:nRep
            if foldTrn(kTr)
                for k = 1:nCls
                    dataOut = [dataOut ; dataLn(addrFold(k,kTr) : ...
                                                addrFold(k,kTr+1)-1, :, k)];
                    tagsOut = [tagsOut ; tagsLn(addrFold(k,kTr) : ...
                                                addrFold(k,kTr+1)-1, k)];
                end
            end
        end
        sizeOut = util_sizeofDb(tagsOut);
        tagsOut = tagsOut-1;
        
        %-- Fichier de sortie - Entra�nement
        nPatterns = size(dataOut, 1);
        nameFile  = sprintf('../database/%s/add%dBlocks%s/%s%d^t%d.db',...
                         nameDb, nBlocks, nameEx, nameLn, t, (rExt-1)*nRep + r);
		if ~BATCH || t == nBlocks
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
        end
        %------------------------- Validation (Dt^v) --------------------------%
        dataOut = [];   tagsOut = [];
        [x addr] = max( foldVal );
        sizeEv   = floor( sizeEt/2 );    %-- External db

        for k = 1:nCls
            dataOut = [dataOut ; dataEt(1:sizeEv(k),:,k) ; ...
                         dataLn(addrFold(k,addr) : addrFold(k,addr+1)-1, :, k)];
            tagsOut = [tagsOut ; tagsEt(1:sizeEv(k),k) ; ...
                         tagsLn(addrFold(k,addr) : addrFold(k,addr+1)-1, k)];
        end
        tagsOut = tagsOut-1;

        %-- Fichier de sortie - Validation
        nPatterns = size(dataOut, 1);
        nameFile  = sprintf('../database/%s/add%dBlocks%s/%s%d^v%d.db',...
                         nameDb, nBlocks, nameEx, nameLn, t, (rExt-1)*nRep + r);
		if ~BATCH || t == nBlocks
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
		end
        %------------------- Performance evaluation (Dt^p) --------------------%
        dataOut = [];   tagsOut = [];
        [x addr] = max(foldPfm);
        sizeEp = sizeEt-sizeEv;
  
        for k = 1:nCls
            dataOut = [dataOut ; dataEt(sizeEv(k)+1:sizeEt(k),:,k) ; ...
                         dataLn(addrFold(k,addr) : addrFold(k,addr+1)-1, :, k)];
            tagsOut = [tagsOut ; tagsEt(sizeEv(k)+1:sizeEt(k),k) ; ...
                         tagsLn(addrFold(k,addr) : addrFold(k,addr+1)-1, k)];
        end
        tagsOut = tagsOut-1;

        %-- Fichier de sortie - Performance evaluation
        nPatterns = size(dataOut, 1);
        nameFile= sprintf('../database/%s/add%dBlocks%s/%s%d^f%d.db',...
                        nameDb, nBlocks, nameEx, nameLn, t, (rExt-1)*nRep + r);
		if ~BATCH || t == nBlocks
        util_writeDb(nPatterns, nFeatures, nCls, dataOut, tagsOut, nameFile);
		end
    end
end
fprintf('\n');
end
fprintf('/*-------------- End - Learn ---------------*/\n');
