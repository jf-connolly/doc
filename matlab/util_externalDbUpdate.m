function [dataRd, tagsRd, dataEt, tagsEt] = util_externalDbUpdate...
                                                (dataEt, tagsEt, dataNw, tagsNw)

maxSizeEb = size(dataEt,1);
nFeatures = size(dataEt,2);
nCls      = size(dataEt,3);

%-- Data set sizes
sizeNw = util_sizeofDb(dataNw);
sizeEt = util_sizeofDb(dataEt);

%------------------------------- External Db (E) ------------------------------%
nTrans = zeros(nCls,1);  nUpdate = zeros(nCls,1);  nReturned = zeros(nCls,1);
maxsizeEt = size(dataEt,1);

for k = 1:nCls
    %-- Number of patterns transfered, used for update, and returned
    nTrans(k)    = min([ maxsizeEt-sizeEt(k), sizeNw(k) ]); 
    nUpdate(k)   = min([ sizeNw(k)-nTrans(k), sizeEt(k) ]);
    nReturned(k) = sizeNw(k)-nTrans(k)-nUpdate(k);
end

%-- Data returned to D_t (initialization)
dataRd = zeros(max(nReturned), nFeatures, nCls);
tagsRd = zeros(max(nReturned), nCls);

for k = 1:nCls
    %-- Update
    if nUpdate(k)
        %-- 1) Push down the patterns.
        dataEt(1:sizeEt(k)-nUpdate(k),:,k) = dataEt(nUpdate(k)+1:sizeEt(k),:,k);
        tagsEt(1:sizeEt(k)-nUpdate(k),k)   = tagsEt(nUpdate(k)+1:sizeEt(k),k);
        
        %-- 2) Update the old patterns.
        dataEt(sizeEt(k)-nUpdate(k)+1:sizeEt(k),:,k) = dataNw(1:nUpdate(k),:,k);
        tagsEt(sizeEt(k)-nUpdate(k)+1:sizeEt(k),k)   = tagsNw(1:nUpdate(k),k);
    end

    %-- Transfer: New data -> E_b.
    if nTrans(k)
        dataEt(sizeEt(k)+1:sizeEt(k)+nTrans(k), :, k) = ...
                                  dataNw(nUpdate(k)+1:nUpdate(k)+nTrans(k),:,k);
        tagsEt(sizeEt(k)+1:sizeEt(k)+nTrans(k), k)    = ...
                                  tagsNw(nUpdate(k)+1:nUpdate(k)+nTrans(k),k);
        sizeEt(k) = sizeEt(k) + nTrans(k);
    end
    
    %-- Data returned to D_t
    if nReturned(k)
        dataRd(1:nReturned(k), :, k) = ...
                                 dataNw(nUpdate(k)+nTrans(k)+1:sizeNw(k), :, k);
        tagsRd(1:nReturned(k), k)   = ...
                                 tagsNw(nUpdate(k)+nTrans(k)+1:sizeNw(k), k);
    end
end

















