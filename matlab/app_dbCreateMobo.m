%----------------- Lecture et infos de la base d'entrainement -----------------%
fprintf('/*------ Mobo train and test creation ------*/\n');
rand('twister',sum(100*clock))

nClsAfr    = 25;
sizeClsAfr = 200;
nItems     = 15;  

fprintf('/*-- Data reading - Learn \n');

nameFile = sprintf('../database/mobof32/learn.db');
dbLearn  = util_readDb(nameFile);

nameFile = sprintf('../database/mobof32/test.db');
dbTest   = util_readDb(nameFile);


nameFile = sprintf('../database/mobof64/test.db');
[nCls nWalks nViews sizeCls info addr] = util_readInfoMobo();

%-- Check
if addr(size(addr,2))-1 ~= size(dataLoad,1),  fprintf('/*-- Error\n');  end

%----------------------------- S�paration des Cls -----------------------------%
fprintf('/*-- Class separation - Learn\n');

%-- New data with 10 classes
[tmp index] = sort(sizeCls,'descend');
nPatAfr = 0;
% for c = 1:nClsAfr,    nPatAfr = nPatAfr + sizeClsLoad(index(c));    end
% 
% dataLn = zeros(ceil(nPatAfr/2), 576);    tagsLn = zeros(ceil(nPatAfr/2), 1);
% dataTt = zeros(ceil(nPatAfr/2), 576);    tagsTt = zeros(ceil(nPatAfr/2), 1);

dataTp = zeros(2*nItems, nFeatures);
tagsTp = zeros(2*nItems,1);
dataLn = zeros(nItems*nViews*nWalks*nClsAfr, nFeatures);
tagsLn = zeros(1,nItems*nViews*nWalks*nClsAfr);
dataTt = zeros(nItems*nViews*nWalks*nClsAfr, nFeatures);
tagsTt = zeros(1,nItems*nViews*nWalks*nClsAfr);
nSamples = 0;
r=[];
for k = 1:nClsAfr
    for j = 1:nViews  %-- Separate the views and keep the walks together
      
        %-- Check of every thing is ok
        flagMiss = zeros(1,nWalks);
        for i = 1:nWalks
            if info(i,j,index(k)) < 2*nItems,      flagMiss(i) = 1;   end
        end
        
        %-- Keep the walks together
        for i = 1:nWalks
            
            if ~flagMiss(i)
                a      = (index(k)-1) *nWalks *nViews + (i-1) *nViews + j;
                aStart = addr(a);
                aEnd   = addr(a) +2*nItems -1;
            else
                %-- Where to get stuff
                flagGet = zeros(1,nWalks);
                for ig = 1:nWalks
                    if info(ig,j,index(k)) >= 4*nItems && ig~=i
                        flagGet(ig) = 1;
                    end
                end
                %-- If we can't get any patterns
                if ~sum(flagGet), fprintf('/*-- Error - split\n');   end
                
                %-- Choose another source randomly
                r = ceil(nWalks.*rand(1,1));
                while ~flagGet(r),   r = ceil(nWalks.*rand(1,1));   end
                                
                a      = (index(k)-1) *nWalks *nViews + (r-1) *nViews + j;
                aStart = addr(a) +2*nItems;
                aEnd   = addr(a) +4*nItems -1;
            end

            %-- Temp
            dataTp = dataLoad(aStart:aEnd, :);
            tagsTp = zeros(2*nItems,1) -1 + k;
            
            %-- Learn -&- test data sets
            for n = 1:nItems
                nSamples = nSamples +1;

                dataLn(nSamples, :) = dataTp(2*n-1, :);
                tagsLn(1, nSamples) = tagsTp(2*n-1);
                
                dataTt(nSamples, :) = dataTp(2*n, :);
                tagsTt(1, nSamples) = tagsTp(2*n);
            end

        end
    end
end

fprintf('/*-- Db creation - Learn \n');
nPatterns = size(dataLn, 1);
nameFile  = sprintf('../database/mobo64/learn.db');
util_createDb(nPatterns, nFeatures, nClsAfr, dataLn, tagsLn, nameFile);

fprintf('/*-- Db creation - Test  \n');
nPatterns = size(dataTt, 1);
nameFile  = sprintf('../database/mobo64/test.db');
util_createDb(nPatterns, nFeatures, nClsAfr, dataTt, tagsTt, nameFile);
