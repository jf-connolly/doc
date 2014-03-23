function [dbase] = util_readDb(nameFile)

dbase = struct('nPatterns',[],'nFeatures',[],'nClasses',[],'data',[],'tags',[]);

f = fopen(nameFile,'r');

    dbase.nPatterns = fread(f, 1, 'int32');%   nombre d'exemples    --%
    dbase.nFeatures = fread(f, 1, 'int32');%   nombre de dimensions --%
    dbase.nClasses  = fread(f, 1, 'int32');%   nombre de classes    --%

    dbase.data = zeros(dbase.nPatterns, dbase.nFeatures); 
    dbase.tags = zeros(dbase.nPatterns, 1);

    for(i=1:dbase.nPatterns)
        dbase.data(i,:) = fread(f, dbase.nFeatures, 'single');%les patrons --%
    end

    dbase.tags = fread(f, dbase.nPatterns, 'int32');% �tiquettes
    
fclose(f);
