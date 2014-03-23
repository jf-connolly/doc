function [dbase] = util_readHeader(nameFile)

dbase = struct('nPatterns',[],'nFeatures',[],'nClasses',[],'data',[],'tags',[]);

f = fopen(nameFile,'r');

    dbase.nPatterns = fread(f, 1, 'int32');
    dbase.nFeatures = fread(f, 1, 'int32');
    dbase.nClasses  = fread(f, 1, 'int32');
    
fclose(f);