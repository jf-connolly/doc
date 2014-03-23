function [] = util_CheckAndCreateFolder(nameDb, nameTest)

%-- File opening trial
nameFile = sprintf('../database/%s/%s/D1^t1.db', nameDb, nameTest);
fInc = fopen(nameFile,'r');

nameFile = sprintf('../database/%s/%s/B1^t1.db', nameDb, nameTest);
fBth = fopen(nameFile,'r');

if fInc == -1 && fBth == -1
    %-- If it does not exist, it is created
    nameFile = sprintf('../database/%s', nameDb);
    mkdir( nameFile, nameTest);
end

%-- Else we close the opened file(s)
if fInc ~= -1,  fclose(fInc);  end
if fBth ~= -1,  fclose(fBth);  end



