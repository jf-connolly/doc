%------------------------------------------------------------------------------%
%-- P2 synthetique problem with hold out validation (balanced update scenario)
%------------------------------------------------------------------------------%
%-- DEFINE
HOV = 1;    TEST = 2;
%-- END DEFINE
%---------------------------------- Settings ----------------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
typeData      = HOV;    %- dataset type: {HOV, TEST}
sizeCls       = 300;    %- Nb. of patterns per class
nBlocks       = 1;      %- Number of blocks
nReplications = 10;
nameDb        = 'cis';  %- Database name
nameSc        = 'upd';

%------------------------------------------------------------------------------%
fprintf('\n/*------- Generating data for learning (CIS) -------*/\n');

%-- Verification if the folder exists
nameFile = sprintf('upd%dBlocks%s', nBlocks, 'Ext');
util_CheckAndCreateFolder(nameDb, nameFile);

%-- Addresses to access the right data at the right time
nClasses = 2;
addrCls  = round(0 : max(sizeCls)/(nBlocks) : max(sizeCls));

%------------------------------------------------------------------------------%
%---------------------------------- Equations ---------------------------------%
%-- Circle with an area of 0.5
%     (x-0.5)^2 + (y-0.5)^2 = r^2, where
r = 0.398942;

%------------------------------------------------------------------------------%
%---------------------------- For ten replications ----------------------------%
for rp=1:nReplications
    %-- Initialization
    nameFile   = [];
    nPatterns  = max(sizeCls)*6;
    
    class1     = zeros(sizeCls,2);
    class2     = zeros(sizeCls,2);
    sizeClass1 = 0;
    sizeClass2 = 0;

    data = rand(nPatterns, 2);
   
    %--- Cr�ation des classes et du vecteur d'�tiquettes --%
    for n = 1:nPatterns
        %-- If we have all our data, go to the data sets creation
        if sizeClass1 + sizeClass2 >= sizeCls *2
            break;
        end
                    
        %-- Else, we fill the data base
        rTest = (data(n,1)-0.5)^2 + (data(n,2)-0.5)^2;
        
        %-- Class 1            
        if rTest < r^2
            %-- Add pattern
            if sizeClass1 < max(sizeCls)
                sizeClass1           = sizeClass1 + 1;
                class1(sizeClass1,:) = data(n,:);
            end
        %-- Class 2
        else
            %-- Add pattern
            if sizeClass2 < max(sizeCls)
                sizeClass2           = sizeClass2 + 1;
                class2(sizeClass2,:) = data(n,:);
            end
        end
    end

    %-- Filling the data base, tags, and other informations.
    tags = [ zeros(size(class1,1),1) ; ones(size(class2,1),1) ];
    data = [class1; class2];
    
    %--------------------------------------------------------------------------%
    %----------------------------- Test data bases ----------------------------%
    if typeData == TEST
        nameFile = sprintf('../database/%s/test.db', nameDb);

        [nPatterns, nFeatures] = size(data);
        util_writeDb(nPatterns, nFeatures, nClasses, data, tags, nameFile);
            
    %--------------------------------------------------------------------------%
    %---------------------------- Learn data bases ----------------------------%
    elseif typeData == HOV
        for(t = 1:nBlocks)

            %-------------------- Batch learning data sets --------------------%
            addrVal = round(addrCls(t+1)*2/3);
            addrFit = round(addrCls(t+1)*5/6);
            addrEnd = sizeCls;
            
            %-- Train - Data
            dataOut = [ data(1         : addrVal         , :)  ; ...
                        data(1+sizeCls : addrVal +sizeCls, :) ];
            tagsOut = [ zeros(addrVal, 1) ; ones(addrVal, 1) ];

            %-- Train - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nameFile = sprintf('../database/%s/upd%dBlocksExt/B%d^t%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Validation - Data
            dataOut = [ data(addrVal           + 1 : addrFit,           :)  ;...
                        data(addrVal + addrEnd + 1 : addrFit + addrEnd, :) ];
            tagsOut = [ zeros(addrFit-addrVal, 1) ; ones(addrFit-addrVal, 1) ];

            %-- Validation - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nClasses = max(tagsOut)+1;

            nameFile = sprintf('../database/%s/upd%dBlocksExt/B%d^v%d.db', ...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Fitness - Data
            dataOut = [ data(addrFit           + 1 : addrEnd,           :)  ;...
                        data(addrFit + addrEnd + 1 : addrEnd + addrEnd, :) ];
            tagsOut = [ zeros(addrEnd-addrFit, 1) ; ones(addrEnd-addrFit, 1) ];

            %-- Fitness - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nClasses = max(tagsOut)+1;

            nameFile = sprintf('../database/%s/upd%dBlocksExt/B%d^f%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %----------------- Incremental learning data sets -----------------%
            addrStt = round(addrCls(t));
            addrVal = addrStt + round( (addrCls(t+1) - addrStt) *2/3 );
            addrFit = addrStt + round( (addrCls(t+1) - addrStt) *5/6 );
            addrEnd = addrCls(t+1);

            %-- Train - Data
            dataOut = [ data(addrStt +1           : addrVal         , :)  ; ...
                        data(addrStt +1 +sizeCls  : addrVal +sizeCls, :) ];
            tagsOut = [ zeros(addrVal-addrStt, 1) ; ones(addrVal-addrStt, 1) ];

            %-- Train - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nameFile = sprintf('../database/%s/upd%dBlocksExt/D%d^t%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Validation - Data
            dataOut = [ data(addrVal           + 1 : addrFit,           :)  ;...
                        data(addrVal + addrEnd + 1 : addrFit + addrEnd, :) ];
            tagsOut = [ zeros(addrFit-addrVal, 1) ; ones(addrFit-addrVal, 1) ];

            %-- Validation - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nClasses = max(tagsOut)+1;

            nameFile = sprintf('../database/%s/upd%dBlocksExt/D%d^v%d.db', ...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Fitness - Data
            dataOut = [ data(addrFit           + 1 : addrEnd,           :)  ;...
                        data(addrFit + addrEnd + 1 : addrEnd + addrEnd, :) ];
            tagsOut = [ zeros(addrEnd-addrFit, 1) ; ones(addrEnd-addrFit, 1) ];

            %-- Fitness - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nClasses = max(tagsOut)+1;

            nameFile = sprintf('../database/%s/upd%dBlocksExt/D%d^f%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);
        end
    else
        fprintf('/*-- Error : data type\n')
    end
end %-- for : nReplications --%

fprintf('/*-------------------- End - CIS -------------------*/\n');
