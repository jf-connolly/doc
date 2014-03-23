%------------------------------------------------------------------------------%
%-- P2 synthetique problem with hold out validation (balanced update scenario)
%------------------------------------------------------------------------------%
%-- DEFINE
HOV = 1;    TEST = 2;
A   = 1;    B    = 2;
%-- END DEFINE
%---------------------------------- Settings ----------------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
typeData      = HOV;   %- dataset type: {HOV, TEST}
sizeCls       = 1000;  %- Nb. of patterns per class
nBlocks       = 3;     %- Number of blocks
nReplications = 10;
nameDb        = '3N';  %- Database name
nameSc        = 'upd';
overlap       = '01';%        Prob d'erreur : {01, 13, 25} --%
geometry      = B;%       {A (trianbgle), B (sandwich)} --%

%------------------------------------------------------------------------------%
fprintf('\n/*----- Generating data for learning (3N) ----*/\n');

%-- Verification if the folder exists
nameFile = sprintf('upd%dBlocks%s', nBlocks, '');
util_CheckAndCreateFolder(nameDb, nameFile);

%-- Addresses to access the right data at the right time
nClasses = 3;
addrCls  = round(0 : sizeCls*nClasses/nBlocks : sizeCls*nClasses);

%------------------------------------------------------------------------------%
%---------------------------------- Equations ---------------------------------%
rho = 0.0;%                                  Corr�lation --%

switch overlap%                  D�finition de sigma^2 --%
    case '01'
        sigmaVar = 1;
    case '13'
        sigmaVar = 1.734^2;
    case '25'
        sigmaVar = 2.737^2;
    otherwise
        sigmaVar = 1;
        fprintf('/*-- Erreur : prob\n');
end

sigma_x = sigmaVar;%                         Variance x --%
sigma_y = sigmaVar;%                         Variance y --%
rho_xy  = rho * sqrt(sigma_x) * sqrt(sigma_y);
sigma   = [ sigma_x , rho_xy ;...
            rho_xy , sigma_y ];

%------------ Definition des centres et du nom ------------%
if(geometry == A)
    mu1        = [0, 0];
    mu2        = [5.116, 0];
    mu3        = [5.116/2, 5.116*sqrt(3)/2];
    sigmaVar3  = [ 1 , rho_xy ; rho_xy , 1 ];
elseif(geometry == B)
    mu1        = [0, 0];
    mu_temp    = 6.886;
    mu2        = [mu_temp/sqrt(2), mu_temp/sqrt(2)]*2;
    mu3        = [mu_temp/sqrt(2), mu_temp/sqrt(2)];
    sigma_temp = 1.885;
    sigmaVar3  = [ sigma_temp^2 , rho_xy ; rho_xy , sigma_temp^2 ];
else
    fprintf('/*-- Erreur : Variation\n');
end

%------------------------------------------------------------------------------%
%---------------------------- For ten replications ----------------------------%
for rp=1:nReplications
    %-- Initialization
    nPatterns = max(sizeCls)*3;

    class1 = mvnrnd(mu1, sigma,     nPatterns/3);
    class2 = mvnrnd(mu2, sigma,     nPatterns/3);
    class3 = mvnrnd(mu3, sigmaVar3, nPatterns/3);

    %-- Filling the data base, tags, and other informations.
    tags = [ zeros(size(class1,1),1)    ; ...
             zeros(size(class2,1),1)+1  ; ...
             zeros(size(class3,1),1)+2 ];
    data = [class1; class2; class3];


    %-- Normalisation
    data(:,1) = ( data(:,1)      - min(data(:,1)) ) / ...
                ( max(data(:,1)) - min(data(:,1)) );
    data(:,2) = ( data(:,2)      - min(data(:,2)) ) / ...
                ( max(data(:,2)) - min(data(:,2)) );

    %--------------------------------------------------------------------------%
    %----------------------------- Test data bases ----------------------------%
    if typeData == TEST
        for t = 1:nBlocks
            nameFile = sprintf('../database/%s/upd%dBlocks/T%d^%d.db',...
                              nameDb, nBlocks, t, rp);
           
            [nPatterns, nFeatures] = size(data);
            util_writeDb(nPatterns, nFeatures, nClasses, data, tags, nameFile);
        end
            
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
            tagsOut = [ tags(1         : addrVal            )  ; ...
                        tags(1+sizeCls : addrVal +sizeCls   ) ];

            %-- Train - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nameFile = sprintf('../database/%s/upd%dBlocks/B%d^t%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Validation - Data
            dataOut = [ data(addrVal           + 1 : addrFit,           :)  ;...
                        data(addrVal + addrEnd + 1 : addrFit + addrEnd, :) ];
            tagsOut = [ tags(addrVal           + 1 : addrFit             )  ;...
                        tags(addrVal + addrEnd + 1 : addrFit + addrEnd   ) ];

            %-- Validation - File creation
            [nPatterns,nFeatures] = size(dataOut);

            nameFile = sprintf('../database/%s/upd%dBlocks/B%d^v%d.db', ...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Fitness - Data
            dataOut = [ data(addrFit           + 1 : addrEnd,           :)  ;...
                        data(addrFit + addrEnd + 1 : addrEnd + addrEnd, :) ];
            tagsOut = [ tags(addrFit           + 1 : addrEnd             )  ;...
                        tags(addrFit + addrEnd + 1 : addrEnd + addrEnd   ) ];

            %-- Fitness - File creation
            [nPatterns,nFeatures] = size(dataOut);

            nameFile = sprintf('../database/%s/upd%dBlocks/B%d^f%d.db',...
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
            tagsOut = [ tags(addrStt +1           : addrVal            )  ; ...
                        tags(addrStt +1 +sizeCls  : addrVal +sizeCls   ) ];

            %-- Train - File creation
            [nPatterns,nFeatures] = size(dataOut);
            nameFile = sprintf('../database/%s/upd%dBlocks/D%d^t%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Validation - Data
            dataOut = [ data(addrVal           + 1 : addrFit,           :)  ;...
                        data(addrVal + addrEnd + 1 : addrFit + addrEnd, :) ];
            tagsOut = [ tags(addrVal           + 1 : addrFit             )  ;...
                        tags(addrVal + addrEnd + 1 : addrFit + addrEnd   ) ];

            %-- Validation - File creation
            [nPatterns,nFeatures] = size(dataOut);

            nameFile = sprintf('../database/%s/upd%dBlocks/D%d^v%d.db', ...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);

            %-- Fitness - Data
            dataOut = [ data(addrFit           + 1 : addrEnd,           :)  ;...
                        data(addrFit + addrEnd + 1 : addrEnd + addrEnd, :) ];
            tagsOut = [ tags(addrFit           + 1 : addrEnd          )  ;...
                        tags(addrFit + addrEnd + 1 : addrEnd + addrEnd) ];

            %-- Fitness - File creation
            [nPatterns,nFeatures] = size(dataOut);

            nameFile = sprintf('../database/%s/upd%dBlocks/D%d^f%d.db',...
                               nameDb, nBlocks, t, rp);
            util_writeDb(nPatterns, nFeatures, nClasses, dataOut, tagsOut, ...
                          nameFile);
        end
    else
        fprintf('/*-- Error : data type\n')
    end
end %-- for : nReplications --%

fprintf('/*----------------- End - 3N -----------------*/\n');
