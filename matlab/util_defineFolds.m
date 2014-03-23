function [addrFold] = util_defineFolds(nBlocs, nRep, nClasses, sizeCls)

%------------------------------- Initialisation -------------------------------%
%-- Cr�� un matrice 3D du nombre de donn�es selon : 
%-- la classe, le bloc et le folds.
%-- i^e ligne     : i^e classe              Folds
%-- j^e colonne   : j^e fold                __________    .
%-- k^e dimension : k^e bloc       Classes | Bloc 2   |  .
%--                                     ___|______    | .
%--                                    | Bloc 1   |___|
%--                                    |          |
%--                                    |__________|

%------------------ Nb de donn�es par bloc pour chaque classe -----------------%
nData_Bloc       = transpose(sizeCls)/nBlocs;
nData_Bloc_Modif = zeros(nClasses,nBlocs);
for(c=1:nClasses)
    for(b=1:nBlocs)
        nData_Bloc_Modif(c,b) = round(nData_Bloc(c) + 0.51 - b/nBlocs);
    end
end

%------------- Nb de donn�es par fold pour chaque bloc et classe --------------%
nData_Temp = nData_Bloc_Modif/nRep;
% dataFold   = zeros(nBlocs, nClasses, nRep);
for(k=1:nBlocs)
    for(i=1:nClasses)
        for(j=1:nRep)
            dataFold(i,j,k) = round(nData_Temp(i,k) + 0.51 - j/nRep);
        end
    end
end

%-------------------------------- Les adresses --------------------------------%
addrFoldTmp = cumsum(dataFold, 2);

% addrFold    = [ones(nClasses,1), addrFoldTmp(:,:,b)];
% mask = and(addrFold,addrFold);
% addrFold = [];

%-- Adresses definition
% addrFold = zeros(nBlocs, nClasses, nRep);
for b = 1:nBlocs
    addrFold(:,:,b) = [ones(nClasses,1), addrFoldTmp(:,:,b)+1];
%     mask = and(addrFold,addrFold);
%     addrFold(:,:,b)
    
    if b ~= 1
        for k = 1:nRep+1
            addrFold(:,k,b) = addrFold(:,k,b) + addrFold(:,nRep +1,b-1)-1;
        end
    end
end

%------------------------------------ Test ------------------------------------%
if (addrFold(:, nRep +1, nBlocs)-1 == transpose(sizeCls))
%     fprintf('/*-- Definition corrrecte!!! \n');
%     fprintf('/*-- (A prononcer comme dans les Pierrafeu) \n');
else
    fprintf('/*-- C`est nul... \n\n');
end
