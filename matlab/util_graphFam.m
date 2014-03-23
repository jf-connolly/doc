function [mapping grid] = util_mappingFam(fam, resolution)

%-- Initialization
grid    = 0:resolution:1;
T       = zeros(1,fam.sizeNn);
mapping = zeros(size(grid,2), size(grid,2));

%----------------- Mapping of the predictions in the 2D space -----------------%
for posY = 1 +1 : size(y,2) -1
    for posX = 1 +1 : size(x,2) -1
        
        %-- Data line vector (one point)
        A = [ grid(posX) grid(posY) 1-grid(posX) 1-grid(posY) ];
        
        %-- Tj = |min(A,Wij)| / (alpha+|Wij|)
        for j = 1:sizeNn
            normAW =  sum( min( A, fam.W(j,:) ) );    %-- |min(A,Wij)|
            T(j) = normAW / fam.normWalpha(j);        %-- T(j)
        end

        %-- Winning node 
        [tempMax J] = max(T);
        
        %-- Chossen class
        [tempMax ck] = max( fam.Wab(J,:) );
        
        %-- Assign it to the class matrix
        mapping(posX, posY) = ck;
    end %** for : posX **%
end %** for : posY **%