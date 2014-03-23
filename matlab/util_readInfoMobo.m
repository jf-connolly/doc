function [mobo] = util_readInfoMobo()

nameFile = sprintf('../database/mobof32/infoMobo.txt');


f = fopen(nameFile,'r');

    tmp = fgetl(f);  tmp = fgetl(f);  tmp = fgetl(f);
    tmp = fgetl(f);  tmp = fgetl(f);

    tmp = fread(f,49);   nClasses = fscanf(f, '%d', 1);
    tmp = fread(f,49);   nWalks   = fscanf(f, '%d', 1);
    tmp = fread(f,49);   nViews   = fscanf(f, '%d', 1);

    mobo = struct( ...
                'nViews',   nViews, ...
                'nWalks',   nWalks, ...
                'nClasses', nClasses, ...
                'nSequences', 0, ...
                'sizeCls',  zeros(1, nClasses), ...
                'info',     zeros(nViews*nWalks*nClasses,1), ...
                'addr',     zeros(nViews*nWalks*nClasses + 1,1) ...
             );
    
    tmp = fgetl(f);
    tmp = fgetl(f);
    tmp = fgetl(f);
    tmp = fgetl(f);
    
    infom = zeros(nViews,nWalks,nClasses);
    counter = 0;
    for j = 1:nWalks
        for i = 1:nViews
            for k = 1:nClasses
                counter     = counter + 1;
                infom(i,j,k) = fscanf(f, '%d', 1);
            end
        end
        tmp = fgetl(f);
    end
    mobo.nSequences = counter;

    %-- Adjust for the TEST sequence
    mobo.info = reshape(floor(infom/2),nWalks*nViews*nClasses,1);
    
    %-- Addresses of each sequences of the TEST data base!!!!!!!!!!
    mobo.addr = [0 ; cumsum(mobo.info)];
fclose(f);
