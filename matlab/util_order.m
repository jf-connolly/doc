function [orderAll, success] = util_order( NEW, nRepExt, nBlocks, nameDb )

nBlocksTot = sum(nBlocks); 
nameFile   = sprintf('../database/%s/orders.txt', nameDb);

if NEW
    %-- Original vector
    nClasses   = nBlocks;
    orderInit  = 1:nClasses;
%     nClasses   = size(nBlocks,1);
%     addrCls    = [0, transpose(cumsum(nBlocks))];
%     orderInit  = zeros (1, nBlocksTot);
%     for c = 1: nClasses
%         orderInit(addrCls(c)+1:addrCls(c+1)) = ones(1,nBlocks(c))*c;
%     end

    %-- Shuffles everything
    orderAll = repmat(orderInit, nRepExt, 1);
    for k = 1:nRepExt
        i = nBlocksTot;
        while i > 0
            %-- Random integers: nRandom = {1:i}
            nRandom = ceil(rand(1)*i);
            %-- Swap
            temp                = orderAll(k,nRandom);
            orderAll(k,nRandom) = orderAll(k,i);
            orderAll(k,i)       = temp;
            %-- D'hu
            i = i-1;
        end
    end

    %-- Save and print it
    fprintf('/*-- Class presentation orders:\n');
    f = fopen(nameFile,'w');

    %-- Writes it in a file and prints on screen
    for k = 1:nRepExt
        fprintf('   External replication %d - ', k);
        for i = 1:nBlocksTot
            fprintf(f, '%d ', orderAll(k,i));   fprintf('%d ', orderAll(k,i));
        end
        
        fprintf(f, '\n');   fprintf('\n');
    end

else
    %-- Read old class presentation order
    f = fopen(nameFile,'r');
    if f == -1
        fprintf('/*-- No order defined!!!\n');
        success = 0;
    else
        %-- ... if able
        orderAll = transpose(fscanf(f, '%d', [nBlocksTot, nRepExt]));
        nData    = size(orderAll,2);
        
        %-- Prints on screen
        for k = 1:nRepExt
            fprintf('   External replication %d - ', k);
            for i = 1:nData,   fprintf('%d ', orderAll(k,i));   end
            fprintf('\n');
        end
    end
end
    
%-- Closure
fclose(f);
success = 1;
% fprintf('\n');
