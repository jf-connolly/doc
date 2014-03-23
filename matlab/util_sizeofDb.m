function sizeDb = util_sizeofDb(data)

nClasses = size(data,3);
sizeDb   = zeros(1,nClasses);

for k = 1:nClasses
    for i = 1:size(data,1)
        if data(i,1,k)
            sizeDb(k) = sizeDb(k) + 1;
        end
    end
end