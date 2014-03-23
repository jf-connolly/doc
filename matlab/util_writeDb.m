function []=util_writeDb(nPatterns, nFeatures, nClasses, data, tags, nameFile)
                    
f = fopen(nameFile,'w');

    fwrite(f, nPatterns, 'int32');
    fwrite(f, nFeatures, 'int32');
    fwrite(f, nClasses,  'int32');

    for(i=1:nPatterns)
        fwrite(f, data(i,:), 'single');
    end

    fwrite(f, tags,'int32');

    %-- For test purposes
    % Size_of_File = nPatterns *( 8*nFeatures +4 ) +12;
fclose(f);
