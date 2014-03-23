%------------------------------------------------------------------------------%
%-- Shows info & chosen faces of real face data bases
%------------------------------------------------------------------------------%

%---------------------------------- Settings ----------------------------------%
%     __
%    |  |
%    |  |
%   _|  |_
%   \    /
%    \  /
%     \/
%
nameDb  = 'mobof64';  %-- {'cnrc', 'mobo'}
nameSc  = 'upd';
nBlocks =  26;
ext     = '';
nameLn  = 'B';
nameTp  = 't';
t       =  26;
r       =  1;


graph   = 0;       %-- {0,1}
pattern = 3;

fprintf('\n/*---------- Start - syn db graph -----------*/\n');

%------------------------------ Reading the data ------------------------------%
%-- Reading

if     nameLn == 'T'
    nameFile = sprintf('../database/%s/test.db', nameDb);
elseif nameLn == 'L'
    nameFile = sprintf('../database/%s/learn.db', nameDb);
else
    nameFile = sprintf('../database/%s/%s%dBlocks%s/%s%d^%s%d.db',...
                        nameDb, nameSc, nBlocks, ext, nameLn, t, nameTp, r);
end

dbase = util_readDb(nameFile);

sizeCls = zeros(dbase.nClasses, 1);
for (i=1:dbase.nPatterns)
    sizeCls(dbase.tags(i)+1,1) = sizeCls(dbase.tags(i)+1,1) + 1;
end
fprintf('\n/*-- Number of patterns - %d\n\n', dbase.nPatterns);

fprintf('/*-- Number of patterns per class - \n   ');
for i = 1:dbase.nClasses
    fprintf('%d ', sizeCls(i));
end

fprintf('\n\n');

%---------------------------- Showing a given face ----------------------------%
if graph
  
    p = pattern    
%     for p = 1 : dbase.nPatterns
    close all;
    
    fprintf('/*-- Classe %d\n', dbase.tags(pattern)+1);
%     fprintf('/*-- Classe %d\n', dbase.tags(p)+1);

    im=[];
    for(i=1:sqrt(dbase.nFeatures))
        for(j=1:sqrt(dbase.nFeatures))
            im(i,j) = dbase.data(pattern, sqrt(dbase.nFeatures)*(i-1)+j);
            im(i,j) = dbase.data(p, sqrt(dbase.nFeatures)*(i-1)+j);
        end
    end
    
    im(:,:,2) = im;
    im(:,:,3) = im(:,:,1);
    fig_1 = figure(1);
    posScreen = [-1400 600 200 200];
    set(fig_1, 'Position', posScreen, 'PaperPosition', [0 0 3 3]);
    axes('Position', [0,0,1,1])
    g = image(im);
    
    set(gca, 'XTick', [], 'YTick', []);
    
%     nameFile = sprintf('../figures/databases/%s/test/%d(class%d).png', ...
%                        nameDb, p, dbase.tags(p)+1);
%     saveas(fig_1, nameFile);
    bp=1;
%     end
end
fprintf('/*----------- End - real db graph -----------*/\n\n');
