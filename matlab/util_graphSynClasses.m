function [] = util_graphSynClasses(nameDb)

if(nameDb(1:2) == 'p2')
    x = 10 *(0:0.001:1);
    y = [ ( 2 *sin(x) +5 )                 /10; ...
          ( (x-2).^2  +1 )                 /10; ...
          ( -0.1*(x).^2 +0.6*sin(4*x) +8)  /10; ...
          ( ((x-10).^2)/2 + 7.902 )        /10 ];
      
    x = x/10;
    p = plot(x,y(1,:),'k', x,y(2,:),'k',x,y(3,:),'k', x,y(4,:),'k');
    
% elseif strcmp(nameDb(1:2), '3_') || strcmp(nameDb(1:2),'2_')
%     x1 = 0   : 0.001 : 0.5;
%     y1 = 0.7387 - x1/sqrt(3);
%     x2 = 0.5 : 0.001 : 1;
%     y2 = 0.1613 + x2/sqrt(3);
%     y3 = 0   : 0.001 : 0.45;
%     x3 = ones(1, size(y3,2))/2;
%     p = plot(x1,y1,'k', x2,y2,'k', x3,y3,'k');
% elseif(nameDb(1:2) == '3b' | nameDb(1:2) == '2_')
%     x1 = 0     : 0.001 : 0.666;
%     x2 = 0.333 : 0.001 : 0.999;
%     y1 = -x1 + 0.666;
%     y2 = -x2 + 0.666*2;
%     p = plot(x1,y1,'k', x2,y2,'k');
elseif(nameDb(1:2) == 'ci')
    %-- Eq. du cercle : (x-0.5)^2 + (y-0.5)^2 = r^2           
    r = 0.398942;
    x  = 0.5-r : 0.001 : 0.5+r+0.001;
    y1 = sqrt( r^2 - (x-0.5).^2 ) + 0.5;
    y2 = -sqrt( r^2 - (x-0.5).^2 ) + 0.5;
    p = plot(x,y1,'k',x,y2,'k');
% elseif(nameDb(1:2) == '2N')
%     x = 0:0.001:1;
%     y1 = 1-x;
%     p = plot(x,y1,'k');
% elseif(nameDb(1:2) == 'xo')
%     if(size(nameDb,2) == 8)
%         x = 0:0.001:1;
%         y = 1-x;
%         p = plot(x,x,'k', x,y,'k');
%     else
%         x = 0:0.001:1;
%         y = ones(1,size(x,2))*0.5;
%         p = plot(x,y,'-k', y,x,'-k');
%     end
else
    fprintf('/-- Erreur : nameDb\n')
end

set(p,'LineWidth', 1, 'LineStyle', '--');

