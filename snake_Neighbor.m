function neighbor = snake_Neighbor(sizeX,sizeY,x,y,V)
%
% 得到轮廓上某点[x,y]周围的点的坐标（8个），[sizeX,sizeY]是图像大小
%
% -----------------------------------------------------------------------------------------------------------

% 邻域的大小，往x方向两边扩展nX个象素，往y方向两边扩展ny个象素
R = 3;
nX = R;
nY = R;

neighbor = [];
neighbor(size(neighbor,1)+1,:) = [x,y];
for xi=-nX:1:nX
    if x+xi>=1 && x+xi<=sizeX
        for yi=-nY:1:nY
            if y+yi>=1 && y+yi<=sizeY
                if xi^2+yi^2 <= R^2 % 邻域是一个圆
                    if length(find(V(find(V(:,1)==x+xi),2)==y+yi))==0
                        neighbor(size(neighbor,1)+1,:) = [x+xi,y+yi];
                    end
                end
            end
        end
    end
end
