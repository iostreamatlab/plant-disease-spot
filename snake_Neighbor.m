function neighbor = snake_Neighbor(sizeX,sizeY,x,y,V)
%
% �õ�������ĳ��[x,y]��Χ�ĵ�����꣨8������[sizeX,sizeY]��ͼ���С
%
% -----------------------------------------------------------------------------------------------------------

% ����Ĵ�С����x����������չnX�����أ���y����������չny������
R = 3;
nX = R;
nY = R;

neighbor = [];
neighbor(size(neighbor,1)+1,:) = [x,y];
for xi=-nX:1:nX
    if x+xi>=1 && x+xi<=sizeX
        for yi=-nY:1:nY
            if y+yi>=1 && y+yi<=sizeY
                if xi^2+yi^2 <= R^2 % ������һ��Բ
                    if length(find(V(find(V(:,1)==x+xi),2)==y+yi))==0
                        neighbor(size(neighbor,1)+1,:) = [x+xi,y+yi];
                    end
                end
            end
        end
    end
end
