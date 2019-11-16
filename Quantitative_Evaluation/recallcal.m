function [precision,recall,fscore,iou,precision1,recall1] =recallcal(trueline,extractline,area)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明


cc= extractline(:,1)==extractline(:,3)&extractline(:,2)==extractline(:,4);
extractline(cc,:)=[];

lengthsum=0;
lengthtrue=0;
counter=0;
iouunion=0;
for i=1:size(trueline,1)

    tl1=single(trueline(i,1:2));
    tl2=single(trueline(i,3:4));
%     close all
%  if norm(tl1-tl2)<60
%      continue;
%  end
%     plot([tl1(1) tl2(1)],[tl1(2) tl2(2)],'color','k','linewidth',2);
%     hold on
    lengthtrue= lengthtrue+norm(tl2-tl1);
%     
%             plot([ tl1(1),tl2(1)],640-[tl1(2),tl2(2)],...
%                 'Color','blue','LineWidth',4);
%             hold on
    
    v1=[tl1 0];
    v2=[tl2 0];
    pts1=[extractline(:,1:2) zeros(size(extractline,1),1)];
    
    pts2=[extractline(:,3:4) zeros(size(extractline,1),1)];
%     distance2=point_to_line_distance(pt2, v1, v2);
    distance1=point_to_line_distance((pts1+pts2)./2, v1, v2);

    %距离判断 中心点小于一个像素
    idx=find(distance1<=1.5);
    if isempty(idx)
        continue;
    end
   
    lines=extractline(idx,:);
    lentl=norm(tl1-tl2);
   
    for j=1:size(lines,1)
        %线段向量的夹角要小于5度
        u=[lines(j,3:4)-lines(j,1:2) 0];
        v=[tl2-tl1 0];
        
        pt1=lines(j,1:2);
        pt2=lines(j,3:4);
        
        distance1=point_to_line_distance([(pt1+pt2)./2 0], v1, v2);
        distance2=point_to_line_distance((v1+v2)./2, [pt1 0], [pt2 0]);
        if distance1>1.001&&distance2>1.0001
            continue
        end
        %夹角判断
        cosang2line=abs(dot(u,v)/(norm(u)*norm(v)));
        if  cosang2line<cos(pi/36)
            continue;
        end
        
        
 
        lenext=norm(pt1-pt2);
        %长度应该大于11
    
        if norm(pt1-pt2)<3||lenext/lentl<area%||lenext/lentl>1/area
%             continue;
        end
        
        %算出投影
        pt1= ProjPoint( pt1,[tl1 tl2] );
        pt2=ProjPoint( pt2,[tl1 tl2] );
        dottl_p1=(dot(tl1-pt1,tl2-pt1)<=0);%pt1 是否在tl内部
        dottl_p2=(dot(tl1-pt2,tl2-pt2)<=0);%pt2 是否在tl内部
        dotp_tl1=dot(pt1-tl1,pt2-tl1)<=0;%tl1 是否在pt内部
        dotp_tl2=dot(pt1-tl2,pt2-tl2)<=0;%tl2 是否在pt内部
        %1 完全包含groundtruth
        if dotp_tl1&&dotp_tl2
            intersec=lentl;
%              plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
%              'total in1'
        elseif  dottl_p1&& dottl_p2%被gd完全包含
            intersec=lenext;
%              plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
%              'total in2'
        elseif dottl_p1&&~dottl_p2&&dotp_tl1%pt1 内 pt2 外 tl1内
            intersec=norm(pt1-tl1);
         
% %              plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
% %              'cross'
        elseif dottl_p1&&~dottl_p2&&dotp_tl2%pt1 内 pt2 外 tl2内
            intersec=norm(pt1-tl2);
%                 plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
%              'cross'
        elseif ~dottl_p1&&dottl_p2&&dotp_tl2%pt2 内 pt1 外 tl2内
            intersec=norm(pt2-tl2);
%                 plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
%              'cross'
        elseif ~dottl_p1&&dottl_p2&&dotp_tl1%pt2 内 pt1 外 tl1内
            intersec=norm(pt2-tl1);
%                 plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','g','linewidth',2);
%              'cross'
        else
%              plot([pt1(1) pt2(1)],[pt1(2) pt2(2)],'color','r','linewidth',2);
%              j
            continue;
        end
        if intersec/lentl<area||intersec/lenext<area
            continue;
        end
        iouunion=iouunion+lentl+lenext-intersec;
        lengthsum=lengthsum+intersec;
        counter=counter+1;
    end
end
      
totallength=sum(sqrt((extractline(:,1)-extractline(:,3)).^2+...
    (extractline(:,2)-extractline(:,4)).^2));
precision=lengthsum/totallength;
recall=lengthsum/lengthtrue;
fscore=2*(precision*recall)/(precision+recall);
iou=lengthsum/(totallength+lengthtrue-lengthsum);
precision1=counter/size(extractline,1);
recall1=counter/size(trueline,1);

end

%找到垂点的发射点
function [error,errorl]=findprjPt(pt1,pt2,ptend,pp)
    %反投回去
    distance1=point_to_line_distance([ptend 0], [pt1 0], [pt2 0]);
    distance2=point_to_line_distance([pp 0], [pt1 0], [pt2 0]);
    error=(distance1+distance2)/2;
    errorl=error*norm(ptend-pp);
end
function distance=point_to_line_distance(pt, v1, v2)
%Calculate distance between a point and a line in 2D or 3D.
% syntax:
% distance = point_to_line(pt, v1, v2)
% pt is a nx3 matrix with xyz coordinates for n points
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
%
% 2D input is extended to 3D by setting all z-values to 0.
%
% The actual calculation is a slightly edit version of this line:
% distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)
% (this line only works for a single 3D point)
%
% Example input:
% v1 = [0,0,0];
% v2 = [3,0,0];
% pt = [0,5,0;0,10,0;0,15,0];
% distance = point_to_line_distance(pt, v1, v2);
%
% Compatibility:
% Matlab: should work on all releases (tested on R2017b, R2012b and R6.5)
% Octave: tested on 4.2.1
% OS:     should work cross-platform
%
% Version: 1.2
% Date:    2017-12-29
% Author:  H.J. Wisselink
% Email=  'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})


%parse inputs (changed to be compatible with R6.5)
if nargin~=3
    error('Incorrect number of inputs, expected 3.');
end
if ~isnumeric(pt) || ~any(size(pt,2)==[2 3]) || any(size(pt)==0)
    error('First input (pt) is not numeric or has an incorrect shape.')
end
if ~isnumeric(v1) || numel(v1)~=size(pt,2)
    error(['Second input (v1) is not numeric or has an incorrect size.' char(5+5),...
        'Expected 1x3 or 1x2, which should match the first input.'])
end
if ~isnumeric(v2) || numel(v2)~=size(pt,2)
    error(['Third input (v2) is not numeric or has an incorrect size.' char(5+5),...
        'Expected 1x3 or 1x2, which should match the first input.'])
end

%prepare inputs
v1=v1(:)';%force 1x3 or 1x2
if length(v1)==2,v1(3)=0;end%extend 1x2 to 1x3 if needed
v2=v2(:)';%force 1x3 or 1x2
if length(v2)==2,v2(3)=0;end%extend 1x2 to 1x3 if needed
v1_ = repmat(v1,size(pt,1),1);
v2_ = repmat(v2,size(pt,1),1);

%actual calculation
a = v1_ - v2_;
b = pt - v2_;
distance = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
%this is equivalent to the following line for a single point
%distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)
end

function proj_point = ProjPoint( point,line_p )
x1 = line_p(1);
y1 = line_p(2);
x2 = line_p(3);
y2 = line_p(4);

x3 = point(1);
y3 = point(2);

yk = ((x3-x2)*(x1-x2)*(y1-y2) + y3*(y1-y2)^2 + y2*(x1-x2)^2) / (norm([x1-x2,y1-y2])^2);
xk = ((x1-x2)*x2*(y1-y2) + (x1-x2)*(x1-x2)*(yk-y2)) / ((x1-x2)*(y1-y2));


if x1 == x2
    xk = x1;
end

if y1 == y2
    xk = x3;
end

proj_point = [xk,yk];

end
function Pr=getSpPoint(Point,Line)
    % getSpPoint(): find Perpendicular on a line segment from a given point
    x1=Line(1);
    y1=Line(2);
    x2=Line(3);
    y2=Line(4);
    x3=Point(1);
    y3=Point(2);

    px = x2-x1;
    py = y2-y1;
    dAB = px*px + py*py;

    u = ((x3 - x1) * px + (y3 - y1) * py) / dAB;
    x = x1 + u * px;
    y = y1 + u * py;

    Pr=[x,y];

end

% vv=int32(abs(pt1-pt2));
%     if vv(1)>vv(2)
%         %以x方向构造linspace
%         xx=linspace(pt1(1),pt2(1),vv(1));
%         yy=linspace(pt1(2),pt2(2),vv(1));
%         %寻找x对应的y
%         y=yy(xx<=ptend(1)+1|xx>=ptend(1)-1);
%         y=y(1);
%         x=ptend(1);
%         %找到以后 构造yy
%         lv=norm(pp-ptend);
%         yy1=linspace(pt1(2),y,lv);
%         yy2=linspace(pp(2),ptend(2),lv);
%         error=sum(abs(yy1-yy2))/lv;
%     else
%         %以y方向构造linspace
%         xx=linspace(pt1(1),pt2(1),vv(2));
%         yy=linspace(pt1(2),pt2(2),vv(2));
%         %寻找y对应的x
%         x=xx(yy<=ptend(2)+1|yy>=ptend(2)-1);
%         x=x(1);
%         y=ptend(2);
%         
%         %找到以后 构造xx
%         lv=norm(pp-ptend);
%         yy1=linspace(pt1(1),x,lv);
%         yy2=linspace(pp(1),ptend(1),lv);
%         error=sum(abs(yy1-yy2))/lv;
%     end
function [x4,y4]=PointPro2Line(x3,y3,x1,y1,x2,y2)
    k = ((y2-y1) * (x3-x1) - (x2-x1) * (y3-y1)) / ((y2-y1)^2 + (x2-x1)^2);
    x4 = x3 - k * (y2-y1);
    y4 = y3 + k * (x2-x1);
end

function [interSect,flag] = coincide(line1,line2)
%需要先将line1 line2 转换成一维
% line1是一个1*2的矩阵，分别代表线段的两端。line2同上。
% interSect是两个两个线段重合的两个端点，
% 如果两个线段没有重合的点，那么interSect = -1
% flag = -1，两个线段相离
% flag = 0，两个线段有交集
% 以上

% % 两个线段
% a = [2 5];
% b = [3 4];
% % 所有线段的端点，并进行了排序

a = [norm(line1(1:2)) norm(line1(3:4))];
b = [norm(line2(1:2)) norm(line2(3:4))];
c = [a b];
c = sort(c);

flag = -1;
result = -1;
index = 1;
for i = 1:length(c)
    if( ((c(i) >= a(1) ) && ( c(i) <= a(2) ))   &&  ( ( c(i)>=b(1)  ) &&  (c(i)<=b(2))   ) )
        result(index) = c(i);
        index = index+1;
        flag = 0;
    end
end

% for debug
if( flag == -1 )
%     disp 两个线段相离!
    interSect = result;
else
    % 把重复的端点剔除
    interSect = unique(result);
%     disp 重合的线段是：
%     interSect
end
end
