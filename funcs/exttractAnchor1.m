function anchor=exttractAnchor1(grad,ang,anchortt,anchorbuffer)

anchorthreshold=1;%角点阈值
[m,n]=size(ang);
anchor=zeros(m,n);

%根据种子点计算
for i=anchorbuffer+1:m-anchorbuffer
    for j=anchorbuffer+1:n-anchorbuffer
        if grad(i,j)<anchorthreshold
            continue;
        end
        isanchor=1;
        isanchorgrad=1;
        if is_go_horizon(ang(i,j))
            for ii=i-anchorbuffer:i+anchorbuffer
                if ii==i
                    continue;
                end
                if grad(i,j)-grad(ii,j)<anchortt
                     isanchorgrad=0;
                    break;
                end
            end
            
            for jj=j-anchorbuffer:j+anchorbuffer
                if jj==j
                    continue;
                end
                if abs(angdiff(ang(i,j),ang(i,jj)))>pi/8
                    isanchor=0;
                    break;
                end
            end
        else
            for jj=j-anchorbuffer:j+anchorbuffer
                if jj==j
                    continue;
                end
                if grad(i,j)-grad(i,jj)<anchortt
                    isanchorgrad=0;
                    break;
                end
            end
            
             for ii=i-anchorbuffer:i+anchorbuffer
                if ii==i
                    continue;
                end
                if abs(angdiff(ang(i,j),ang(ii,j)))>pi/8
                    isanchor=0;
                    break;
                end
            end
        end
        if isanchorgrad==1&&isanchor==1
            anchor(i,j)=2;
        elseif isanchorgrad==1
             anchor(i,j)=1;
        end
    end
end