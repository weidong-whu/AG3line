function compres=recall(method,aerar,groundpath,manhatanres)


reslist=dir(manhatanres);
truelist=dir(groundpath);
compres=zeros(size(truelist,1),4);

for i=94%1:size(truelist,1)
    namei=truelist(i).name;
    if size(namei,2)<10
        continue;
    end
%     for m=1:size(problempics,1)
%         probname=problempics(m);
%         probname=probname{1:8};
%         probname=probname(1:8);
%         isprob=1;
%         if strcmp(probname(1:8),namei(1:8))
%              isprob=0;
%              break;
%         end
%     end
%     if isprob==0
%         continue;
%     end

    
    reallines=importdata([groundpath namei]);
    
    
    %寻找对应的提取结果
    for j=1:size(reslist,1)
        namej=reslist(j).name;

        if size(namej,2)<8+size(method,2)
            continue;
        end

        if strcmp(namei(1:8),namej(1:8))&&...
                strcmp(namej(9:8+size(method,2)),method)
            extractlines=importdata([manhatanres namej]);
              [precision,recall,fscore,iou]  =recallcal(reallines,extractlines,aerar);  
            compres(i,:)= [precision,recall,fscore,iou];
      
         end

    end
end
