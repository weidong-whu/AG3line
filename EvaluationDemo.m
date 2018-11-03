close all
clear
clc
addpath('funcs');
fscores=zeros(14,6);
precisions=zeros(14,6);
recalls=zeros(14,6);
counter=0;
for i=5%1:8 intersection ratio is set as 0.5, you can test all by circle
    aerar=i/10;
%     comprest=recall('.t',aerar,validatetype);
%     compresmcm=recall('MCM',aerar,validatetype);
%     compresppht=recall('ppht',aerar,validatetype);
%     compreslinlet=recall('linelet',aerar,validatetype);
%     compresed=recall('ed',aerar,validatetype);
%     compreslsd=recall('lsd',aerar,validatetype);
%     continue;
    counter=counter+1;
    compres=recall('.t',aerar);
    compres(compres(:,1)==0,:)=[];
    precisions(counter,1)=mean(compres(:,1));
    recalls(counter,1)=mean(compres(:,2));
    fscores(counter,1)=mean(compres(:,3));
    
    compres=recall('MCM',aerar);
    compres(compres(:,1)==0,:)=[];
    precisions(counter,2)=mean(compres(:,1));
    recalls(counter,2)=mean(compres(:,2));
    fscores(counter,2)=mean(compres(:,3));
    
    compres=recall('ppht',aerar);
    compres(compres(:,1)==0,:)=[];
    precisions(counter,3)=mean(compres(:,1));
    recalls(counter,3)=mean(compres(:,2));
    fscores(counter,3)=mean(compres(:,3));
    
    compres=recall('linelet',aerar);
    compres(compres(:,1)==0,:)=[];
    precisions(counter,4)=mean(compres(:,1));
    recalls(counter,4)=mean(compres(:,2));
    fscores(counter,4)=mean(compres(:,3));
    
    compres=recall('ed',aerar);
    compres(compres(:,1)==0,:)=[];
    precisions(counter,5)=mean(compres(:,1));
    recalls(counter,5)=mean(compres(:,2));
    fscores(counter,5)=mean(compres(:,3));
    
    
    
    compresppht=recall('lsd',aerar);
    compresppht(compresppht(:,1)==0,:)=[];
    precisions(counter,6)=mean(compresppht(:,1));
    recalls(counter,6)=mean(compresppht(:,2));
    fscores(counter,6)=mean(compresppht(:,3));
    
end


