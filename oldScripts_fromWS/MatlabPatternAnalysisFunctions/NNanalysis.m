 
%modified April 24, 2017

function [d d1 d2]=NNanalysis(x,y,z)


if size(x,1)<size(x,2)
    x=x';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of average min distance between cells
d=zeros(size(x,1),size(x,1));
for i=1:size(x,1);
    for j=1:size(x,1);
        d(i,j)=norm([x(i)-x(j) y(i)-y(j) z(i)-z(j)]);
    end
end



for i=1:size(x,1);%nearest neighbour distance
    dtemp=sort(d(i,:));
    d1(i)=dtemp(2);
    d2(i)=dtemp(3);
end
d1mean=mean(d1);
end