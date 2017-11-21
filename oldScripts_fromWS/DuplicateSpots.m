
%Update May 2, 2017

clear all
close all
filename='D:\Spatial Analysis Nico Laure\temp nico\ImarisObjects';

load([filename,'.mat'])
% 
% N=20;
% 
% 




a=who('-regexp', '1_Position');
b=who('-regexp', '1_Radius');
c=who('-regexp', '1_Edges');
e=who('-regexp', '1_Time');

for spotnameid=1:length(a)
    clear d        
    
    eval(['position=',a{spotnameid},';']);
    eval(['radius=',b{spotnameid},';']);
    eval(['edges=',c{spotnameid},';']);
    eval(['time=',e{spotnameid},';']);
    
    if ~isempty(position)

        
        % position=position(1:10,:);
        % position=[position;position];
        
        Nspots=size(position,1);
        
        x=position(:,1);y=position(:,2);z=position(:,3);
        Duplicata=x.*0;
        id=x'.*0;
        
        
        for i=1:Nspots
            % for i=1:10
            id=id.*0;id(i)=1;
            if Duplicata(i)==0
                for j=1:Nspots;
                    d(j)=norm([x(i)-x(j) y(i)-y(j) z(i)-z(j)]);
                end
                
                if length(find(d==0))>2%in case we have triplets or more
                    find(d==0)
                end
                Duplicata(find(d==0&id==0))=1;
            end
        end
        
        
        
        
        
        
        position=[x(Duplicata==0) y(Duplicata==0) z(Duplicata==0)];
        time=time(Duplicata==0);
        radius=radius(Duplicata==0,:);
        edges=[];
        %display
        a{spotnameid}
        Nduplicata=sum(Duplicata);
        NspotsNEW=size(position,1);
        [Nspots NspotsNEW Nduplicata]% for each spots object, we display the initial number of spots, the final and number of duplicates
        
        
        eval([a{spotnameid},'=position;']);
        eval([b{spotnameid},'=radius;']);
        eval([c{spotnameid},'=edges;']);
        eval([e{spotnameid},'=time;']);
        % vObjectsName=strcat(vObjectsName,' No Duplicates');
    else
        eval([a{spotnameid},'=single([0 0 0])']);
        eval([b{spotnameid},'=single([0 0 0])']);
        eval([e{spotnameid},'=int32(0);']);
    end
end

save([filename,'DuplicateCorrected'], '-regexp', 'v');
