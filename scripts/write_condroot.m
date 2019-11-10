function write_condroot(file,K,name)
% Write a new condroot named(file+name) with the new K strucuture
% parameters

if nargin <3
   name='mod'; 
end


fidin=fopen(file,'r');
if fidin ==-1
    error (['Trouble while opening ',file]);
end

fidout=fopen([file(1:end-3),'_',name,'.in'],'w');
lines=[7:9,14:16];
line=1;
compt=1;
while ~feof(fidin)
   s = fgetl(fidin);
   if mean(line == lines)>0
       Kalt=[];
       for i=1:size(K{compt},1)
          Kalt=[Kalt,K{compt}(i,:)]; 
       end
       fprintf(fidout,'%s\r\n',num2str(Kalt));
       compt=compt+1;
   else
       fprintf(fidout,'%s\r\n',s);
   end
   line=line+1;
end

fclose(fidin);
fclose(fidout);

%delete(file);
%copyfile([file,'mod'],file);
%delete([file,'mod']);

end