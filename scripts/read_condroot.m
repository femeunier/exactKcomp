function [K,mat]=read_condroot(file)
% read condroot (only conductivities)


fidin=fopen(file,'r');
if fidin ==-1
    error (['Trouble while opening ',file]);
end

lines= [7,8,9,15,16,17];
line=1;
compt=1;
while ~feof(fidin)
   s = fgetl(fidin);
   if mean(line == lines)>0
       temp=str2num(s);
%        temp
       K{compt}=[temp(1:2:end)',temp(2:2:end)'];
       mat(compt)=size(K{compt},1);
       compt=compt+1;
   end
   line=line+1;
end
if isempty(K{compt-1})
    fclose(fidin);
    fidin=fopen(file,'r');
    lines= [7,8,9,14,15,16];
    line=1;
    compt=1;
    while ~feof(fidin)
       s = fgetl(fidin);
       if mean(line == lines)>0
           temp=str2num(s);
    %        temp
           K{compt}=[temp(1:2:end)',temp(2:2:end)'];
           mat(compt)=size(K{compt},1);
           compt=compt+1;
       end
       line=line+1;
    end
    
end

fclose(fidin);

end