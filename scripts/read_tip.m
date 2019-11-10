function info_tip=read_tip(rootsys,initial_tip,n_tip)

fid=fopen(rootsys,'r');

if fid==-1
   error('reading ',rootsys, 'impossible')
end
C=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',initial_tip);


C1=C{1};
Closest=dsearchn(C1,[1:n_tip]');
for i=2:length(Closest)
    if Closest(i)<=(Closest(i-1))
        Closest(i)=(Closest(i+1)+Closest(i-1))/2;
    end
end

for i=1:n_tip
    for j=1:9
        info_tip(i,j)=C{j}(Closest(i));
    end
    for j=1:2
        info_tip(i,j+9)=C{j}(Closest(i)+1);
    end
end
fclose(fid);


end