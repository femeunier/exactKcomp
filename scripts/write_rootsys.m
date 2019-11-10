function write_rootsys(info_tip,seg_info,age,ax,file)

fidout=fopen(file,'w+');
fprintf(fidout,'%s\r\n','Time:');
fprintf(fidout,'\t\t%s\r\n',num2str(age));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Number of seeds');
fprintf(fidout,'\t\t%s\r\n',num2str(1));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','ID, X and Y coordinates of the seeds (one per line)');
fprintf(fidout,'\t\t%s\r\n',num2str([1 0 0]));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Root DM, shoot DM, leaf area:');
fprintf(fidout,'\t\t%s\r\n',num2str([0 0 0]));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Average soil strength and solute concentration experienced by root system:');
fprintf(fidout,'\t\t%s\r\n',num2str([0 0]));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Total # of axes:');
fprintf(fidout,'\t\t%s\r\n',num2str(ax));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Total # of branches, including axis(es):');
fprintf(fidout,'\t\t%s\r\n',num2str(size(info_tip,1)));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Total # of segment records:');
fprintf(fidout,'\t\t%s\r\n',num2str(size(seg_info,1)));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','segID#    x          y          z      prev or  br#  length   surface  mass');
fprintf(fidout,'%s\r\n','origination time');

for i=1:size(seg_info,1)
    fprintf(fidout,'\t\t%s\r\n',num2str(seg_info(i,1:end-1)));
    fprintf(fidout,'\t\t%s\r\n',num2str(seg_info(i,end))); 
end

fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','Total # of growing branch tips:');
fprintf(fidout,'\t\t%s\r\n',num2str(size(info_tip,1)));
fprintf(fidout,'%s\r\n','');
fprintf(fidout,'%s\r\n','tipID#    xg          yg          zg      sg.bhd.tp. ord  br#  tot.br.lgth. axs#');
fprintf(fidout,'%s\r\n','overlength  # of estblished points');
fprintf(fidout,'%s\r\n','time of establishing (-->)');


for i=1:size(info_tip,1)
    fprintf(fidout,'\t\t%s\r\n',num2str(info_tip(i,1:end-2)));
    fprintf(fidout,'\t\t%s\r\n',num2str(info_tip(i,end-1:end))); 
end

fclose(fidout);
end