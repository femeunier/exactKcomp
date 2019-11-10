clc ; clear ;

Nprimary = 4 ; % Number of primaries
lseg = 0.5; % Root segment length
Lprimary = 50; % Primary length
alpha =pi/4 ; % 0 --> pi/2, 0 = All vertical, pi/2 = All horizontal
r=0.05; % Primary radius
age=60; % Root System age

beta = linspace(0,(2*pi-(2*pi)/Nprimary),Nprimary); 
Nseg_total = 1;
seg=[];
tip=[];
for i = 1:Nprimary
    Nseg=Lprimary/lseg;
    
    xtip=Lprimary*sin(alpha)*sin(beta(i));
    ytip=Lprimary*sin(alpha)*cos(beta(i));
    ztip=-Lprimary*cos(alpha);
    
    x1=lseg*sin(alpha)*sin(beta(i));
    y1=lseg*sin(alpha)*cos(beta(i));
    z1=-lseg*cos(alpha);
    
    xseg=linspace(x1,xtip,Nseg);
    yseg=linspace(y1,ytip,Nseg);
    zseg=linspace(z1,ztip,Nseg);
    
    seg_num=(Nseg_total:(Nseg_total+Nseg-1))';
    seg_root=[seg_num,[xseg',yseg',zseg'],[0;seg_num(1:end-1)],ones(Nseg,1),i*ones(Nseg,1),lseg*ones(Nseg,1),...
        2*pi*r*lseg*ones(Nseg,1),zeros(Nseg,1),linspace(0,age,Nseg)'];
    tip_root=[i,[xseg(end),yseg(end),zseg(end)-0.01],seg_num(end),1,i,Lprimary,0,0,0];
    
    seg=[seg;seg_root];
    tip=[tip;tip_root];
    Nseg_total = Nseg_total+Nseg;    
end

RS=[pwd,'/in/RootSys_simple_fibrous'];
write_rootsys(tip,seg,age,1,RS);

figure
draw_root(RS)