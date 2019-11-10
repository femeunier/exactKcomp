clc ; clear ;

Nprimarycrown = 3 ; % Number of primaries on a crown
Ncrown=4; %number of crowns
Nprimary=Nprimarycrown*Ncrown;


lseg = 0.5; % Root segment length
Lprimary = 50; % Primary length
alpha = linspace(0,(pi/2-pi/2/Ncrown),Ncrown); % 0 --> pi/2, 0 = All vertical, pi/2 = All horizontal
r=0.05; % Primary radius
age=60; % Root System age

beta = linspace(0,(2*pi-(2*pi)/Nprimarycrown),Nprimarycrown); 
Nseg_total = 1;
seg=[];
tip=[];
for j = 1:Ncrown
    for i = 1:Nprimarycrown
    Nseg=Lprimary/lseg;
    
    xtip=Lprimary*sin(alpha(j))*sin(beta(i));
    ytip=Lprimary*sin(alpha(j))*cos(beta(i));
    ztip=-Lprimary*cos(alpha(j));
    
    x1=lseg*sin(alpha(j))*sin(beta(i));
    y1=lseg*sin(alpha(j))*cos(beta(i));
    z1=-lseg*cos(alpha(j));
    
    xseg=linspace(x1,xtip,Nseg);
    yseg=linspace(y1,ytip,Nseg);
    zseg=linspace(z1,ztip,Nseg);
    
    seg_num=(Nseg_total:(Nseg_total+Nseg-1))';
    seg_root=[seg_num,[xseg',yseg',zseg'],[0;seg_num(1:end-1)],ones(Nseg,1),((j-1)*Nprimarycrown+i)*ones(Nseg,1),lseg*ones(Nseg,1),...
        2*pi*r*lseg*ones(Nseg,1),zeros(Nseg,1),linspace(0,age,Nseg)'];
    tip_root=[(j-1)*Nprimarycrown+i,[xseg(end),yseg(end),zseg(end)-0.01],seg_num(end),1,(j-1)*Nprimarycrown+i,Lprimary,0,0,0];
    
    seg=[seg;seg_root];
    tip=[tip;tip_root];
    Nseg_total = Nseg_total+Nseg;    
    end
end

RS=[pwd,'/in/RootSys_simple_fibrousmultialfa'];
write_rootsys(tip,seg,age,1,RS);

figure
draw_root(RS)