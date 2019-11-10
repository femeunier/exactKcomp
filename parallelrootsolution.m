ra=0.1;
rs=1;

krs=(1/(ra+rs)+1/(2*ra+rs)+1/(3*ra+rs));

suf=zeros(3,1);
suf(1)=1/(krs*(rs+ra));
suf(2)=1/(krs*(rs+2*ra));
suf(3)=1/(krs*(rs+3*ra));


%kcomp=zeros(3,1);
%kcomp(1)=1/(ra+rs+1/(1/(2*ra+rs)+1/(3*ra+rs)));
%kcomp(2)=1/(2*ra+rs+1/(1/(ra+rs)+1/(3*ra+rs)));
%kcomp(3)=1/(3*ra+rs+1/(1/(ra+rs)+1/(2*ra+rs)))

%sufcomp=zeros(3,3);
%sufcomp(1,1)=1;
%sufcomp(1,2)=-(1-(rs+ra)*kcomp(1))/(2*ra+rs)/kcomp(1);
%sufcomp(1,3)=-(1-(rs+ra)*kcomp(1))/(3*ra+rs)/kcomp(1);

%sufcomp(2,1)=-(1-(rs+2*ra)*kcomp(2))/(ra+rs)/kcomp(2);
%sufcomp(2,2)=1;
%sufcomp(2,3)=-(1-(rs+2*ra)*kcomp(2))/(3*ra+rs)/kcomp(2);

%sufcomp(3,1)=-(1-(rs+3*ra)*kcomp(3))/(ra+rs)/kcomp(3);
%sufcomp(3,2)=-(1-(rs+3*ra)*kcomp(3))/(2*ra+rs)/kcomp(3);
%sufcomp(3,3)=1;


kcomp=zeros(3,1);
kcomp(1)=1/(ra+rs+1/(1/(2*ra+rs)+1/(3*ra+rs)));
kcomp(2)=1/(2*ra+rs+1/(1/(ra+rs)+1/(3*ra+rs)));
kcomp(3)=1/(3*ra+rs+1/(1/(ra+rs)+1/(2*ra+rs)))

sufcomp=zeros(3,3);
sufcomp(1,1)=1;
sufcomp(2,1)=-(1-(rs+ra)*kcomp(1))/(2*ra+rs)/kcomp(2);
sufcomp(3,1)=-(1-(rs+ra)*kcomp(1))/(3*ra+rs)/kcomp(3);

sufcomp(1,2)=-(1-(rs+2*ra)*kcomp(2))/(ra+rs)/kcomp(1);
sufcomp(2,2)=1;
sufcomp(3,2)=-(1-(rs+2*ra)*kcomp(2))/(3*ra+rs)/kcomp(3);

sufcomp(1,3)=-(1-(rs+3*ra)*kcomp(3))/(ra+rs)/kcomp(1);
sufcomp(2,3)=-(1-(rs+3*ra)*kcomp(3))/(2*ra+rs)/kcomp(2);
sufcomp(3,3)=1;

% Doussan solution

IM=...
    [1,0,0,-1,0,0;...
0,1,0,0,-1,0;...
0,0,1,0,0,-1;...
0,0,0,1,0,0;...
0,0,0,0,1,0;...
0,0,0,0,0,1];

IMT=cat(2,[-1;-1;-1;0;0;0],transpose(IM));

K=diag([1/ra 1/(2*ra) 1/(3*ra) 1/rs 1/rs 1/rs]);
C=IM*K*IMT;


%couvreur solution

kcomp_couvreur=kcomp./(1-suf)./suf;
sufcomp_couvreur=sufcomp.*(1-suf)+ones(3).*suf';

% psi soil distribution 1

psicolar=-1;
psis1=zeros(3,1);

%psis(1)=-1/(suf(1));
%psis(2)=0;
%psis(3)=1/(suf(3));

psis1(1)=-1.0;
psis1(2)=0.0;
psis1(3)=1.0;

psismean1=sum(psis1.*suf);

deltapsi1=psis1-psismean1;

qtot=krs*(psismean1-psicolar);

sink1=zeros(3,1);

sink1=suf*qtot+kcomp.*(sufcomp*deltapsi1);




psix=-inv(C(1:3,2:4))*(C(1:3,5:7)*psis1+C(1:3,1)*psicolar);


sinkdoussan1=K(4:6,4:6)*(psis1-psix);
sinkcouvreur1=suf*qtot+kcomp_couvreur.*suf.*(sufcomp_couvreur*deltapsi1);

sinkcouvreurapprox1=suf*qtot+krs*suf.*deltapsi1;


%%%%%%%%%
%%%% psi soil distribution 2

%psicolar=-1;
%psis=zeros(3,1);

%psis(1)=-1/(suf(1));
%psis(2)=0;
%psis(3)=1/(suf(3));

%psismean2=sum(psis2.*suf);

%deltapsi2=psis2-psismean2;

%qtot2=krs*(psismean2-psicolar);

%sink2=zeros(3,1);

%sink2=suf*qtot2+kcomp.*(sufcomp*deltapsi2);




%psix=-inv(C(1:3,2:4))*(C(1:3,5:7)*psis2+C(1:3,1)*psicolar);


%sinkdoussan2=K(4:6,4:6)*(psis2-psix);
%sinkcouvreur2=suf*qtot2+kcomp_couvreur.*suf.*(sufcomp_couvreur*deltapsi2);

%sinkcouvreurapprox2=suf*qtot2+krs*suf.*deltapsi2;


figure (1);
title('\fontsize{16} kx=10, kr=1, parallel root')
yyaxis left;
plot(-1*(1:3),-suf*qtot+sink1,'b',-1*(1:3),-suf*qtot+sinkcouvreurapprox1,'+c');
xlabel('\fontsize{16} depth');
ylabel('\fontsize{16} delta sink')
yyaxis right;
plot(-1*(1:3),psis1,'r-');
ylabel('\fontsize{16} soil water potential')
legend('\fontsize{16} exact-parallel root','Couvreur approx')
ax = gca; % current axes
ax.FontSize = 12;

figure (2);
title('\fontsize{16} kx=10, kr=1, parallel root')
yyaxis left;
plot(-1*(1:3),sink1,'b',-1*(1:3),sinkcouvreurapprox1,'+c','MarkerSize',10);
xlabel('\fontsize{16} depth');
ylabel('\fontsize{16} sink')
ylim([-0.2 1.8])
yyaxis right;
plot(-1*(1:3),psis1,'r-');
ylabel('\fontsize{16} soil water potential')
legend('\fontsize{16} exact-parallel root','\fontsize{16} Couvreur approx')
ax = gca; % current axes
ax.FontSize = 12;

figure (3);
plot(-1*(1:3),suf,'b');
title('\fontsize{16} kx=10, kr=1, parallel root')
xlabel('\fontsize{16} depth');
ylim([0.2 0.6]);
ylabel('\fontsize{16} SUF');
legend('\fontsize{16} exact-parallel root','\fontsize{16} Couvreur approx')
ax = gca; % current axes
ax.FontSize = 12;
