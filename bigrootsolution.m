
ra=0.1;
rs=1;

krs=(ra/3+1/(1/rs+1/(ra/2+1/(1/rs+1/(ra+rs)))))^-1;

suf=zeros(3,1);
suf(1)=1/(krs*rs)-ra/(3*rs);
suf(2)=1/(krs*rs)-1/rs*(ra/3+ra/2*(1-1/(rs*krs)+ra/(3*rs)));
suf(3)=1-suf(1)-suf(2);


%kcomp=zeros(3,1);
%kcomp(1)=1/(rs+ra/2+1/(1/rs+1/(rs+ra)));
%kcomp(2)=1/(rs+1/(1/(ra/2+rs)+1/(ra+rs)));
%kcomp(3)=1/(rs+ra+1/(1/rs+1/(rs+ra/2)));

%sufcomp=zeros(3,3);
%sufcomp(1,1)=1;
%sufcomp(1,2)=-1/(kcomp(1)*rs)*(1-kcomp(1)*(rs+ra/2));
%sufcomp(1,3)= 1/(kcomp(1)*rs)*(1-kcomp(1)*(2*rs+ra/2));

%sufcomp(2,1)=-1/(1+(ra/2+rs)/(ra+rs));
%sufcomp(2,2)=1;
%sufcomp(2,3)=-1/(1+(ra+rs)/(ra/2+rs));

%sufcomp(3,1)=1/(kcomp(3)*rs)*(1-kcomp(3)*(2*rs+ra));
%sufcomp(3,2)=-1/(kcomp(3)*rs)*(1-kcomp(3)*(rs+ra));
%sufcomp(3,3)=1;


kcomp=zeros(3,1);
kcomp(1)=1/(rs+ra/2+1/(1/rs+1/(rs+ra)));
kcomp(2)=1/(rs+1/(1/(ra/2+rs)+1/(ra+rs)));
kcomp(3)=1/(rs+ra+1/(1/rs+1/(rs+ra/2)));

sufcomp=zeros(3,3);
sufcomp(1,1)=1;
sufcomp(2,1)=-1/(kcomp(2)*rs)*(1-kcomp(1)*(rs+ra/2));
sufcomp(3,1)= 1/(kcomp(3)*rs)*(1-kcomp(1)*(2*rs+ra/2));

sufcomp(1,2)=-kcomp(2)/kcomp(1)/(1+(ra/2+rs)/(ra+rs));
sufcomp(2,2)=1;
sufcomp(3,2)=-kcomp(2)/kcomp(3)/(1+(ra+rs)/(ra/2+rs));

sufcomp(1,3)=1/(kcomp(1)*rs)*(1-kcomp(3)*(2*rs+ra));
sufcomp(2,3)=-1/(kcomp(2)*rs)*(1-kcomp(3)*(rs+ra));
sufcomp(3,3)=1;

% Doussan solution


IM=...
    [1,-1,0,-1,0,0;...
0,1,-1,0,-1,0;...
0,0,1,0,0,-1;...
0,0,0,1,0,0;...
0,0,0,0,1,0;...
0,0,0,0,0,1];

IMT=cat(2,[-1;0;0;0;0;0],transpose(IM));

K=diag([3/ra 2/ra 1/ra 1/rs 1/rs 1/rs]);
C=IM*K*IMT;

% Couvreur solution

kcomp_couvreur=kcomp./(1-suf)./suf;
sufcomp_couvreur=sufcomp.*(1-suf)+ones(3).*suf';

% solutions for different psi soil distributions

psicolar=-1;
psis1=zeros(3,1);


%distribution 1
%psis1(1)=-1/(suf(1));
%psis1(2)=0;
%psis1(3)=1/(suf(3));

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

%%%%%%%
%distribution 2
%psis2=zeros(3,1);
%psis2(1)=-1/(suf(1));
%psis2(2)=1/(suf(2)+suf(3));
%psis2(3)=1/(suf(2)+suf(3));

%psismean2=sum(psis2.*suf);

%deltapsi2=psis2-psismean2;

%qtot=krs*(psismean2-psicolar);

%sink2=zeros(3,1);

%sink2=suf*qtot+kcomp.*(sufcomp*deltapsi2);




%psix=-inv(C(1:3,2:4))*(C(1:3,5:7)*psis2+C(1:3,1)*psicolar);
%sinkdoussan2=K(4:6,4:6)*(psis2-psix);
%sinkcouvreur2=suf*qtot+kcomp_couvreur.*suf.*(sufcomp_couvreur*deltapsi2);

%sinkcouvreurapprox2=suf*qtot+krs*suf.*deltapsi2;


figure (1);
title('kx=10, kr=1, big root')
yyaxis left;
plot(-1*(1:3),-suf*qtot+sink1,'b',-1*(1:3),-suf*qtot+sinkcouvreurapprox1,'c-');
xlabel('depth');
ylabel('delta sink')
yyaxis right;
plot(-1*(1:3),psis1,'r-');
ylabel('soil water potential')
legend('exact-big root','Couvreur approx')

figure (2);
title('\fontsize{16} kx=10, kr=1, big root')
yyaxis left;
plot(-1*(1:3),sink1,'b',-1*(1:3),sinkcouvreurapprox1,'c-');
xlabel('\fontsize{16} depth');
ylabel('\fontsize{16} sink')
ylim([-0.2 1.8])
yyaxis right;
plot(-1*(1:3),psis1,'r-');
ylabel('\fontsize{16} soil water potential')
legend('\fontsize{16} exact-big root','\fontsize{16} Couvreur approx')
ax = gca; % current axes
ax.FontSize = 12;

figure (3);
plot(-1*(1:3),suf,'b');
title('\fontsize{16} kx=10, kr=1, big root')
xlabel('depth');
ylabel('SUF')
ylim([0.2 0.6]);
legend('\fontsize{16} exact-big root; Couvreur approx')
ax = gca; % current axes
ax.FontSize = 12;



