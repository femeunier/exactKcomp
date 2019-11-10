grass=open('grass.mat');
SUFgrass=grass.SUFupscale;
Krsgrass=grass.Krs;
C6grass=grass.C6upscale;

maize=open('maize.mat');
SUFmaize=maize.SUFupscale;
Krsmaize=maize.Krs;
C6maize=maize.C6upscale;
depths=maize.depths;
ndepthsgrass=length(SUFgrass);
ndepthsmaize=length(SUFmaize);
dz=depths(1)-depths(2);


Tpot=0.5; %cm/d
fracmaize=0.6;
fracgrass=1-fracmaize;

maizeplants=10; %plants/m²
grassplants=10000/25;%plants/m²



Tmaize=Tpot*fracmaize*10000/maizeplants; %cm³/d
Tgrass=Tpot*fracgrass*10000/grassplants; %cm³/d

psisoil=-3100-depths*30;

psigrasseff=SUFgrass'*psisoil(1:ndepthsgrass);
psimaizeeff=SUFmaize'*psisoil(1:ndepthsmaize);

psileafmaize=psimaizeeff-Tmaize/Krsmaize;
psileafgrass=psigrasseff-Tgrass/Krsgrass;


sinkgrass=(C6grass*(psisoil(1:ndepthsgrass)-ones(ndepthsgrass,1)*psigrasseff)...
    +Krsgrass.*SUFgrass.*(psigrasseff-psileafgrass))*grassplants/10000/dz; %100 cm/m/(10^6 cm³/m³)

sinkmaize=(C6maize*(psisoil(1:ndepthsmaize)-ones(ndepthsmaize,1)*psimaizeeff)...
    +Krsmaize.*SUFmaize.*(psimaizeeff-psileafmaize))*maizeplants/10000/dz;





figure (1);
yyaxis left;
plot(depths(1:ndepthsgrass),sinkgrass,depths,sinkmaize,'--')
ylabel('sink term (1/d)')
yyaxis right;
plot(depths(1:ndepthsgrass),SUFgrass,depths,SUFmaize,'--')
ylabel('SUF')
legend('grass','maize');
xlabel('depth cm')

sinkgrass=cat(1,sinkgrass,zeros(ndepthsmaize-ndepthsgrass,1));
sinktotal_exact=sinkgrass+sinkmaize;

SUFgrass=cat(1,SUFgrass,zeros(ndepthsmaize-ndepthsgrass,1));

SUF_approx=fracmaize*SUFmaize+fracgrass*SUFgrass;
psieff_approx=SUF_approx'*psisoil;
Krs_approx=(Krsgrass*grassplants+Krsmaize*maizeplants)/10000;
psileaf_approx=psieff_approx-Tpot/Krs_approx;

sinktotal_approx=(Krs_approx.*SUF_approx.*(psisoil-psieff_approx)+...
    Krs_approx.*SUF_approx.*(psieff_approx-psileaf_approx))/dz;

figure (2);
plot(depths,sinktotal_exact,depths,sinktotal_approx,'--')
ylabel('sink term (1/d)')
legend('exact','approx');
xlabel('depth cm')

