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


sinkgrass=C6grass*(psisoil(1:ndepthsgrass)-ones(ndepthsgrass,1)*psigrasseff)...
    +Krsgrass.*SUFgrass.*(psigrasseff-psileafgrass);

sinkmaize=C6maize*(psisoil(1:ndepthsmaize)-ones(ndepthsmaize,1)*psimaizeeff)...
    +Krsmaize.*SUFmaize.*(psimaizeeff-psileafmaize);

figure (1);
yyaxis left;
plot(depths(1:ndepthsgrass),sinkgrass,depths,sinkmaize,'g--')
yyaxis right;
plot(depths(1:ndepthsgrass),SUGgrass,depths,SUFmaize,'--')

