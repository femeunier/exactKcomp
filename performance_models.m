clc ; clear ; 

WD = pwd;

addpath(genpath(WD));

% RootSys = [WD,'/in/RootSys_singleroot'];  % single root
% RootSys = [WD,'/in/RootSys_fibrous'];  % Fibrous root system of Valentin
% RootSys = [WD,'/in/RootSys_taproot'];  % Taproot 
RootSys = [WD,'/in/RootSys_simple_fibrous'];  % Simple fibrous system

% CR=[WD,'/in/CondRoot_singleroot.in'];
% CRtemp=[WD,'/in/CondRoot_singleroot_mod.in'];

CR=[WD,'/in/CondRoot.in'];
CRtemp=[WD,'/in/CondRoot_mod.in'];

K=read_condroot(CR);
alpha_init=K{1}(1,2)/K{4}(1,2);

alpha = logspace(-6,-1,101);

for ialpha=1:length(alpha)
    disp(ialpha/length(alpha))
    
    Ktemp=K;
    for i = 1:3
        Ktemp{i}=[Ktemp{i+3}(:,1),Ktemp{i+3}(:,2)*alpha(ialpha)];
    end
    write_condroot(CR,Ktemp,'mod');
    [DeltaSparallel,DeltaSbigroot] = calc_model_differences(RootSys,CRtemp);
    RMSE_par(ialpha) = sqrt(sum(DeltaSparallel.^2)/(length(DeltaSparallel-1))) ;
    RMSE_BR(ialpha) = sqrt(sum(DeltaSbigroot.^2)/(length(DeltaSbigroot-1))) ;
    
end

figure
plot(alpha,RMSE_par)
hold on
plot(alpha,RMSE_BR)
set(gca,'yscale','log','xscale','log')
legend('Parallel','Big root')
ylabel('RMSE model')
xlabel('k_r/k_x ratio')

