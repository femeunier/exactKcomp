function [str,p,m]=maize_ref()

     p = cell(1,1);

    %-----------------------------------------
    % Tap root
    %-----------------------------------------

    % Parameters editable by the GUI [mean, sd]
    p{1}.r = [1.5, 0.2]; % Initial elongation rate (cm/day)
    p{1}.a = [0.1, 0.01]; % Root radius (cm)  
    p{1}.theta = [0., 0.09]; % Insertion angle(rad)
    p{1}.lb = [2, 0.2]; % Length of basal zone (cm) proximal part       
    p{1}.la = [15, 0.15]; % Length of apical zone (cm) distal part                   
    p{1}.ln = [1, 0.1]; % Length between laterals (cm)
    p{1}.nob = [200, 3]; % Maximal number of laterals (1)

    % Additional parameters (predefined values, not altered by GUI)
    p{1}.name = 'Tap root'; % Name of the root type
    p{1}.color = [1, 0, 0]; % Color of the root (rgb)
    type = 1; % Type of tropism (0 plagio-, 1 gravi-, 2 exo-, 3 chemo-/hydrotropism)
    N = 1.5; % Strength of tropism
    sigma = 0.2; % Expected change of root tip heading (rad/cm)
    p{1}.tropism = [type, N, sigma]; % Root tropism
    p{1}.dx = 1; % Spatial resolution along root axis (cm)
    p{1}.rlt = [Inf, 0]; % Maximal root life time [mean, std] (days)
    p{1}.gf = 1; % Type of growth function (1: negative exponential, 2: linear)

    % Preset by GUI 
    p{1}.successor = [2, 0.4 ; 3, 0.6]; % Laterals [type,probability;type,probability;...]
    p{1}.sef = @(x) 1; % Scale root elongation funcntion
    p{1}.sbpf = @(x) 1; % Scale lateral branching probability function
    p{1}.saf = @(x) 1; % Scale lateral branching angle function

    %-----------------------------------------
    % Short Laterals 
    %-----------------------------------------

    % Parameters editable by the GUI [mean, sd]
    p{2}.r = [1., 0.13]; % Initial elongation rate (cm/day)
    p{2}.a = [0.05, 0.005]; % Root radius (cm)  
    p{2}.theta = [1.22173, 0.12]; % Insertion angle (rad)
    p{2}.lb = [1, 0.1]; % Length of basal zone (cm)        
    p{2}.la = [1, 0.1]; % Length of apical zone (cm)                    
    p{2}.ln = [1, 0.1]; % Length between laterals (cm)
    p{2}.nob = [1, 0.1]; % Maximal number of laterals (1)

    % Additional parameters (predefined values)
    p{2}.name = 'Lateral'; % Name of the root type
    p{2}.color = [0, 1, 0]; % Color of the root (rgb)
    type = 0; % Type of tropism (0 plagio-, 1 gravi-, 2 exo-, 3 chemo-/hydrotropism)
    N = 1; % Strength of tropism
    sigma = 0.3; % Expected change of root tip heading (rad/cm)
    p{2}.tropism = [type, N, sigma]; % Root tropism
    p{2}.dx = 1; % Spatial resolution along root axis (cm)
    p{2}.rlt = [Inf, 0]; % Maximal root life time [mean, std] (days)
    p{2}.gf = 1; % Type of growth function (1: negative exponential, 2: linear)

    % Preset by GUI 
    p{2}.successor = [0 0]; % Laterals [type,probability;type,probability;...]
    p{2}.sef = @(x) 1; % Scale root elongation funcntion
    p{2}.sbpf = @(x) 1; % Scale lateral branching probability function
    p{2}.saf = @(x) 1; % Scale lateral branching angle function
    
    %-----------------------------------------
    % Long laterals 
    %-----------------------------------------

    % Parameters editable by the GUI [mean, sd]
    p{3}.r = [1., 0.13]; % Initial elongation rate (cm/day)
    p{3}.a = [0.05, 0.005]; % Root radius (cm)  
    p{3}.theta = [1.22173, 0.2]; % Insertion angle (rad)
    p{3}.lb = [2, 0.2]; % Length of basal zone (cm)        
    p{3}.la = [3, 0.3]; % Length of apical zone (cm)                    
    p{3}.ln = [1, 0.1]; % Length between laterals (cm)
    p{3}.nob = [30, 3]; % Maximal number of laterals (1)

    % Additional parameters (predefined values)
    p{3}.name = 'Long Lateral'; % Name of the root type
    p{3}.color = [0, 1, 0]; % Color of the root (rgb)
    type = 0; % Type of tropism (0 plagio-, 1 gravi-, 2 exo-, 3 chemo-/hydrotropism)
    N = 1; % Strength of tropism
    sigma = 0.3; % Expected change of root tip heading (rad/cm)
    p{3}.tropism = [type, N, sigma]; % Root tropism
    p{3}.dx = 1; % Spatial resolution along root axis (cm)
    p{3}.rlt = [Inf, 0]; % Maximal root life time [mean, std] (days)
    p{3}.gf = 1; % Type of growth function (1: negative exponential, 2: linear)

    % Preset by GUI 
    p{3}.successor = [0 0]; % Laterals [type,probability;type,probability;...]
    p{3}.sef = @(x) 1; % Scale root elongation funcntion
    p{3}.sbpf = @(x) 1; % Scale lateral branching probability function
    p{3}.saf = @(x) 1; % Scale lateral branching angle function

    %-----------------------------------------
    % Monocot plant data (editable by the GUI)
    %-----------------------------------------

    plantingdepth = 0; % (cm)

%     Basal roots (Seminal)
    basal_first = [1000, 0.]; % First occurence [mean, std] (days)
    basal_delay = [0., 0.]; % Interim time [mean, std] (days)
    basal_max = 0; % Maximal number of basal roots (1)

%     Shoot borne roots
    sb_first = [1000, 1.8]; % Emergence time of first shoot borne root (days)
    sb_delay = [1.5, 0.15]; % Time delay between the emergence of shoot borne roots (days)
    sb_nCR = [2, 0.2]; % Number of shoot borne roots per root crown (1)
    sb_delayRC = [3, 0.3]; % Time delay between the emergence of root crowns (days)
    sb_dzRC = [0.1, 0.01]; % Distance between root crowns along the shoot (cm)

    % No rhizotron
    dfunc=@(x) x(3);
    p = completeParameters(p,dfunc);

    set(0,'RecursionLimit',10000); % sometimes needed
    str = createMonocotRS(plantingdepth,basal_first,basal_delay,basal_max,0.1,...
                    sb_first,sb_delay,sb_nCR,sb_delayRC,sb_dzRC); % create initial string


    m.plantingdepth=plantingdepth;
    m.firstB=basal_first;
    m.delayB=basal_delay;
    m.maxB=basal_max;
    m.nCR=sb_nCR;
    m.firstCR=sb_first;
    m.delayCR=sb_delay;
    m.delayRC=sb_delayRC;
    m.dzRC=sb_dzRC;
end