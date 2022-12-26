function A_tot = compute_A_continous_integral(Pest_stages,trap_on)
    %Input =  Parameters for each stage
    % trap_on = 1 (installed), 0 (non installed) 
    %Output =  Continuous matrix A
    
    mort_trap = 0.4;
    
    % We have 9 states
    A_tot = [-Pest_stages(1).growth-Pest_stages(1).death 0 0 0 0 0 0 Pest_stages(8).birth 0;... %Egg (1)
        Pest_stages(1).growth -Pest_stages(2).growth-Pest_stages(2).death 0 0 0 0 0 0 0;... %L1 (2)
        0 Pest_stages(2).growth -Pest_stages(3).growth-Pest_stages(3).death 0 0 0 0 0 0;... %L2 (3)
        0 0 Pest_stages(3).growth -Pest_stages(4).growth-Pest_stages(4).death 0 0 0 0 0;... %L3 (4)
        0 0 0 Pest_stages(4).growth -Pest_stages(5).growth-Pest_stages(5).death 0 0 0 0;...%P (5)
        0 0 0 0 Pest_stages(5).sex_ratio*Pest_stages(6).growth -Pest_stages(6).growth-Pest_stages(6).death-trap_on*mort_trap 0 0 0;...%AM (6)
        0 0 0 0 (1-Pest_stages(6).sex_ratio)*Pest_stages(7).growth 0 -(Pest_stages(7).mate-Pest_stages(7).death) 0 0;%NMF (7)
        0 0 0 0 0 0 Pest_stages(7).mate-Pest_stages(7).death -Pest_stages(8).growth-Pest_stages(8).death 0; %MF (8)
        0 0 0 0 0 trap_on 0 0 0];  %Trap (9)

    
end