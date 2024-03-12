clear; close all; clc;
run parameters
close all
clc
%%

h_c = 7.94*10^(-8);
h_m = 6.30957344480193e-8;
ic_fum_m = 0.00047226;

ic = [ic_glc_c, ic_pyr_c, ic_lac_c, ic_nadh_c, ic_atp_c, ic_adp_c, h_c,...
    ic_pyr_m,ic_nadh_m, ic_acoa_m, h_m, ic_cit_m, ic_akg_m, ic_atp_m,...
    ic_adp_m, ic_suc_m,ic_fum_m, ic_mal_m,ic_oaa_m, p.qh2_m, p.cred_i];


Tspan = 0:0.01:5; 

options = odeset('RelTol',5e-14,'AbsTol',5e-14);


o2 = 2.6e-5;
o2s = 1 * [1 , 0.7];

for i = 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Integrate   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [t, U] = ode15s(@(t,x) test_2(t,x,p,o2s(i)),Tspan,ic,options);
    
    glc_c = U(:,1);
    pyr_c = U(:,2);
    lac_c = U(:,3);
    nadh_c= U(:,4);
    atp_c = U(:,5);
    adp_c = U(:,6);
    h_c   = U(:,7);
    pH_c = -log10(U(:,7));
    
    pyr_m = U(:,8);
    nadh_m= U(:,9);
    acoa_m= U(:,10);
    h_m   = U(:,11);
    cit_m = U(:,12);
    akg_m = U(:,13);
    atp_m = U(:,14);
    adp_m = U(:,15);
    suc_m = U(:,16);
    fum_m = U(:,17);
    mal_m = U(:,18);
    oaa_m = U(:,19);
    qh2_m = U(:,20);
    cred_i= U(:,21);

    nad_c = p.nadtot_c - nadh_c;
    nad_m = p.nadtot_m - nadh_m;
    rsn_m = nad_m./nadh_m;
    rpn_c = nadh_c./nad_c;
    psn_m = adp_m/atp_m;
    q_m = p.Qtot - qh2_m;
    cox_i = p.Ctot - cred_i;
    
    
    %% Metabolites Plot
    
    figure(1)
    subplot(2,2,1)
    ax = gca;
    plot(t,pH_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('pH_c')
    xlabel('Time (min)')

    subplot(2,2,2)
    ax = gca;
    plot(t,nadh_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[NADH]_c')
    xlabel('Time (min)')

    subplot(2,2,3)
    ax = gca;
    plot(t,lac_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Lactate]_c')
    xlabel('Time (min)')

    subplot(2,2,4)
    ax = gca;
    plot(t,pyr_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Pyruvate]_c')
    xlabel('Time (min)')
    
    
    %% cytosolic Fluxes Plot
    
    gGLUT = 1; % Normalisation Constant. 
    Jglut = glut1(gGLUT,p);

    gGLYC = 13835058055282163712/121280368632803425;
    Jglyc = glyc(gGLYC,glc_c ,nad_c,adp_c,pyr_c,atp_c,nadh_c,p,h_c);


    gMCT4 = 648518346341351424/7243212075608875;
    Jmct4 = mct4(gMCT4,lac_c,p);


    gLDH = 5188146770730811392/78746425648772125;
    Jldh = ldh(gLDH,pyr_c,nadh_c,nad_c,lac_c,p,h_c);


    gMPC = 32997/49625; 
    Jmpc = mpc(gMPC,pyr_c, p,h_c) * o2s(i);



    gATP=432345564227567616/144152400490079875;
    Jatp = atp(gATP,atp_c,adp_c,p,h_c);

    gNHE1 = 324259173170675712/810735649012457375;
    Jnhe1c = nhe1c(gNHE1,h_c,p);

     
    gMDH = 102905318/84375;
    Jmdh = mdh(gMDH,mal_m,rsn_m,rpn_c,p);
    
    figure(2)
    subplot(2,2,1)
    ax = gca;
    plot(t,Jmpc,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('MPC')
    xlabel('Time (min)')

    subplot(2,2,2)
    ax = gca;
    plot(t,Jldh,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('LDH')
    xlabel('Time (min)')

    subplot(2,2,3)
    ax = gca;
    plot(t,2 * Jglyc,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('Glyc')
    xlabel('Time (min)')
    
    subplot(2,2,4)
    ax = gca;
    plot(t,Jmdh ,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('MDH')
    xlabel('Time (min)')
    
    
    %% Mitochondrial Fluxes Plot
    
    figure(3)
    subplot(2,2,1)
    ax = gca;
    plot(t,pyr_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Pyruvate]_m')
    xlabel('Time (min)')

    subplot(2,2,2)
    ax = gca;
    plot(t,nadh_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[NADH]_m')
    xlabel('Time (min)')
    
    subplot(2,2,3)
    ax = gca;
    plot(t,suc_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Succinate]_m')
    xlabel('Time (min)')
    
    subplot(2,2,4)
    ax = gca;
    plot(t, fum_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Fumarate]_m')
    xlabel('Time (min)')
    
end