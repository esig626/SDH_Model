clear; close all; clc;
run parameters
close all
clc
%%

h_c = 7.94*10^(-8);
%h_m = 6.30957344480193e-8;
h_m = 1.58*10^(-8);
ic_fum_m = 0.00047226;

ic = [ic_glc_c, ic_pyr_c, ic_lac_c, ic_nadh_c, ic_atp_c, ic_adp_c, h_c,...
    ic_pyr_m,ic_nadh_m, ic_acoa_m, h_m, ic_cit_m, ic_akg_m, ic_atp_m,...
    ic_adp_m, ic_suc_m,ic_fum_m, ic_mal_m,ic_oaa_m, p.qh2_m, p.cred_i,160];


Tspan = 0:0.01:5; 

options = odeset('RelTol',5e-14,'AbsTol',5e-14);


o2 = 2.6e-5;

o2ss = linspace(1,0.005,500);

% %% SS test
% SS = zeros(22,length(o2ss));
% 
% for i = 1:length(o2ss)
% 
%     [t, U] = ode15s(@(t,x) test_2(t,x,p,o2ss(i)),Tspan,ic,options);
% 
%     SS(:,i) = U(end,:);
%     i
% end
% 
% 
% %% Experiment 1. Hypoxia - Benchmarking Model Dynamics.
% close all
% 
%     Ox = o2 * o2ss./o2;
%     glc_c = SS(1,:);
%     pyr_c = SS(2,:);
%     lac_c = SS(3,:);
%     nadh_c= SS(4,:);
%     atp_c = SS(5,:);
%     adp_c = SS(6,:);
%     h_c   = SS(7,:);
%     pyr_m = SS(8,:)/SS(8,1);
%     nadh_m= SS(9,:)/SS(9,1);
%     acoa_m= SS(10,:);
%     h_m   = SS(11,:);
%     cit_m = SS(12,:);
%     akg_m = SS(13,:);
%     atp_m = SS(14,:);
%     adp_m = SS(15,:);
%     suc_m = SS(16,:)/SS(16,end);
%     fum_m = SS(17,:)/SS(17,1);
%     mal_m = SS(18,:);
%     oaa_m = SS(19,:);
%     qh2_m = SS(20,:);
%     cred_i= SS(21,:);
%     dPsi_m= SS(22,:);
%     pH_c = -log10(h_c);
%     pH_m = -log10(h_m);
%     Psi_x = -0.65*dPsi_m;
%     Psi_i = 0.35*dPsi_m;
%  	SO = 0;
%     H2O2 = 0; % Hydrogen Peroxide 
%     
%     
%     nad_c = p.nadtot_c - nadh_c;
%     nad_m = p.nadtot_m - nadh_m;
%     rsn_m = nad_m./nadh_m;
%     rpn_c = nadh_c./nad_c;
%     psn_m = adp_m/atp_m;
%     q_m = p.Qtot - qh2_m;
%     cox_i = p.Ctot - cred_i;
%     pi_m = p.pi_m;
%     k_m = 0.14;
%     
%     
%     
%     
%     
%     figure(1)
%     subplot(2,2,1)
%     ax = gca;
%     plot(Ox,pyr_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[Pyruvate]_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,2)
%     ax = gca;
%     plot(Ox,nadh_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[NADH]_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,3)
%     ax = gca;
%     plot(Ox,suc_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[Succinate]_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,4)
%     ax = gca;
%     plot(Ox,fum_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[Fumarate]_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     figure(2)
%     subplot(2,2,1)
%     ax = gca;
%     plot(Ox,dPsi_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('\Delta \Psi_{m}')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,2)
%     ax = gca;
%     plot(Ox,Psi_x,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('\Psi_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,3)
%     ax = gca;
%     plot(Ox,Psi_i,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('\Psi_i')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,4)
%     ax = gca;
%     plot(Ox, pH_m,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('pH_m')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     
%     figure(3)
%     subplot(2,2,1)
%     ax = gca;
%     plot(Ox,pH_c,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('pH_c')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,2)
%     ax = gca;
%     plot(Ox,nadh_c,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[NADH]_c')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,3)
%     ax = gca;
%     plot(Ox,lac_c,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[Lactate]_c')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')
%     
%     subplot(2,2,4)
%     ax = gca;
%     plot(Ox,pyr_c,'LineWidth',2)
%     hold on
%     ax.LineWidth = 2;
%     ax.FontSize = 19;
%     ylabel('[Pyruvate]_c')
%     xlabel('O_2')
%     set(gca, 'XDir', 'reverse')    
%     
%     
%     
    
    %% SDH
    

% Obtain the IronSulphur Cluster's Fluxes and demonstrate reduction as a
% means of visualising SDH-B silencing. I've set them as AA -- AN. Now Need
% to look at transition diagram ad multiply them by E0 --- E4.

% Once this is clear. Then need to implement volume model, and copy the
% dinamic proposed in the metformin model. 

close all
te = [1, 0.5];
tee = 0 * [0, 8e-11];
bb = [1, 0.7];

for i = 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Integrate   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    o2 = 2.6e-5 * te(i);
    p.atpenin = 1 * tee(i);
    [t, U] = ode15s(@(t,x) test_2(t,x,p,te(i),bb(i)),Tspan,ic,options);
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
    dPsi_m= U(:,22);
    pH_m = -log10(U(:,11));
    nad_c = p.nadtot_c - nadh_c;
    nad_m = p.nadtot_m - nadh_m;
    rsn_m = nad_m./nadh_m;
    rpn_c = nadh_c./nad_c;
    psn_m = adp_m/atp_m;
    q_m = p.Qtot - qh2_m;
    cox_i = p.Ctot - cred_i;
    Psi_x = -0.65*dPsi_m;
    Psi_i = 0.35*dPsi_m;
    SO = 0;
    H2O2 = 0;
    pi_m = p.pi_m;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Cytosolic Concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1)
    subplot(2,3,1)
    ax = gca;
    plot(t,pyr_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Pyruvate]_c')
    xlabel('\tau (ND Time)')
    
    subplot(2,3,2)
    ax = gca;
    plot(t,lac_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Lactate]_c')
    xlabel('\tau (ND Time)')
    
    subplot(2,3,3)
    ax = gca;
    plot(t,pH_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('pH_c')
    xlabel('\tau (ND Time)')
    
    subplot(2,3,4)
    ax = gca;
    plot(t,nadh_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Glucose]_c')
    xlabel('\tau (ND Time)')
    
    subplot(2,3,5)
    ax = gca;
    plot(t,atp_c,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[atp]_c')
    xlabel('\tau (ND Time)')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Mitochondrial Concentrations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	figure(2)
    subplot(2,4,1)
    ax = gca;
    plot(t,pyr_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Pyruvate]_m')
    xlabel('\tau (ND Time)')

    subplot(2,4,2)
    ax = gca;
    plot(t,nadh_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[NADH]_m')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,3)
    ax = gca;
    plot(t,suc_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Succinate]_m')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,4)
    ax = gca;
    plot(t,Psi_x,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('\Psi_m')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,5)
    ax = gca;
    plot(t,Psi_i,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('\Psi_i')
    xlabel('\tau (ND Time)')
   	
    
    subplot(2,4,6)
    ax = gca;
    plot(t,dPsi_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('\Delta \Psi')
    xlabel('\tau (ND Time)')
    
    % pH_c + 0.7
    subplot(2,4,7)
    ax = gca;
    plot(t, pH_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('pH_m')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,8)
    ax = gca;
    plot(t, atp_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[atp]_m')
    xlabel('\tau (ND Time)')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ETC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dG_H = dgh(p,dPsi_m,h_m);
    
    x_C1 = 153122387330596864/8276219508788375;
    J_C1 = etc1(x_C1,nadh_m,nad_m,dG_H,q_m,qh2_m,h_m,p);
    
    x_C3 = 46393490751890289084270077739532288/5651539874953727234232488159903;
    J_C3 = etc3(x_C3,h_m,qh2_m,q_m,pi_m,dG_H,dPsi_m,cox_i,cred_i,p);

    x_C4 = 61963218215226261/247546104476283548;
    dG_C4op = p.dG_C4o - 2 .* p.RT * log(h_m ./ 1e-7)- p.RT./2.*log(o2./1);
    OX_4 =  1 ./ (1 + p.k_O2./o2 ) .* cred_i ./ p.Ctot;
    E_4 = exp(-( dG_C4op + 2 .* dG_H)./(2 .* p.RT));
    NE_4 = exp( p.F .* dPsi_m ./ p.RT);
    J_C4 = x_C4 .* OX_4 .* (E_4 .* cred_i - cox_i .* NE_4);
    
    x_F1 = 10625680370827264/341789307421930625;
    J_F1 = atpsynth(x_F1,pi_m,adp_m, atp_m, dG_H,p); 
    
    x_ANT = 8619889686787129344/1108999006519302125;
    J_ANT = ant(x_ANT,adp_c,atp_c,Psi_i,atp_m,adp_m,Psi_x,p);
    

    [J1,J2, AA, AB, AC, AD, AE, AF, AG, AH, AI, AJ, AK, AL, AM, AN]= ...
        sdh2(bb(i),SO,H2O2,o2,h_m,q_m,qh2_m,suc_m,fum_m,p);

    % Succinate oxidation rate
    gSDH = 1196268651020288/127126864613234625;
    Jsuc = gSDH * J1;
    %solve(Jakgd == Jsuc,gSDH)

    % QH2 reduction 
    gQH2 =0.1779;
    Jqh2 = gQH2 * J2; 
    
    figure(3)
    subplot(2,4,1)
    ax = gca;
    plot(t,J_C1,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('C_I')
    xlabel('\tau (ND Time)')

    subplot(2,4,2)
    ax = gca;
    plot(t,J_C3,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('C_{III}')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,3)
    ax = gca;
    plot(t,J_C4,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('C_{IV}')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,4)
    ax = gca;
    plot(t,mal_m,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('[Malate]_m')
    xlabel('\tau (ND Time)')
    
    subplot(2,4,5)
    ax = gca;
    plot(t,Jsuc,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('C_{II}')
    xlabel('\tau (ND Time)')
    
    % - 7.25
    subplot(2,4,6)
    ax = gca;
    plot(t,J_F1,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('ATPSynthase')
    xlabel('\tau (ND Time)')
clc


    J_ISC_1 = (AF + AG + AH)/4;
    J_ISC_3 = (AA + AB + AC + AD)/4;
    
    SDH_B = (J_ISC_1 + J_ISC_3) ;
    figure(4)
    subplot(1,3,1)
    ax = gca;
    plot(t,J_ISC_3/J_ISC_3(1:1),'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('FADH + O_2 + (B_1)-> FAD + H^+')
    xlabel('\tau (ND Time)')
    
    subplot(1,3,2)
    ax = gca;
    plot(t,J_ISC_1/J_ISC_1(1:1)* bb(i) /2,'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('SQ^-[4Fe-3S]_{red} + 2H -> QH_2 + [4Fe-3S]_{ox}')
    xlabel('\tau (ND Time)')
    
    subplot(1,3,3)
    ax = gca;
    plot(t,SDH_B/(SDH_B(1:1)/ bb(i)*2),'LineWidth',2)
    hold on
    ax.LineWidth = 2;
    ax.FontSize = 19;
    ylabel('SDH_B')
    xlabel('\tau (ND Time)')

end
    
    
    
    
  