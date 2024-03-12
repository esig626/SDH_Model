function [J1,J2]=sdh(SO,H2O2,o2,h_m,q_m,qh2_m,suc_m,fum_m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron Transport Chain (Complex II - SDH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Binding Polynomials
BP_Q= ( 1 + (q_m/p.KD_Q) + (qh2_m / p.KD_QH2)+(p.atpenin/p.KD_atpenin));
BP_FAD=(1 + (suc_m / p.KD_SUC) + (fum_m / p.KD_FUM));    
A_5=(1+p.atpenin/p.KD_atpenin2)/(1+p.beta_atpenin*p.atpenin/p.KD_atpenin2);

% Quinol reduction
kf_QH2 = p.kf_QH2 * A_5;
                    
% Succinate oxidation 
kf_SUC = p.kf_SUC * A_5 / (1 + h_m / p.KD_H);

% My definitions
H_FADH  = (1 + h_m / p.K_FADH);
H_FADH2 = (1 + h_m / p.K_FADH2);
H_iFADH = (1 + p.K_FADH / h_m);

% Flavins 
kf_FADH  = p.kf_FADH / BP_FAD / H_FADH;
kf_FADH2 = p.kf_FADH2 / BP_FAD / H_FADH2;

% Peroxide
kf_H2O2  = p.kf_H2O2 / BP_FAD / H_FADH2;

%% Midpoint potential pH corrections.
% FAD potential as a pH dependence
% One electron donation wrt proton concentration
Em_FAD_FADH = p.Em0.FAD_FADH + p.RTF * log(h_m * H_iFADH);
Em_FADH_FADH2 = p.Em0.FADH_FADH2 + p.RTF * log(h_m * H_FADH2 / H_FADH);

% Two electrons donation wrt proton concentration (stepwise)
Em_FAD_FADH2 = p.Em0.FAD_FADH2 + p.RTF/2 * log(h_m^2 * H_FADH2 * H_iFADH/H_FADH);

% Em of SQ/QH2 couple at a given pH
Em_SQ_QH2 = 2 * p.Em0.Q_QH2 + 2 * p.RTF * log(h_m) - p.Em0.Q_SQ; 

% Direct 2 e- donation
Em_Q_QH2 = p.Em0.Q_QH2 + p.RTF * log(h_m);

% Redox of Fumarate to Succinate                             
Em_Fum_Suc = p.Em0.Fum_Suc + p.RTF * log(h_m);


% Midpoint potentials for electron transfer of bound states

% Electron transfer free energy (mV)
Emb_SQ_QH2 = Em_SQ_QH2 - p.RTF * log(p.KD_QH2);

% Electron transfer free energy (mV)
Emb_Q_SQ = p.Em0.Q_SQ + p.RTF * log(p.KD_Q); 

% Electron transfer free energy (mV)
Emb_Fum_Suc = Em_Fum_Suc - p.RTF/2 * log(p.KD_SUC) + p.RTF/2 * log(p.KD_FUM); 

% Equilibrium Constants: 
Keq_SO_ISC3   = exp((p.Em.O2_SO  - p.Em.ISC3) / p.RTF );
Keq_SO_SQ     = exp((p.Em.O2_SO  - Emb_Q_SQ) / p.RTF );
Keq_SO_FADH   = exp((p.Em.O2_SO - Em_FAD_FADH) / p.RTF );
Keq_SO_FADH2  = exp((p.Em.O2_SO - Em_FADH_FADH2) / p.RTF );
Keq_H2O2_FADH2= exp(2/p.RTF * (p.Em.O2_H2O2 - Em_FAD_FADH2));
Keq_QH2_ISC3  = exp((Emb_SQ_QH2  - p.Em.ISC3) / p.RTF ); % should be bound state Em for SQ_QH2
Keq_FUM_FADH2 = exp(2/p.RTF * (Em_FAD_FADH2 - Emb_Fum_Suc)); % should be Em.FAD - Emb.Suc

% % Since the Em of ISC2 is quite low, the model assumes ISC1 as second e- source
% Keq.GEA_ISC   = exp(2/p.RTF * (p.Em.GEA - (p.Em.ISC3 + p.Em.ISC1))); 


%% Boltzman redox poise potentials
% compute Eh for Q/SQ reaction
Eh_FAD_FADH = Em_FAD_FADH;
Eh_FADH_FADH2 = Em_FADH_FADH2;
Eh_ISC1 = p.Em.ISC1;
Eh_ISC2 = p.Em.ISC2;
Eh_ISC3 = p.Em.ISC3;    
Eh_Q_SQ = Emb_Q_SQ + p.RTF * log(q_m / p.KD_Q / BP_Q);

%% ΔG = -n F * ΔE, (Free Energy - Gibbs, F is Faraday's Constant)
% Here, n is the number of electrons transferred in the reaction.

DfG1=-[Eh_FAD_FADH;Eh_FADH_FADH2;Eh_ISC1;Eh_ISC2;Eh_ISC3;Eh_Q_SQ]/p.RTF;

%% Substate fraction calculations. 
% E0
s0 = 1;
% E1
E1i = DfG1(p.substates.s1,:); 
s1  = exp(-E1i) / sum(exp(-E1i),1);
% E2
E2i = sum(DfG1(p.substates.s2),2);          
s2  = exp(-E2i) / sum(exp(-E2i));
% E3
E3i = sum(DfG1(p.substates.s3),2); 
s3  = exp(-E3i) / sum(exp(-E3i));
% E4
E4i = sum(DfG1(p.substates.s4),2); 
s4  = exp(-E4i) / sum(exp(-E4i),1);

%% State Transitions
% This calculation helps us understand the likelihood of the enzyme 
% being in a specific substate.

% E1 substates used in state transitions
s1_FADH = sum(s1(p.s1.FADH,:),1);
s1_FAD  = sum(s1(p.s1.FADox,:),1);
s1_SQ   = sum(s1(p.s1.SQ,:),1);
s1_SQempty = sum(s1(p.s1.SQempty,:),1);
s1_SQempty_ISC3ox = sum(s1(p.s1.SQempty_ISC3ox,:),1);
% s1_ISC1ox_ISC3ox = sum(s1(p.s1.ISC1ox_ISC3ox,:),1);
s1_ISC3 = sum(s1(p.s1.ISC3,:),1);
%s1_ISC3ox = sum(s1(p.s1.ISC3ox,:),1);

% E2 substates used in state transitions
s2_FAD = sum(s2(p.s2.FADox,:),1);
s2_FADH = sum(s2(p.s2.FADH,:),1);
s2_FADH2 = sum(s2(p.s2.FADH2,:),1);
s2_SQ = sum(s2(p.s2.SQ,:),1);
s2_SQempty = sum(s2(p.s2.SQempty,:),1);
s2_ISC3 = sum(s2(p.s2.ISC3_SQempty,:),1);
%s2_ISC3ox = sum(s2(p.s2.ISC3ox,:),1);
% s2_ISC1_ISC3 = sum(s2(p.s2.ISC1_ISC3,:),1); 
s2_SQ_ISC3 = sum(s2(p.s2.SQ_ISC3,:),1);
s2_SQempty_ISC3ox = sum(s2(p.s2.SQempty_ISC3ox,:),1);
% s2_ISC1ox_ISC3ox = sum(s2(p.s2.ISC1ox_ISC3ox,:),1);

% E3 substates used in state transitions
s3_FAD = sum(s3(p.s3.FADox,:),1);
s3_FADH = sum(s3(p.s3.FADH,:),1);
s3_FADH2 = sum(s3(p.s3.FADH2,:),1);
s3_SQ = sum(s3(p.s3.SQ,:),1);
s3_SQempty = sum(s3(p.s3.SQempty,:),1);
s3_ISC3 = sum(s3(p.s3.ISC3_SQempty,:),1);
%s3_ISC3ox = sum(s3(p.s3.ISC3ox,:),1);
% s3_ISC1_ISC3 = sum(s3(p.s3.ISC1_ISC3,:),1);
s3_SQ_ISC3 = sum(s3(p.s3.SQ_ISC3,:),1);
% s3_ISC1ox_ISC3ox = sum(s3(p.s3.ISC1ox_ISC3ox,:),1);
% s3_SQempty_ISC3ox = sum(s3(p.s2.SQempty_ISC3ox,:),1);

% E4 substates used in state transitions
%s4_FAD    = sum(s4(p.s4.FADox,:),1);
s4_FADH = sum(s4(p.s4.FADH,:),1);
s4_FADH2 = sum(s4(p.s4.FADH2,:),1);
s4_SQ = sum(s4(p.s4.SQ,:),1);
% s4_SQempty  = sum(s4(p.s4.SQempty,:),1);
s4_ISC3 = sum(s4(p.s4.ISC3_SQempty,:),1);
% s4_ISC3ox   = sum(s4(p.s4.ISC3ox,:),1);
% s4_ISC1_ISC3 = sum(s4(p.s4.ISC1_ISC3,:),1);
s4_SQ_ISC3 = sum(s4(p.s4.SQ_ISC3,:),1);
% s4_ISC1ox_ISC3ox = sum(s4(p.s4.ISC1ox_ISC3ox,:),1);
% s4_SQempty_ISC3ox = sum(s4(p.s4.SQempty_ISC3ox,:),1);

%% State transition details and rates

% 1) E0->E1
% SO + Q -> O2 + SQ
% FAD + H + SO -> FADH + O2

k01_SO_ISC3 = (p.kf_SO2 / Keq_SO_ISC3) * (1 /BP_Q) * SO * s0;
k01_SO_Q = (p.kf_SO / Keq_SO_SQ) * (q_m / p.KD_Q / BP_Q) * SO * s0;
k01_SO_FAD = (kf_FADH / Keq_SO_FADH) * SO * s0;
k0_1 = k01_SO_Q + k01_SO_FAD + k01_SO_ISC3;

% 2) E1->E0
% O2 + SQ -> SO + Q
% FADH + O2 -> FAD + H + SO

k10_O2_ISC3 = p.kf_SO2 * o2 * s1_ISC3 / BP_Q;
k10_O2_SQ = p.kf_SO * o2 * s1_SQ;
k10_O2_FADH = kf_FADH * o2 * s1_FADH;
k1_0 = k10_O2_SQ  + k10_O2_FADH + k10_O2_ISC3;

% 3) E0 -> E2                       
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ- + [4Fe_3S]red + 2H

k02_H2O2_FAD = (kf_H2O2./Keq_H2O2_FADH2) * H2O2 * s0;
k02_SUC_FAD =  kf_SUC * suc_m / p.KD_SUC / BP_FAD * s0;
k02_QH2_SQempty_ISC3ox =  (kf_QH2 /  Keq_QH2_ISC3) * qh2_m / p.KD_QH2 / BP_Q  * s0;
k0_2 = k02_H2O2_FAD + k02_SUC_FAD + k02_QH2_SQempty_ISC3ox; 

% 4) E2->E0
% O2 + FADH2 -> H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox

k20_O2_FADH2 = kf_H2O2 * o2 * s2_FADH2;
k20_FUM_FADH2 = (kf_SUC / Keq_FUM_FADH2) * fum_m / p.KD_FUM / BP_FAD * s2_FADH2; 
k20_SQ_ISC3_QH2 =  kf_QH2 * s2_SQ_ISC3;
k2_0 = k20_O2_FADH2 + k20_FUM_FADH2 + k20_SQ_ISC3_QH2; 

% 5) E1->E2 
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k12_SO_ISC3 = (p.kf_SO2./Keq_SO_ISC3) * 1 / BP_Q * SO * s1_SQempty;
k12_SO_Q = (p.kf_SO / Keq_SO_SQ) * q_m / p.KD_Q / BP_Q * SO * s1_SQempty;
k12_SO_FAD = (kf_FADH / Keq_SO_FADH) * SO * s1_FAD;
k12_SO_FADH = (kf_FADH / Keq_SO_FADH2) * SO * s1_FADH;
k1_2 = k12_SO_Q + k12_SO_FAD + k12_SO_FADH + k12_SO_ISC3;

% 6) E2->E1
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H + SO
% O2 + FADH2 -> FADH + H+ SO


k21_O2_ISC3 = p.kf_SO2 * o2 * s2_ISC3 / BP_Q;
k21_O2_SQ = p.kf_SO * o2 * s2_SQ;
k21_O2_FADH = kf_FADH * o2 * s2_FADH;
k21_O2_FADH2 = kf_FADH2 * o2 * s2_FADH2;
k2_1 = k21_O2_SQ  + k21_O2_FADH + k21_O2_FADH2 + k21_O2_ISC3;

% 7) E1->E3
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ + [4Fe_3S]red + 2H

k13_H2O2_FAD = (kf_H2O2 / Keq_H2O2_FADH2) * H2O2 * s1_FAD;
k13_SUC_FAD =  kf_SUC * suc_m / p.KD_SUC / BP_FAD * s1_FAD; 
k13_QH2_SQempty_ISC3ox = (kf_QH2 / Keq_QH2_ISC3) * qh2_m / p.KD_QH2 / BP_Q * s1_SQempty_ISC3ox;
k1_3 = k13_H2O2_FAD + k13_SUC_FAD + k13_QH2_SQempty_ISC3ox; 

% 8) E3->E1
% O2 + FADH2 ->  H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox

k31_O2_FADH2 = kf_H2O2 * o2 * s3_FADH2;
k31_FUM_FADH2 = (kf_SUC / Keq_FUM_FADH2) * fum_m / p.KD_FUM /BP_FAD * s3_FADH2;
k31_SQ_ISC3_QH2 =  kf_QH2 * s3_SQ_ISC3;
k3_1 = k31_O2_FADH2 + k31_FUM_FADH2 + k31_SQ_ISC3_QH2;

% 9) E2->E3
% SO + ISC3ox -> O2 + ISC3red
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k23_SO_ISC3 = (p.kf_SO2./Keq_SO_ISC3) * 1 / BP_Q * SO * s2_SQempty;
k23_SO_Q = (p.kf_SO / Keq_SO_SQ) * q_m / p.KD_Q / BP_Q * SO * s2_SQempty;
k23_SO_FAD = (kf_FADH / Keq_SO_FADH) * SO * s2_FAD;
k23_SO_FADH = (kf_FADH / Keq_SO_FADH2) * SO * s2_FADH;
k2_3 =   k23_SO_Q + k23_SO_FADH + k23_SO_FAD + k23_SO_ISC3;

% 10) E3 -> E2
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H + SO
% O2 + FADH2 -> FADH + H + SO

k32_O2_ISC3 = p.kf_SO2 * o2 * s3_ISC3 / BP_Q; 
k32_O2_SQ = p.kf_SO * o2 * s3_SQ;
k32_O2_FADH = kf_FADH * o2 * s3_FADH;
k32_O2_FADH2 = kf_FADH2 * o2 * s3_FADH2;
k3_2= k32_O2_SQ + k32_O2_FADH + k32_O2_FADH2 + k32_O2_ISC3;

% 11) E2->E4
% H2O2 + FAD -> O2 + FADH2
% FAD + Succinate -> FADH2 + fumarate
% QH2 + [4Fe_3S]ox -> SQ + [4Fe_3S]red + 2H

k24_H2O2_FAD = (kf_H2O2./Keq_H2O2_FADH2) .* H2O2 .* s2_FAD;
k24_SUC_FAD = kf_SUC.* suc_m./p.KD_SUC./BP_FAD .* s2_FAD; 
k24_QH2_SQempty_ISC3ox = (kf_QH2./Keq_QH2_ISC3) .* qh2_m./p.KD_QH2./BP_Q .* s2_SQempty_ISC3ox;
k2_4 = k24_H2O2_FAD + k24_SUC_FAD + k24_QH2_SQempty_ISC3ox;

% 12) E4->E2
% O2 + FADH2 ->  H2O2 + FAD
% FADH2 + fumarate -> FAD + Succinate
% SQ-[4Fe_3S]red + 2H -> QH2 + [4Fe_3S]ox

k42_O2_FADH2 = kf_H2O2 * o2 * s4_FADH2;
k42_FUM_FADH2 = (kf_SUC / Keq_FUM_FADH2) * fum_m / p.KD_FUM / BP_FAD * s4_FADH2; 
k42_SQ_ISC3_QH2 =  kf_QH2 * s4_SQ_ISC3;
k4_2 = k42_O2_FADH2 + k42_FUM_FADH2 + k42_SQ_ISC3_QH2;

% 13) E3->E4
% SO + Q -> O2 + SQ
% SO + FAD + H -> FADH + O2
% SO + FADH + H -> FADH2 + O2

k34_SO_ISC3 = (p.kf_SO2 / Keq_SO_ISC3) * 1 / BP_Q * SO * s3_SQempty; 
k34_SO_Q = (p.kf_SO / Keq_SO_SQ) * q_m /p.KD_Q /BP_Q * SO * s3_SQempty; 
k34_SO_FAD = (kf_FADH / Keq_SO_FADH) * SO * s3_FAD;
k34_SO_FADH = kf_FADH2 *SO * s3_FADH;
k3_4 =   k34_SO_Q + k34_SO_FAD + k34_SO_FADH + k34_SO_ISC3;

% 14) E4->E3
% O2 + SQ -> SO + Q
% O2 + FADH -> FAD + H + SO
% O2 + FADH2 -> FADH + H + SO

k43_O2_ISC3 = p.kf_SO2 * o2 * s4_ISC3 / BP_Q;
k43_O2_SQ = p.kf_SO * o2 * s4_SQ;
k43_O2_FADH = kf_FADH * o2 * s4_FADH;
k43_O2_FADH2 =kf_FADH2 * o2 * s4_FADH2;
k4_3 = k43_O2_SQ + k43_O2_FADH + k43_O2_FADH2 + k43_O2_ISC3;


%% Electron occupancy steady states (explicit solution)

D = (k0_2 * k1_0 * k2_4 * k3_1 + k0_1 * k1_2 * k2_4 * k3_1 + k0_1 * k1_3 * k2_0 * k3_4 + k0_2 * k1_0 * k2_4 * k3_2 + k0_1 * k1_2 * k2_4 * k3_2 + k0_1 * k1_3 * k2_1 * k3_4 + k0_2 * k1_0 * k2_3 * k3_4 + k0_2 * k1_2 * k2_4 * k3_1 + k0_1 * k1_2 * k2_3 * k3_4 + k0_1 * k1_3 * k2_4 * k3_2 + k0_2 * k1_0 * k2_4 * k3_4 + k0_2 * k1_2 * k2_4 * k3_2 + k0_2 * k1_3 * k2_1 * k3_4 + k0_1 * k1_2 * k2_4 * k3_4 + k0_1 * k1_3 * k2_3 * k3_4 + k0_2 * k1_2 * k2_3 * k3_4 + k0_2 * k1_3 * k2_4 * k3_2 + k0_1 * k1_3 * k2_4 * k3_4 + k0_2 * k1_2 * k2_4 * k3_4 + k0_2 * k1_3 * k2_3 * k3_4 + k0_2 * k1_3 * k2_4 * k3_4 + k0_1 * k1_3 * k2_0 * k4_2 + k0_1 * k1_3 * k2_0 * k4_3 + k0_1 * k1_3 * k2_1 * k4_2 + k0_2 * k1_0 * k2_3 * k4_2 + k0_1 * k1_2 * k2_3 * k4_2 + k0_1 * k1_3 * k2_1 * k4_3 + k0_2 * k1_0 * k2_3 * k4_3 + k0_2 * k1_3 * k2_1 * k4_2 + k0_1 * k1_2 * k2_3 * k4_3 + k0_1 * k1_3 * k2_3 * k4_2 + k0_2 * k1_0 * k2_4 * k4_3 + k0_2 * k1_2 * k2_3 * k4_2 + k0_2 * k1_3 * k2_1 * k4_3 + k0_1 * k1_2 * k2_4 * k4_3 + k0_1 * k1_3 * k2_3 * k4_3 + k0_2 * k1_2 * k2_3 * k4_3 + k0_2 * k1_3 * k2_3 * k4_2 + k0_1 * k1_3 * k2_4 * k4_3 + k0_2 * k1_2 * k2_4 * k4_3 + k0_2 * k1_3 * k2_3 * k4_3 + k0_2 * k1_3 * k2_4 * k4_3 + k0_2 * k1_0 * k3_1 * k4_2 + k0_1 * k1_2 * k3_1 * k4_2 + k0_2 * k1_0 * k3_1 * k4_3 + k0_2 * k1_0 * k3_2 * k4_2 + k0_1 * k1_2 * k3_1 * k4_3 + k0_1 * k1_2 * k3_2 * k4_2 + k0_2 * k1_0 * k3_2 * k4_3 + k0_2 * k1_2 * k3_1 * k4_2 + k0_1 * k1_2 * k3_2 * k4_3 + k0_1 * k1_3 * k3_2 * k4_2 + k0_2 * k1_0 * k3_4 * k4_2 + k0_2 * k1_2 * k3_1 * k4_3 + k0_2 * k1_2 * k3_2 * k4_2 + k0_1 * k1_2 * k3_4 * k4_2 + k0_1 * k1_3 * k3_2 * k4_3 + k0_2 * k1_2 * k3_2 * k4_3 + k0_2 * k1_3 * k3_2 * k4_2 + k0_1 * k1_3 * k3_4 * k4_2 + k0_2 * k1_2 * k3_4 * k4_2 + k0_2 * k1_3 * k3_2 * k4_3 + k0_2 * k1_3 * k3_4 * k4_2 + k0_1 * k2_0 * k3_1 * k4_2 + k0_1 * k2_0 * k3_1 * k4_3 + k0_1 * k2_0 * k3_2 * k4_2 + k0_1 * k2_1 * k3_1 * k4_2 + k0_1 * k2_0 * k3_2 * k4_3 + k0_1 * k2_1 * k3_1 * k4_3 + k0_1 * k2_1 * k3_2 * k4_2 + k0_2 * k2_1 * k3_1 * k4_2 + k0_1 * k2_0 * k3_4 * k4_2 + k0_1 * k2_1 * k3_2 * k4_3 + k0_1 * k2_3 * k3_1 * k4_2 + k0_2 * k2_1 * k3_1 * k4_3 + k0_2 * k2_1 * k3_2 * k4_2 + k0_1 * k2_1 * k3_4 * k4_2 + k0_1 * k2_3 * k3_1 * k4_3 + k0_2 * k2_1 * k3_2 * k4_3 + k0_2 * k2_3 * k3_1 * k4_2 + k0_1 * k2_4 * k3_1 * k4_3 + k0_2 * k2_1 * k3_4 * k4_2 + k0_2 * k2_3 * k3_1 * k4_3 + k0_2 * k2_4 * k3_1 * k4_3 + k1_0 * k2_0 * k3_1 * k4_2 + k1_0 * k2_0 * k3_1 * k4_3 + k1_0 * k2_0 * k3_2 * k4_2 + k1_0 * k2_1 * k3_1 * k4_2 + k1_0 * k2_0 * k3_2 * k4_3 + k1_0 * k2_1 * k3_1 * k4_3 + k1_0 * k2_1 * k3_2 * k4_2 + k1_2 * k2_0 * k3_1 * k4_2 + k1_0 * k2_0 * k3_4 * k4_2 + k1_0 * k2_1 * k3_2 * k4_3 + k1_0 * k2_3 * k3_1 * k4_2 + k1_2 * k2_0 * k3_1 * k4_3 + k1_2 * k2_0 * k3_2 * k4_2 + k1_0 * k2_1 * k3_4 * k4_2 + k1_0 * k2_3 * k3_1 * k4_3 + k1_2 * k2_0 * k3_2 * k4_3 + k1_3 * k2_0 * k3_2 * k4_2 + k1_0 * k2_4 * k3_1 * k4_3 + k1_2 * k2_0 * k3_4 * k4_2 + k1_3 * k2_0 * k3_2 * k4_3 + k1_3 * k2_0 * k3_4 * k4_2);
E0 = (k1_0 * k2_0 * k3_1 * k4_2 + k1_0 * k2_0 * k3_1 * k4_3 + k1_0 * k2_0 * k3_2 * k4_2 + k1_0 * k2_1 * k3_1 * k4_2 + k1_0 * k2_0 * k3_2 * k4_3 + k1_0 * k2_1 * k3_1 * k4_3 + k1_0 * k2_1 * k3_2 * k4_2 + k1_2 * k2_0 * k3_1 * k4_2 + k1_0 * k2_0 * k3_4 * k4_2 + k1_0 * k2_1 * k3_2 * k4_3 + k1_0 * k2_3 * k3_1 * k4_2 + k1_2 * k2_0 * k3_1 * k4_3 + k1_2 * k2_0 * k3_2 * k4_2 + k1_0 * k2_1 * k3_4 * k4_2 + k1_0 * k2_3 * k3_1 * k4_3 + k1_2 * k2_0 * k3_2 * k4_3 + k1_3 * k2_0 * k3_2 * k4_2 + k1_0 * k2_4 * k3_1 * k4_3 + k1_2 * k2_0 * k3_4 * k4_2 + k1_3 * k2_0 * k3_2 * k4_3 + k1_3 * k2_0 * k3_4 * k4_2)./D;
E1 = (k0_1 * k2_0 * k3_1 * k4_2 + k0_1 * k2_0 * k3_1 * k4_3 + k0_1 * k2_0 * k3_2 * k4_2 + k0_1 * k2_1 * k3_1 * k4_2 + k0_1 * k2_0 * k3_2 * k4_3 + k0_1 * k2_1 * k3_1 * k4_3 + k0_1 * k2_1 * k3_2 * k4_2 + k0_2 * k2_1 * k3_1 * k4_2 + k0_1 * k2_0 * k3_4 * k4_2 + k0_1 * k2_1 * k3_2 * k4_3 + k0_1 * k2_3 * k3_1 * k4_2 + k0_2 * k2_1 * k3_1 * k4_3 + k0_2 * k2_1 * k3_2 * k4_2 + k0_1 * k2_1 * k3_4 * k4_2 + k0_1 * k2_3 * k3_1 * k4_3 + k0_2 * k2_1 * k3_2 * k4_3 + k0_2 * k2_3 * k3_1 * k4_2 + k0_1 * k2_4 * k3_1 * k4_3 + k0_2 * k2_1 * k3_4 * k4_2 + k0_2 * k2_3 * k3_1 * k4_3 + k0_2 * k2_4 * k3_1 * k4_3)./D;
E2 = (k0_2 * k1_0 * k3_1 * k4_2 + k0_1 * k1_2 * k3_1 * k4_2 + k0_2 * k1_0 * k3_1 * k4_3 + k0_2 * k1_0 * k3_2 * k4_2 + k0_1 * k1_2 * k3_1 * k4_3 + k0_1 * k1_2 * k3_2 * k4_2 + k0_2 * k1_0 * k3_2 * k4_3 + k0_2 * k1_2 * k3_1 * k4_2 + k0_1 * k1_2 * k3_2 * k4_3 + k0_1 * k1_3 * k3_2 * k4_2 + k0_2 * k1_0 * k3_4 * k4_2 + k0_2 * k1_2 * k3_1 * k4_3 + k0_2 * k1_2 * k3_2 * k4_2 + k0_1 * k1_2 * k3_4 * k4_2 + k0_1 * k1_3 * k3_2 * k4_3 + k0_2 * k1_2 * k3_2 * k4_3 + k0_2 * k1_3 * k3_2 * k4_2 + k0_1 * k1_3 * k3_4 * k4_2 + k0_2 * k1_2 * k3_4 * k4_2 + k0_2 * k1_3 * k3_2 * k4_3 + k0_2 * k1_3 * k3_4 * k4_2)./D;
E3 = (k0_1 * k1_3 * k2_0 * k4_2 + k0_1 * k1_3 * k2_0 * k4_3 + k0_1 * k1_3 * k2_1 * k4_2 + k0_2 * k1_0 * k2_3 * k4_2 + k0_1 * k1_2 * k2_3 * k4_2 + k0_1 * k1_3 * k2_1 * k4_3 + k0_2 * k1_0 * k2_3 * k4_3 + k0_2 * k1_3 * k2_1 * k4_2 + k0_1 * k1_2 * k2_3 * k4_3 + k0_1 * k1_3 * k2_3 * k4_2 + k0_2 * k1_0 * k2_4 * k4_3 + k0_2 * k1_2 * k2_3 * k4_2 + k0_2 * k1_3 * k2_1 * k4_3 + k0_1 * k1_2 * k2_4 * k4_3 + k0_1 * k1_3 * k2_3 * k4_3 + k0_2 * k1_2 * k2_3 * k4_3 + k0_2 * k1_3 * k2_3 * k4_2 + k0_1 * k1_3 * k2_4 * k4_3 + k0_2 * k1_2 * k2_4 * k4_3 + k0_2 * k1_3 * k2_3 * k4_3 + k0_2 * k1_3 * k2_4 * k4_3)./D;
E4 = 1 -(E0+E1+E2+E3);



%% SDH Fluxes

% Succinate Reduction
J1 =(k02_SUC_FAD * E0 + k13_SUC_FAD * E1 + k24_SUC_FAD * E2 ...
   - k20_FUM_FADH2 * E2 - k31_FUM_FADH2 * E3 - k42_FUM_FADH2 * E4);

% Ubiquinone Reduction
J2 =(k20_SQ_ISC3_QH2 * E2 + k31_SQ_ISC3_QH2 * E3 ...
               + k42_SQ_ISC3_QH2 * E4 - k02_QH2_SQempty_ISC3ox ...
               * E0 - k13_QH2_SQempty_ISC3ox * E1 ...
               - k24_QH2_SQempty_ISC3ox * E2);
end