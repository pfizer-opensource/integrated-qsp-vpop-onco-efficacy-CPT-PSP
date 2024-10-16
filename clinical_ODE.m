function dydt = clinical_ODE(t,y,data_dictionary_dosing)

% This functions implements the model ODE, returning time derivatives for
% the passed state vector and time. The simulation time is used to
% interpolate the PK and compute the instantaneous drug concentration.

% Input
% t -- time at which to compute model time derivative for ODE
% y -- state vector at which to compute model time derivative for ODE
% data_dictionary_dosing -- dosing dictionary for pk table interpolation

% Output
% dydt -- model time derivative

pars = data_dictionary_dosing.parameters;

% 5 differential variables in pathway model (pRAS, pALK, pERK, pAKT, pS6)
% 6 differential variables in shell & core model (Nprolif,Nnecrotic,Nnec1,Nnec2,Nnec3,Nnec4)
% 3 differential variables in ALKi clinical PK model (Depot_Drug, Cc_Drug, Cp_Drug)
dydt = zeros(5+7+3,1);

%Assign parameters using parameter struct
pGFR          = pars( 1) ; % dimensionless, ALK-independent signaling
Imax_Drug     = pars( 2) ; % dimensionless, maximal pALK inhibition by ALKi
IC50_Drug     = pars( 3) ; % nM, ALKi concentration for half-maximal inhibition of pALK
h_Drug        = pars( 4) ; % dimensionless, Hill coefficient
alpha         = pars( 5) ; % dimensionless, sensitivity of proliferation to changes in pS6
KC50_prolif   = pars( 6) ; % dimensionless, value of pS6 for half maximal inhibtion of proliferation
beta          = pars(7) ; % dimensionless, sensitivity of apoptosis to changes in pAKT
KC50_apop     = pars(8) ; % dimensionless, value of pAKT for half maximal inhibtion of apoptosis
kg0           = pars(9) ; % 1/day, proliferation rate constant
kel           = pars(10) ; % 1/day, elimination of cells in nectrotic core into clearance transit compartment
tau           = pars(11) ; % days, necrotic core clearance transit compartment time constant
kout          = pars(12) ; % 1/day, de-phosphorlation rate constant of proteins
scale_kkill   = pars(13) ; % dimensionless, empirical adjustment factor to translate in-vitro calibrated cell killing to effect on tumor volume
fup_Drug     = pars(14) ; % fraction unbounds of drug in plasma

% Additional parameters
Pi   = 4*atan(1) ;
Vpmc = 1e-8*1e6  ; % volume per million cells - assuming 1e5 cells/uL

%%
%Model variables
% Pathway model
% Phosphorylated proteins
pRAS  = y(1);
pALK   = y(2);
pERK   = y(3);
pAKT   = y(4);
pS6    = y(5);

pRAS0 = 1;
pALK0  = 1;
pERK0  = 1;
pAKT0  = 1;
pS60   = 1;

%  Shell and core model
%  Cellular balances described below:
%  Nprolif   - cells in proliferating outer shell
%  Nnecrotic - cells that have entered the necrotic core
%  Nnec_k    - cells in transit compartment as they are cleared from the necrotic core
Nprolif    = y(5+1);
Nnecrotic  = y(5+2);
Nnec1      = y(5+3);
Nnec2      = y(5+4);
Nnec3      = y(5+5);
Nnec4      = y(5+6);
Nkilled    = y(5+7);

if data_dictionary_dosing.schedule_index == 0
    
    C_cell_Drug   =  0;

elseif data_dictionary_dosing.schedule_index == 1
    
    %%% USE LOOK-UP TABLE FOR Drug CONCENTRATION if available
    
    Drug_PK_table = data_dictionary_dosing.cc_Drug;
    [~, ind]        = unique(Drug_PK_table(:,1));
    Drug_PK_interp = interp1(Drug_PK_table(ind,1),Drug_PK_table(ind,2),t);
    
    C_cell_Drug   = fup_Drug*Drug_PK_interp;
end

%Compute inhibitor effect on within cell phosphorylated signalling protein
%Additional sensitivity effect of pALK on pRAS production
I_Drug     = Imax_Drug*max(1e-20,C_cell_Drug)^h_Drug/(IC50_Drug^h_Drug+max(1e-20,C_cell_Drug)^h_Drug) ;

%Phosphorylation (i.e., "production") of signaling proteins
dpRAS   = ((pALK + pGFR)/(1 + pGFR))*kout*pRAS0 - kout*pRAS ;
dpALK    = (1-I_Drug)*kout*pALK0                               - kout*pALK  ;
dpERK    = pRAS*kout*pERK0                                     - kout*pERK  ;
dpAKT    = pRAS*kout*pAKT0                                     - kout*pAKT  ;
w_pS6    = 1                                                                 ;
dpS6     = (w_pS6*pERK + (1-w_pS6)*pAKT)*kout*pS60              - kout*pS6   ;

% Proliferation rate constant depends on pS6
x_prolif = pS6;
fprolif = (KC50_prolif^alpha+1)*max(1e-20,x_prolif)^alpha/(KC50_prolif^alpha+max(1e-20,x_prolif)^alpha);
% Apoptosis rate depends on inhibition of pAKT
x_apop   = 1-pAKT;
%fapop    = (max(1e-9,x_apop)^beta)/(KC50_apop^beta+max(1e-9,x_apop)^beta);
fapop    = max(1e-20,x_apop)^beta/(KC50_apop^beta+max(1e-20,x_apop)^beta) ;
%%%%
% Define killing rate for connection to tumor model - empirical scaling for in-vitro to in-vivo translation
kkill = (1 - (fprolif-fapop))*scale_kkill;

Vshell     = Vpmc*(Nprolif/1e6); % Vpmc is volume per million cells, so scale Nprolif by 1e6
Vcore      = Vpmc*((Nnecrotic+Nnec1+Nnec2+Nnec3+Nnec4)/1e6); % ... also scale total number of cells in necrotic core by 1e6
Vtumor     = Vshell + Vcore;

Rtumor     = (3*Vtumor/(4*pi))^(1/3);
Rcore      = (3*Vcore /(4*pi))^(1/3);
shell      = Rtumor - Rcore ;

% Constant percent necrotic core, see slides for derivation
knecrotic = (Nnec4/tau + (Vcore/Vtumor)*((1 - kkill)*kg0*Nprolif - Nnec4/tau))/Nprolif;
%%%%%%%

dNprolif   = ((1 - kkill)*kg0  - knecrotic)*Nprolif;
dNnecrotic = knecrotic*Nprolif - kel*Nnecrotic;
dNnec1     = kel*Nnecrotic     - Nnec1/tau;
dNnec2     = (Nnec1            - Nnec2)/tau;
dNnec3     = (Nnec2            - Nnec3)/tau;
dNnec4     = (Nnec3            - Nnec4)/tau;
dNkilled   = kkill*Nprolif     - kel*Nkilled ;

% Return RHS calulations to the solver
dydt = [dpRAS; dpALK; dpERK; dpAKT; dpS6; dNprolif; dNnecrotic; dNnec1; dNnec2; dNnec3; dNnec4; dNkilled];
end