function results = process_model_parameters_and_states(t,x,pars,schedule_index,screentime,time_data,dropout_time,censor_label,at)

% This function process the output of the ODE simulation of the data and
% computes clinically relevant outputs like BPC, pSLD etc. If a PFS/dropout
% time is passed, endpoints will be computed with the specified value
% given.

% Input
% t -- time vector from ODE simulation
% x -- state vector from ODE simulation
% pars -- parameters used for ODE simulation
% schedule_index -- pk schedule index used for simulation
% screentime -- screen time for trial
% time_data -- obsevation times from data
% dropout_time -- specified dropout time to use for computation (optional)
% censor_label -- censoring lable to be assigned (optional w/ dropout time)
% at -- indices for output times

% Output
% results -- results data structure

temp_t = t;

x(:,6:12) = x(:,6:12)*1e9; %% convert cell numbers to correct scale.

if (~isempty(dropout_time) && schedule_index ~= 0)
    t(temp_t > dropout_time) = nan;
    x(temp_t > dropout_time, :) = nan;
end

alpha         = pars( 5) ; % dimensionless, sensitivity of proliferation to changes in pS6
KC50_prolif   = pars( 6) ; % dimensionless, value of pS6 for half maximal inhibtion of proliferation
beta          = pars(7) ; % dimensionless, sensitivity of apoptosis to changes in pAKT
KC50_apop     = pars(8) ; % dimensionless, value of pAKT for half maximal inhibtion of apoptosis
kg0           = pars(9) ; % 1/day, proliferation rate constant
kel           = pars(10) ; % 1/day, elimination of cells in nectrotic core into clearance transit compartment
tau           = pars(11) ; % days, necrotic core clearance transit compartment time constant
scale_kkill   = pars(13) ; % dimensionless, empirical adjustment factor to translate in-vitro calibrated cell killing to effect on tumor volume


results.pRAS        = x(:,1);
results.pALK        = x(:,2);
results.pERK        = x(:,3);
results.pAKT        = x(:,4);
results.pS6         = x(:,5);

results.Nprolif     = x(:,5+1);
results.Nnecrotic   = x(:,5+2);
results.Nnec1       = x(:,5+3);
results.Nnec2       = x(:,5+4);
results.Nnec3       = x(:,5+5);
results.Nnec4       = x(:,5+6);
results.Nkilled     = x(:,5+7);
Vpmc                = 1e-8*1e6  ; % volume per million cells - assuming 1e5 cells/uL
Ntransit            = 4;

%%% some intermediate calculations
x_prolif = results.pS6;
fprolif  = (KC50_prolif.^alpha+1).*max(1e-20,x_prolif).^alpha./(KC50_prolif.^alpha+max(1e-20,x_prolif).^alpha);
x_apop   = 1-results.pAKT;
fapop    = max(1e-20,x_apop).^beta./(KC50_apop.^beta+max(1e-20,x_apop).^beta) ;
results.kkill    = (1 - (fprolif-fapop))*scale_kkill;

% observable model outputs
results.Vshell       = Vpmc*(results.Nprolif/1e6); % Vpmc is volume per million cells, so scale Nprolif by 1e6
results.Vcore        = Vpmc*((results.Nnecrotic+results.Nnec1+results.Nnec2+results.Nnec3+results.Nnec4)/1e6); % ... also scale total number of cells in necrotic core by 1e6
results.Vtumor       = results.Vshell + results.Vcore;

results.Rtumor       = (3*results.Vtumor/(4*pi)).^(1/3);
results.Rcore        = (3*results.Vcore/(4*pi)).^(1/3);
results.shell        = results.Rtumor - results.Rcore ;
results.SLD_mm       = results.Rtumor*2*10;
results.SLD_mm_nadir = results.SLD_mm;
results.diffSLD      = results.SLD_mm - results.SLD_mm_nadir;

results.dSLD                 = 100*(-results.SLD_mm(1)+results.SLD_mm)./results.SLD_mm(1);
results.dSLD_nadir           = 100*(-results.SLD_mm_nadir+results.SLD_mm)./results.SLD_mm_nadir;
results.prct_necrotic_core   = 100*results.Vcore./results.Vtumor;
results.prct_baseline_Rtumor = 100*(-results.Rtumor(1)+results.Rtumor)./results.Rtumor(1);
results.size_baseline_change = ((2*results.Vtumor).^(1/3)-(2*results.Vtumor(1))^(1/3))/((2*results.Vtumor(1))^(1/3))*100;

results.knecrotic             = (results.Nnec4./tau + (results.Vcore./results.Vtumor).*((1 - results.kkill).*kg0.*results.Nprolif - results.Nnec4./tau))./results.Nprolif;
results.prct_necrotic_core_ss(1:size(x,1),1) = 100*(kg0./kel + tau.*kg0.*Ntransit)./(1 + kg0./kel + tau.*kg0.*Ntransit);

therapy_start_ind = find(t == screentime);
time_sim = [0; t(therapy_start_ind+1:end)/7 - screentime/7]; % converting days to weeks and adjusting for gap between screening and C1D1, if any.
inds_to_check = find(ismember(time_sim, time_data));

% calculate best dSLD and PFS time - only at scanning intervals
if length(inds_to_check)==1
    if size(x,1) > 1
        results.best_dSLD(1:size(x,1),1)     = min(results.dSLD([2, size(x,1)]));
    else
        results.best_dSLD(1:size(x,1),1)     = results.dSLD(1);
    end
    results.time_to_best(1:size(x,1),1)  = t(find(results.dSLD == results.best_dSLD,1,'first'));
    if ~isempty(dropout_time)
        results.time_to_pfs(1:size(x,1),1) = dropout_time;
        results.censor_label(1:size(x,1),1) = censor_label;
    else
        results.time_to_pfs(1:size(x,1),1)   = t(end);
        results.censor_label(1:size(x,1),1) = 0;
    end
else
    % recalculating SLD_nadir and dSLD_nadir only at scanning intervals
    results.SLD_mm_nadir(inds_to_check)  = cummin(results.SLD_mm(inds_to_check));
    results.diffSLD(inds_to_check)       = results.SLD_mm(inds_to_check) - results.SLD_mm_nadir(inds_to_check);
    results.dSLD_nadir(inds_to_check)    = 100*(-results.SLD_mm_nadir(inds_to_check)+...
        results.SLD_mm(inds_to_check))./results.SLD_mm_nadir(inds_to_check);
    
    results.best_dSLD(1:size(x,1),1)     = min(results.dSLD(inds_to_check(2:end)));
    results.time_to_best(1:size(x,1),1)  = t(find(results.dSLD == results.best_dSLD,1,'first'));
    if ~isempty(dropout_time)
        results.time_to_pfs(1:size(x,1),1) = dropout_time;
        results.censor_label(1:size(x,1),1) = censor_label;
    else
        target_pfs_time = t(inds_to_check(find(...
            results.dSLD_nadir(inds_to_check) > 20 & results.diffSLD(inds_to_check) > 5, 1, 'first')));
        if ~isempty(target_pfs_time)
            results.time_to_pfs(1:size(x,1),1) = target_pfs_time;
            results.censor_label(1:size(x,1),1) = 0;
        else
            results.time_to_pfs(1:size(x,1),1) = t(end);
            results.censor_label(1:size(x,1),1) = 0;
        end
    end
end

% rescale the cell numbers before storing
% this is needed to ensure correct calculations when 
% process_model_parameters_and_states is called in multiple places
x(:,6:12) = x(:,6:12)/1e9; 

results.x = x;
results.t = t;

doubled = find(results.Vtumor >= 2*results.Vtumor(1),1,'first');
if ~isempty(doubled)
    results.tum_doubling_time(1:size(x,1),1) = t(doubled);
else
    results.tum_doubling_time(1:size(x,1),1) = nan;
end

% get indices corresponding to output times
at = cell2mat(arrayfun(@(x) find(temp_t >= at(x), 1, 'first'), 1:length(at), 'UniformOutput',false));
results = structfun(@(x) x(at,:) , results, 'UniformOutput', false);
results.pars = pars;
end