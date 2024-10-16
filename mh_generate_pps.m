function [p_pp,p_names,p_bnds,p_vary_ind,init_pp,init_names,init_bnds,init_vary_ind,v_pp,v_names,v_bnds,dropout_pp,censor_pp,pp_yield,state_bnds] = ...
    mh_generate_pps(num_pps,param_table,init_table,pd,mu,sigma,screentime,TSTART,TSTOP,...
    cc_Drug,schedule_index,vars_to_score,time_data,dropout_times_data,verbose_fun,censor_label_data)

% This function implements the core Metropolis Hasting algorithm for
% sampling plausible patients. This includes the top-level iterative loop
% for proposing, simulating and scoring candidate patients, as well as the
% dropout and censoring assignment mechanism.

% Input
% num_pps -- number of desired plausible patients
% param_table -- parameter table with values bounds and names
% init_table -- IC table with state values bounds and names
% pd -- distribution object used for scoring plausibles
% mu -- means of gaussian mixture distribution object
% sigma -- covariance of gaussian mixture distribution object
% screentime -- prefixed screening time for trial
% TSTART -- start time for trial
% TSTOP -- stop time for trial
% cc_Drug -- drug pk table
% schedule_index -- indices of pk tables used in trial simulation/selection
% vars_to_score -- variable names used in scoring plausible patients
% time_data -- observations times from actual trial, used for simulation
% dropout_times_data -- observed dropout times from actual trial
% verbose_fun -- function handle for printing verbose feedback
% censor_label_data -- censoring labels for each observed dropout time in
%                      data

% Output
% p_pp -- parameter table (all params) for plausible patients
% p_names -- parameter names of all params
% p_bnds -- parameter bounds for all params
% p_vary_ind -- indices of parameters that were varried in the plaus. pop.
% init_pp -- IC table for plausibl patients
% init_names -- IC state names
% init_bnds -- IC state bounds 
% init_vary_ind -- indices of initial conditions to be varied in plaus. pop.
% v_pp -- table of params and ICS specifically varied in the plaus. pop. 
% v_names -- names for params and ICs varied in the plaus. pop.
% v_bnds -- bounds for params and ICs varied in the plaus. pop.
% dropout_pp -- dropout times for each plausible patient
% censor_pp -- censoring times for each plausible patient
% pp_yield -- yield of plausible patients accepted
% state_bnds -- state bounds for model

%% Load parameters
parameters = param_table.value;
p_names = param_table.name;
p_vary_ind = param_table.isvariable == 1;
p_values = param_table.value(p_vary_ind);
num_p = numel(p_values);

% Load initial conditions (ICs)
initial_conditions = init_table.value;
initial_conditions(6:12) = initial_conditions(6:12)/1e9;
init_names = init_table.name;
init_vary_ind = init_table.isvariable == 1;
init_values = init_table.value(init_vary_ind);
init_values = log10(init_values);

% calculate/fetch the bounds for the parameters and states:
search_mag = 5;
[p_bnds,init_bnds,state_bnds] = ...
    get_p_and_state_bnds(0, parameters,p_names,param_table,initial_conditions,init_names,init_table,search_mag);
init_bnds = log10(init_bnds);

% get cell ratio and shell thickness for intial vector of parameters and ICs
init_results = process_model_parameters_and_states(TSTART,initial_conditions',parameters,schedule_index(1),screentime(1),time_data,[],[],TSTART);
cell_ratio = init_results.prct_necrotic_core./(100-init_results.prct_necrotic_core);
shell = init_results.shell;
% if varying initial conditions, vary only cell ratio.
% else determine initial conditions based on cell ratio and shell thickness
if sum(init_vary_ind) > 0
    v_values = [p_values; init_values; cell_ratio];
    num_vary = numel(v_values);
    v_bnds = [p_bnds(p_vary_ind,:); init_bnds(init_vary_ind,:); state_bnds(2,:)./(100-state_bnds(2,:))];
    v_names = [p_names(p_vary_ind); init_names(init_vary_ind); 'cell_ratio'];
else
    v_values = [p_values; shell; cell_ratio];
    num_vary = numel(v_values);
    v_bnds = [p_bnds(p_vary_ind,:); state_bnds(2,:)./(100-state_bnds(2,:)); state_bnds(1,:)];
    v_names = [p_names(p_vary_ind); 'cell_ratio'; 'shell'];
end

%Alias the scoring function to set some passed parameters
score_func = @(x,dropout_time_vp,only_constraints)mh_score_model(x,state_bnds,schedule_index,screentime,vars_to_score,time_data,dropout_times_data,...
    pd,mu,sigma,dropout_time_vp,only_constraints,censor_label_data); % set the scoring function for the model

%Alias the simulation function to set some passed parameters
simulate_func = @(parameters,initial_conditions,TSTOP,INTERVAL,schedule_index,output_times,only_constraints)simulate_clinical_dosing(TSTART,TSTOP,INTERVAL,...
    [],[],parameters,cc_Drug,schedule_index,screentime,time_data,initial_conditions,output_times,state_bnds,score_func,only_constraints); % set the simulation function for the model

feval(verbose_fun, 'M-H setup completed')

k = 1;
iter = 1;
max_iter = num_pps*1000;

%create matrices to store parameters,ICs, dropouts and censoring for
%plausibes patients that are selected
p_pp       = nan(length(parameters),num_pps);
init_pp    = nan(length(initial_conditions),num_pps);
v_pp       = nan(length(v_values),num_pps);
dropout_pp = nan(num_pps,1);
censor_pp = nan(num_pps,1);
best_dsld  = nan(num_pps,1);
        
% adjust initial conditions and parameters
% to ensure model constraints are satisfied and kel > 0
v = v_values;
parameters(p_vary_ind) = v(1:num_p);
if sum(init_vary_ind) > 0
    initial_conditions(init_vary_ind) = v(num_p+1:end-1);
    [initial_conditions, parameters,v(end)] = adjust_initial_conditions_and_parameters(initial_conditions,parameters,v(end),v_bnds(end,:),1);
else
    [initial_conditions, parameters,v(end-1:end)] = adjust_initial_conditions_and_parameters(initial_conditions,parameters,v(end-1:end),v_bnds(end-1:end,:),1);
end

%run initial conditions to check/populated initial results and score
r1 = simulate_func(parameters, initial_conditions,TSTOP,[state_bnds(3,1)-1,7],schedule_index,{TSTART,nan},0);
assert(~isempty(r1),'not a valid initialisation, change parameters or initial conditions');
s1 = score_func(r1,TSTOP(2),0); % assuming no dropout in the initial vector
feval(verbose_fun,'Initial parameter vector simulation done')

%Main Metropolis Hasting loop
while k <= num_pps && iter <= max_iter

    %Proposal
    % propose new parameter set:
    randv = rand(num_vary,1);
    q = max(v_bnds(:,1),min(v_bnds(:,2),v + (v_bnds(:,2)-v_bnds(:,1))/10.*(2*randv-1)));
    parameters(p_vary_ind) = q(1:num_p);

    %IC randomization and adjustment
    % adjust initial conditions and parameters
    % to ensure model constraints are satisfied and kel > 0
    if sum(init_vary_ind) > 0
        initial_conditions(init_vary_ind) = q(num_p+1:end-1);
        [initial_conditions, parameters,q(end)] = adjust_initial_conditions_and_parameters(initial_conditions,parameters,q(end),v_bnds(end,:),1);
    else
        [initial_conditions, parameters,q(end-1:end)] = adjust_initial_conditions_and_parameters(initial_conditions,parameters,q(end-1:end),v_bnds(end-1:end,:),1);
    end

    %Simulation
    %call simulated to do quick check of constraints
    [~,c] = simulate_func(parameters, initial_conditions,state_bnds(3,1),nan,0,{TSTART},1);
    % if constraints are valid, do full simulation
    if all([c(1:2); ~c(3)])
        r2  = simulate_func(parameters, initial_conditions,TSTOP,[state_bnds(3,1)-1,7],schedule_index,{TSTART,nan},0); %Simulate the model for parameters
        [s2,d2,l2]  = score_func(r2,[],0); %Score the simulation output(r2); d2 contains dropout_time_vp and param constraints, s2 is vec of score
    else
        s2 = 1e256;
        r2 = {process_model_parameters_and_states(TSTART,initial_conditions',parameters,schedule_index(1),screentime(1),time_data,[],[],TSTART)};
    end
    if any(s2 == 1e256)
        s2 = eps; % correct the error checking
    end

    
    %Acceptance
    %generate random threshold
    r = rand();
    % find indices of plausible patient copies that could possibly be accepted.
    to_accept = find(r < s2/s1 & s2 > eps); % Vector testing if the ratio of probabilities is larger than the random number r
    if ~isempty(to_accept) % If any patient is acceptable (i.e. meets threshold)
        % randomly choose one copy as the final PP from all the accepted ones  
        select = randi(length(to_accept),1,1);
        if mod(k,1) == 0
            fprintf('Found: %d PPs\n',k);
        end
        %store details of accepted plausible patient
        p_pp(:,k) = parameters(:); %Save the patient params
        init_pp(:,k) = initial_conditions(:); %ICs
        init_pp(6:12,k) = init_pp(6:12,k)*1e9; %Rescale cell numbers
        v_pp(:,k) = [parameters(p_vary_ind); q(num_p+1:end)]; %plausible population
        dropout_pp(k) = d2(to_accept(select)); %Save dropout time for the accepted Vpatient
        censor_pp(k) = l2(to_accept(select)); 
        best_dsld(k) = r2{2}.dSLD(find(r2{2}.t >= dropout_pp(k),1,'first')); 
        v = q;
        s1 = s2(to_accept(select));
        k = k + 1;
 
    end
    iter = iter + 1;
    parameters = param_table.value;
    initial_conditions = init_table.value;
end

if iter > max_iter
    error('Exceeded maximum iterations');
end

pp_yield = (k-1)/(iter-1);

init_bnds = 10.^(init_bnds); % converting back from log scale to save results correctly.

end % function mh_generate_pps
