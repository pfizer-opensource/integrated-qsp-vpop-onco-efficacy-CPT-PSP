function [results, constraints] = simulate_clinical_dosing(TSTART,TSTOP,INTERVAL,dropout_time,censor_label,parameters,cc_Drug,schedule_index,...
    screentime,time_data,initial_condition_vector,output_times,state_bnds,score_fn,only_constraints)

% This function simulates the clinical ODE model with the loaded pk table
% and process the simulation output using
% process_model_parameters_and_states. Simulation processing also performs
% contraints checks on each simulation.

%Input
% TSTART -- trial start time
% TSTOP -- trial stop time
% INTERVAL -- time interval
% dropout_time -- drop out time for simulation
% censor_label -- censoring label for given simulation
% parameters -- parameter values
% cc_Drug -- pk table data structure
% schedule_index -- list of pk indices to simulate
% screentime -- trial screening tiume
% time_data -- observation times from trial data
% initial_condition_vector -- initial condition vector for simulation
% output_times -- times to compute outputs at
% state_bnds -- bounds for the ODE state vector
% score_fn -- function handle for scoring the simulation
% only_constraints -- boolean indicating that only contraints should be
%                     checked

%Output
% results -- processed results from the simulation
% constraints -- boolean indicating contstrain satisfaction

num_schedules = length(schedule_index);
time = cell(num_schedules,1);
states = cell(num_schedules,1);
results = [];
constraints = zeros(size(state_bnds,1),1);
% score_final = 1e256;
valid = ones(num_schedules,1);

if ~any(isnan(parameters))
    for i = 1:num_schedules 
        data_dictionary_dosing = struct();
        data_dictionary_dosing.parameters = parameters;
        if schedule_index(i) == 1
            data_dictionary_dosing.cc_Drug = cc_Drug{i};
        end
        
        data_dictionary_dosing.schedule_index = schedule_index(i);
        
        if isnan(INTERVAL(i))
            TSIM = [TSTART,TSTOP(i)];
        else
            TSIM = TSTART:INTERVAL(i):TSTOP(i);
            if max(TSIM) ~= TSTOP(i)
                TSIM = [TSIM, TSTOP(i)];
            end
        end
        
        results_init = process_model_parameters_and_states(TSTART,initial_condition_vector',parameters,schedule_index(i),screentime(i),time_data,dropout_time,censor_label,TSTART);
        % Call the ODE solver -
        %         options= odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
        %         options= odeset('RelTol',1e-6,'AbsTol',1e-6,'JPattern',calculate_jacobian_pattern_pk_table(),'Stats','on');
        if schedule_index(i) == 0
            options= odeset('MStateDependence','none', ...
            'MassSingular','no','JPattern',jacobian_pattern(),'Stats','off',...
                'Events',@(t,x)checkdoublingtime(t,x,results_init.Vtumor,parameters,schedule_index(i),screentime(i),time_data,t));
        else
            options= odeset('MStateDependence','none', ...
            'MassSingular','no','JPattern',jacobian_pattern(),'Stats','off');
        end
        
        [t,x] = ode15s(@(t,x) clinical_ODE(t,x,data_dictionary_dosing),TSIM,initial_condition_vector,options);
        
        % Check for negatives -
        idx_n = x<0;
        x(idx_n) = 0.0;
        
        if schedule_index(i) == 0 && x(end,6) < x(1,6)
            valid(i) = 0;
        else
            valid(i) = 1;
        end
        time{i} = t;
        states{i} = x;
        if isnan(output_times{i})
            output_times{i} = t;
        end
    end
    if all(valid)
        results = arrayfun(@(n) process_model_parameters_and_states(time{n},states{n},parameters,schedule_index(n),screentime(n),time_data,dropout_time,censor_label,output_times{n}), ...
            1:num_schedules, 'UniformOutput',false);
        if only_constraints
            [~,~,~,constraints] = feval(score_fn,results,[],only_constraints);
        end
    end
end
end

function [value, isterminal, direction] = checkdoublingtime(t,x,init_vol,pars,schedule_index,screentime,time_data,at)
if schedule_index == 0
    results_t = process_model_parameters_and_states(t,x',pars,schedule_index,screentime,time_data,[],[],at);
    value = (results_t.Vtumor - 2*init_vol) >= 0;
else
    value = 0;
end
isterminal = 1;
direction = 0;
end