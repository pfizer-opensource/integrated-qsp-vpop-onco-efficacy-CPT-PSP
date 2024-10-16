function [initial_conditions, parameters,other_states] = adjust_initial_conditions_and_parameters(initial_conditions,parameters,...
    other_states,bounds,randomize)

    % This function adjusts the initial conditions after new parameters
    % have been proposed. Adjustments enforce several constrains on the
    % shell-and-core model and randomize initial conditions subject to
    % certain assumptions (i.e. necrotic compartments all have same IC).
    % This function will also adjust k_el after IC's have been set.

    % Input
    % initial_conditions -- proposed initial conditions, to be adjusted
    % parameters -- parameters values
    % other_states -- values of cell ratio and thickness, if to be adjusted
    % bounds -- bounds for state variables
    % randomize -- boolean indicating if cell ratio should be randomized

    % Output
    % initial_conditions -- adjusted initial conditions
    % parameters -- adjusted parametr values (note kel is adjusted w.r.tICs)
    % other_states -- adjusted cell ratio and shell thickness

    Ntransit = 4;

    % starting with kel < 0, adjust initial conditions and parameters till kel > 0
    iter = 0;
    while true
        kg0 = parameters(9);
        tau = parameters(11);
        
        initial_conditions(6:12) = initial_conditions(6:12)*1e9;
        
        % adjust intial conditions based on parameters being varied
        switch length(other_states)
            case 1 % varying cell ratio and Nprolif
                total_necrotic = other_states(1)*initial_conditions(6);
            case 2 % varying cell ratio and shell (other states set in this order) and setting Nprolif accordingly.
                initial_conditions(6) = other_states(2)^3/((3*1e-8/(4*pi()))*((1+other_states(1))^(1/3) - other_states(1)^(1/3))^3);
                total_necrotic = other_states(1).*initial_conditions(6);
            otherwise
                error('other_states has more than 2 variables');
        end
        
        % only run this is initial conditions need to be adjusted - allows for
        % reproducibility of results with fixed parameters and ICs 
        % (e.g.: PDC calibration experiments)
        if (total_necrotic ~= sum(initial_conditions(7:11)))
            %% adjust initial conditions based on parameter values to ensure tumour grows without therapy
            max_nnec4 = kg0*initial_conditions(6)*tau;
            min_r = (total_necrotic - Ntransit*max_nnec4)/total_necrotic;
            
            r = min_r + (1-min_r).*rand(1,1); %% assuming that Nnecrotic is at least minr percent of the total necrotic cells
            initial_conditions(7) = r*total_necrotic;
            initial_conditions(8:11) = (total_necrotic - initial_conditions(7))/Ntransit;
        end
        % adjust parameters based on initial conditions
        phi_necrotic = (sum(initial_conditions(7:11)))./(initial_conditions(6) + sum(initial_conditions(7:11)));
        phi_necrotic_ss = phi_necrotic/0.985;
        
        kel = kg0*(phi_necrotic_ss-1)/(Ntransit*tau*kg0*(1-phi_necrotic_ss) - phi_necrotic_ss);
        parameters(10) = kel;
        
        initial_conditions(6:12) = initial_conditions(6:12)/1e9;
        
        if kel > 0
            break;
        else
            %%% (on line 39) since prct_necrotic_core < 100, phi_necrotic_ss < 1, 
            %%% thus, numerator < 0 => kel > 0 when denominator < 0
            min_phi_necrotic = Ntransit*tau*kg0/(1 + Ntransit*tau*kg0)*0.985;
            %%% since bounds(1,2) corresponds to cell_ratio, convert min_phi_necrotic to cell_ratio
            min_cell_ratio = (min_phi_necrotic/(1 - min_phi_necrotic));
            if min_cell_ratio > bounds(1,2) % if no feasible solution is possible, return NAs
                parameters = nan(size(parameters));
                break;
            else
                if (randomize)
                    other_states(1) = min_cell_ratio + (bounds(1,2) - min_cell_ratio).*rand(1,1);
                else
                    other_states(1) = bounds(1,2);
                end
                iter = iter+1;
            end
        end
    end