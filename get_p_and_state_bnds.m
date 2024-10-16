function [p_bnds,init_bnds,state_bnds] = ...
    get_p_and_state_bnds(get_state_bnds_only,varargin)

% This function formulates the lower and upper bounds on the fitting parameters and
% states of the model. This function is entirely model dependent and must
% be modified for each individual model.

% Inputs
% get_state_bnds_only -- boolean indicating if only state bounds are needed
% varargin -- variable length input argument

% Outputs
% p_bnds -- parameterc bounds
% init_bnds -- initial condition bounds
% state_bnds -- model state vector bounds

% Prior Notes:
% NOTE(1): THIS FUNCTION IS HIGHLY MODEL-SPECIFIC AND MUST BE MODIFIED FOR OTHER
% MODELS.
% NOTE(2): This function is a modified version of the original function
% "Input_Ranges".

%% Set bounds on model states - these are bounds on model states that are constrained.

% shell thickness
state_lower(1) = 0.07; % mm
% state_lower(1) = 0.45; % mm
% state_upper(1) = 0.33; % mm
state_upper(1) = 1.3; % mm
state_bestfit(1) = mean([state_lower(1), state_upper(1)]);

% prct_necrotic_core
state_lower(2) = 14.6;
% state_upper(2) = 63.17;
state_upper(2) = 75;
% state_upper(2) = 45;
state_bestfit(2) = mean([state_lower(2), state_upper(2)]);

% tumour doubling time
state_lower(3) = 26;  
state_upper(3) = 2000; 
state_bestfit(3) = mean([state_lower(3), state_upper(3)]);

if get_state_bnds_only
    state_bnds = [state_lower(:) state_upper(:)];
    p_bnds = [];
    init_bnds = [];
else
    p_values    = varargin{1};
    p_names     = varargin{2};
    param_table = varargin{3};
    init_values = varargin{4};
    init_names  = varargin{5};
    init_table  = varargin{6};
    search_mag  = varargin{7};
    
    %% Initalize
    if isempty(p_values) || isempty(p_names)
        p_bestfit = param_table.value;
        p_names = param_table.name;
    else
        p_bestfit = p_values;
    end
    num_p = numel(p_values);
    
    if isempty(init_values) || isempty(init_names)
        init_bestfit = init_table.value;
        init_names = init_table.name;
    else
        init_bestfit = init_values;
    end
    num_init = numel(init_values);
    
    % Initialize arrays:
    p_lower = -1*ones(num_p,1);
    p_upper = -1*ones(num_p,1);
    init_lower = -1*ones(num_init,1);
    init_upper = -1*ones(num_init,1);
    
    %% Set bounds on model states - these are bounds on model states that are varied in the cohort.
    % use bounds from Excel file else
    % apply a uniform scalar to constrain from the best fit value
    for i = 1:length(init_names)
        ind = find(strcmp(init_table.name,init_names(i)));
        if ind
            if ~isempty(init_table.lower_bnd(ind))
                init_lower(i) = init_table.lower_bnd(ind);
            else
                init_lower(i) = init_bestfit(i)/search_mag;
            end
            if ~isempty(init_table.upper_bnd(ind))
                init_upper(i) = init_table.upper_bnd(ind);
            else
                init_upper(i) = init_bestfit(i)*search_mag;
            end
        end
    end
    %% Set bounds on model parameters:
    % use bounds from Excel file else
    % apply a uniform scalar to constrain.
    for i = 1:length(p_names)
        ind = find(strcmp(param_table.name,p_names(i)));
        if ind
            if ~isempty(param_table.lower_bnd(ind))
                p_lower(i) = param_table.lower_bnd(ind);
            else
                p_lower(i) = p_bestfit(i)/search_mag;
            end
            if ~isempty(param_table.upper_bnd(ind))
                p_upper(i) = param_table.upper_bnd(ind);
            else
                p_upper(i) = p_bestfit(i)*search_mag;
            end
        end
    end
    %% Package outputs:
    p_bnds = [p_lower(:) p_upper(:)];
    init_bnds = [init_lower(:), init_upper(:)];
    state_bnds = [state_lower(:) state_upper(:)];
end
end