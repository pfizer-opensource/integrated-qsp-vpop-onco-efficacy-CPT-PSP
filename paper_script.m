%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs all computationa and plotting for the main
% and supplemental figures in the corresponding paper (Braniff et al., 2024)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date last updated: 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set to true to reproduce exact paper figures from saved simulations
plot_only = true;
% set to true if generating new plausibles, otherwise they are loaded
gen_plausibles = false;
%use the original saved mixture fit to clinical data (vs re-fitting)
use_original_mixture = true;
% set to true if overwriting saved plots is desired
save_plots = false;

%% User Defined Settings
% File paths
path_root     = './Data';
data_file     = [path_root '/synthetic_clinical_data.csv'];     % fullpath of the data file (.csv) used for fitting
param_file    = [path_root '/params.xlsx']; % fullpath of the parameter file to use
init_file     = [path_root '/initial_conditions.xlsx'];     % fullpath of the initial conditions file to use
pk_file       = [path_root '/pk_table.xlsx'];
results_path  = ['./Output'];
output_root   = 'new_synthetic_plausible_population';   % the prefix to be used for output file.

% User defined inputs
schedule_flag   = {'none','Drug'};                        % schedule to simulate - none corresponds to no treatment. Provide all schedules names as a cell array.
num_mh_pps         = 10000;                               % Number of plausible patients to attempt to create with MH
data_to_fit     = {'SLD', 'BESTPCHG', 'EEVALUMP'};        % columns of the data file to use for fitting/scoring. Ensure this matches the column names in the data. 
ln_transform    = [1,0,0];                                % should the variable be log-transformed: boolean vector of size data_to_fit.
gm_comp         = 2;                                      % number of components in GM distribution
vars_to_score   = {'SLD_mm','best_dSLD','time_to_pfs'};   % model outcomes to use in the scoring function - names may be different, but correspond to the names in data_to_fit
TSTOP           = [165,500];                              % simulation stop time - vector same size as schedule_flag
screentime      = [0,0];                                  % in days - gap between screening and first dose
starttime       = [0,0];                                  % in days - start of therapy
change          = {};                                     % these arguments will not be loaded from the results file, so values from 'varargin' will be used. (doesn't apply to run_sims_flag = 1)
verbose         = 1;                                      % is true, will print messages on screen for every step.
do_plot         = 0;                                      % to plot countours and histograms after fitting PDF to data

% Set global variables for figure plotting
set(0, 'DefaultAxesFontSize',18)
set(0, 'DefaultErrorBarLineWidth', 1);
set(0, 'DefaultLineLineWidth',2)

% Add paper plotting functions to the path from subfolder
addpath('./PaperPlottingFunctions')

%% Initialization and Data/File Loading

if(gen_plausibles)
    % Parse patient data; measurement and dropout times 
    data = readtable(data_file,'TreatAsEmpty',{'NA'});
    if(use_original_mixture)
        load('Scoring_Gaussian_mixture_fit.mat');
    else
        [pd, mu, sigma] = get_distribution_from_data(data,data_to_fit,ln_transform,gm_comp,do_plot);
        logme('Joint PDF fitted to data', verbose)
    end
    
    % Get time points from data - VPs will be evaluated at these times for best
    % % change, PFS etc. 
    time_data = cell2mat(arrayfun(@(x) str2double(cell2mat(regexp(data.Properties.VariableNames{x}, '\d*', 'match'))),...
        1:size(data,2),'UniformOutput', false));
    time_data(isnan(time_data)) = [];
    
    % get vector of dropout times from the data.
    [dropout_times_data,sort_idx] = sortrows(data.EEVALUMP);
    censor_label_data = cell2mat(arrayfun(@(x)strcmp('NO PFS EVENT', data.EVENTNP_NEW{x}), 1:size(data,1),'UniformOutput', false))';
    censor_label_data = censor_label_data(sort_idx);
    dropout_na = isnan(dropout_times_data);
    dropout_times_data(dropout_na) = [];
    censor_label_data(dropout_na) = [];
    
    %load param and ICs
    default_param_table = readtable(param_file);
    init_table = readtable(init_file);
    
    % load dosing schedule for ALKi drug
    sim_days      = init_table.value(strcmp(init_table.name,'Total_Days'));
    TSTART        = 0;
    if isempty(TSTOP) || any(sim_days >  TSTOP)
        TSTOP         = sim_days;
    end
    init_table = init_table(1:end-1,:);
    
    % Load the PK table
    [schedule_index, cc_Drug, TSTOP, INTERVAL] = load_pk_table(pk_file,schedule_flag,TSTOP,screentime,starttime, @(m)logme(m,verbose));
    schedule_flag_op = sprintf('%s_', schedule_flag{:});
    output_file = sprintf('%s',output_root,num2str(num_mh_pps),'_',schedule_flag_op,datestr(now, '_mmm_dd_yyyy'), datestr(now, '_hh_MM'));

end

%% Select Plausible Patients

if(gen_plausibles)
    % Use M-H algorithm to generate plausible patients based on baseline SLD and mean spider distributions        
    tic
    [p_pp,p_names,p_bnds,p_vary_ind,init_pp,init_names,init_bnds,init_vary_ind,v_pp,v_names,v_bnds,dropout_pp,censor_pp,pp_yield,state_bnds] = mh_generate_pps(num_mh_pps,...
        default_param_table,init_table,pd,mu,sigma,screentime,TSTART,TSTOP,cc_Drug,schedule_index,...
        vars_to_score,time_data,dropout_times_data,@(m)logme(m,verbose),censor_label_data);
    toc
    save(fullfile(results_path, output_file), 'p_pp', 'p_names', 'p_bnds', 'p_vary_ind', ...
        'init_pp', 'init_names', 'init_bnds', 'init_vary_ind', 'v_pp', 'v_names','v_bnds',...
        'dropout_pp','censor_pp','state_bnds','schedule_flag','screentime','starttime',...
        'TSTART','TSTOP','INTERVAL','pd','mu','sigma','-v7.3');
end

%% Simulate Plausible Patients
if(gen_plausibles)
    num_patients = size(p_pp,2);
    if (num_patients == num_mh_pps)
        store_sims = run_vps(num_mh_pps,p_pp,init_pp,TSTART,TSTOP,INTERVAL,cc_Drug,schedule_index,screentime,time_data,state_bnds,dropout_pp,censor_pp);
    else
        warning('only %d patients created out of %d total required', num_patients, num_mh_pps);
    end
    save(fullfile(results_path, output_file), 'store_sims', ...
        'p_pp', 'p_names', 'p_bnds', 'p_vary_ind', ...
        'init_pp', 'init_names', 'init_bnds', 'init_vary_ind', 'v_pp', 'v_names','v_bnds',...
        'dropout_pp','censor_pp','state_bnds', 'schedule_flag', 'screentime', 'starttime',...
        'TSTART','TSTOP','INTERVAL','pd','mu','sigma','-v7.3');
    
    results_file = fullfile(results_path, output_file);
end

%% Loading an Existing Plausible Population 

if(~gen_plausibles)
    %add path the plotting helper functions sub-directory
    addpath('./PaperPlottingFunctions') 
    %set path the plausible population results file
    results_file = './Output/plausible_population_Jan_04_2024_19_38.mat';
    %load plausible pop. results file
    load(results_file)
    %load the clinical data
    data = readtable(data_file,'TreatAsEmpty',{'NA'});
end

% if(save_plots)
%     save('./Scoring_Gaussian_mixture_fit','pd','mu','sigma');
% end

%% Diagnostic Plots for Plausible Patient MH Generation (Supplementary Fig.S2-4)

%parameter names
param_names = ["GFR","IC50_{Drug}","h_{Drug}","\alpha",...
               "KC50_{prolif}","\beta","KC50_{apop}",...
               "k_{g0}","\tau","\rho_{core}","\delta_{shell}"]';
%get parameter values for plausible population
param_values = v_pp;
%get parameter bounds for plausible population
param_bounds = v_bnds;
%generate diagnostic plots for plausible MH run
[mixing_fig,autocorr_fig,pairwise_fig]=mh_diagnostic_plots(param_values,...
                                                           param_bounds,...
                                                           param_names);
%save figures in ./Figures folder (Note: this will overwrite existing)
if(save_plots)
    saveas(pairwise_fig,'./Figures/FigS2_MH_Pairwise_Parameters.png')
    saveas(autocorr_fig,'./Figures/FigS3_MH_Parameter_Autcorr.png')
    saveas(mixing_fig,'./Figures/FigS4_MH_Parameter_Mixing.png')
end

%% Sample a Virutal Population & Compute Probability of Inclusion

if(~plot_only)
    %select a virtual pop from plausibles using Allen/Rieger alg.
    % (non-deterministic sample size), also return probability of inclusion
    [vpop_inclusion_bool,~,~,~,...
     inclusion_ratio,~] = get_vpop(pd,mu,sigma,'data_file',data_file,...
                                                'results_file',results_file);
    %normalize the inclusion ratio so it is a proper probability distribution
    %over the plausible population
    prob_inclusion = inclusion_ratio/sum(inclusion_ratio);
end

%% Downsampling and Plotting a Virtual Population for Comparison to Data (Figure 4)

%desired sample size for single Vpop trial
target_sample_size = 155;%sum(vpop_inclusion_bool);
if(~plot_only)
    %downsample the the virtual population to the desired sample size
    subset_selected = downsample(target_sample_size,vpop_inclusion_bool);
    km_bootstrap_selected_list={};
    power_analysis_weighted_samp_list={};
else
    load('./Vpop_Plotting.mat');
end
%generate paper plots to compare virtual population to observed data
compare_fig = compare_plot(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,...
                          TSTART,TSTOP,INTERVAL,pd,mu,sigma,plot_only,'data',data,'select',subset_selected);
%save figures in ./Figure folder (Note: this will overwrite existing)
if(save_plots)
    saveas(compare_fig,'./Figures/Fig4_Vpop_Example_Comparison.png')
end

%% Simulate Kaplan-Meier PFS Estimate via Vpop Bootstrap Resampling (Figure 5) 

%set number of bootstrapped trials to run for KM interval
num_replicate_trials = 1000;
%set sample size for simulated trials
sample_size = 155;
%boostrap the KM estimate for the Vpop and compare to data, plot figure
[km_uncertainty_fig,...
 km_bootstrap_selected_list] = sim_km_bootstrap(sample_size,num_replicate_trials,...
                                                prob_inclusion,data,store_sims,plot_only,km_bootstrap_selected_list);
%save figure in ./Figures folder (Note: this will overwrite existing)
if(save_plots)
    saveas(km_uncertainty_fig,'./Figures/Fig5_KM_Uncertainty.png')
end

%% Vpop Parameter Analysis (Figure 6)

%boolean indicating an older ver. of Matlab is in use w/out boxchart
old_matlab = false;
%Figure 6A -- Responder Difference and Distributions of Vpop Parameters
%rank-sum test parameters on responder status, rank and box plots top params
[weight_scaled_param,...
 weighted_responder_stat,top_params,...
 layout,parameter_analysis_fig] = responder_rank_params(param_values,...
                                                        param_bounds,...
                                                        param_names,...
                                                        prob_inclusion,...
                                                        store_sims,old_matlab);

%Figure 6B -- Power Analysis for Responder Difference in Parameter Values
%perform power analysis on top ranked responder status via bootstrapping
power_analysis_weighted_samp_list = power_analysis(param_values,param_bounds,prob_inclusion,param_names,top_params,store_sims,layout,plot_only,power_analysis_weighted_samp_list);
%save figures in ./Figures folder (Note: this will overwrite existing)
if(save_plots)
    saveas(parameter_analysis_fig,'./Figures/Fig6_Parameter_Analysis.png')
end

%if(~plot_only)
%    save('./Vpop_Plotting','vpop_inclusion_bool','prob_inclusion','subset_selected','km_bootstrap_selected_list','power_analysis_weighted_samp_list');
%end

%% Helper Functions

function [store_sims] = run_vps(num_vps,p_pp,init_pp,TSTART,TSTOP,INTERVAL,cc_Drug,schedule_index,screentime,time_data,state_bnds,dropout_pp,censor_pp)
    parfor i = 1:num_vps
        %Print patient index being evaluated by the parfor pool
        str= sprintf('evaluating patient %d',i);
        disp(str)
        parameters = p_pp(:,i);
        initial_conditions = init_pp(:,i);
        initial_conditions(6:12) = initial_conditions(6:12)./1e9;
        
        %%%%%%%%%%%%%
        output_times = cell(length(schedule_index),1);
        output_times(:) = {nan(1,1)};
        
        tstart = tic;
        results = simulate_clinical_dosing(TSTART,TSTOP,INTERVAL,dropout_pp(i),censor_pp(i),parameters,cc_Drug,...
            schedule_index,screentime,time_data,initial_conditions,output_times,state_bnds,0,0);
        results{1}.clock = toc(tstart);
        
        clc
        store_sims(i).results = results;
    end 
end

function logme(messagetoLog,verbose)
    if verbose
        fprintf('%s\n', messagetoLog)
    end
end

