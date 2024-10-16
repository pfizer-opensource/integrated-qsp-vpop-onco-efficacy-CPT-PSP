function [select,hist_score,hist_mu,hist_std,p_incl,sf] = get_vpop(pd,mu,sigma,varargin)

% This function implements the Allen/Rieger virtual population selection
% algorithm for generating a virtual population from an existing plausible
% population.

% Input
% pd -- distribution object used for scoring virtuals
% mu -- means of gaussian mixture distribution object for scoring
% sigma -- covariance of gaussian mixture distribution object for scoring
% varargin -- variable length input arguments

% Output
% select -- boolean matrix indicating which of plausibles are selected
% hist_score -- scores from repeated selection used to find best pop
% hist_mu -- mean score from repeated selection used to find best pop
% hist_std -- standard dev. of score from repeated selection for best pop
% p_incl -- probability of inclusion, ratio of data likelihood to plausible density
% sf -- scaling factor

%% Parse the Varargin
p = inputParser;
addParameter(p,'run_sims_flag',1,@(x)ismember(x,[1;2]));
addParameter(p,'scaling_factor',[])
addParameter(p,'do_select',1,@(x) ismember(x,[0,1]));
addParameter(p, 'data_file','data.csv',@(x)ischar(x));
addParameter(p, 'data_to_fit',{'SLD', 'BESTPCHG','EEVALUMP'},@(x)iscell(x));
addParameter(p, 'ln_transform',[1,0,0], @(x) isvector(x));
addParameter(p, 'gm_comp',2, @(x) isnumeric(x));
addParameter(p, 'vars_to_fit',{'SLD_mm','best_dSLD','time_to_pfs'},@(x) iscell(x));
addParameter(p,'results_file','results.mat',@(x)ischar(x));
addParameter(p,'output_root','vpop_results_',@(x)ischar(x));
addParameter(p,'screentime',0,@(x) isnumeric(x));
parse(p,varargin{:});

%% Alias inputs
run_sims_flag  = p.Results.run_sims_flag;  % <--- SETS NEW RUN OR RERUN (1 = NEW VPOP, 2 = PLOT OLD VPOP)
scaling_factor = p.Results.scaling_factor; % acceptance-rejection sampling algorithm specific parameter. Default is empty.
do_select      = p.Results.do_select;      %
data_file      = p.Results.data_file;      % fullpath of the data file (.csv) used for fitting
data_to_fit    = p.Results.data_to_fit;    % columns of the data file to use for fitting/scoring. Ensure this matches the column names in the data.
ln_transform   = p.Results.ln_transform;   % should the variable be log-transformed: boolean vector of size data_to_fit.
gm_comp        = p.Results.gm_comp;        % number of components in GM distribution
vars_to_fit    = p.Results.vars_to_fit;    % model outcomes to use for fitting - names may be different, but correspond to the names in data_to_fit
results_file   = p.Results.results_file;   % fullpath of results file with plausible patients to be used for fitting
output_root    = p.Results.output_root;    % the prefix to be used for output file. num_vps, schedule_flag, date and time are added later.
screentime     = p.Results.screentime;     % in days - gap between screening and first dose - does not affect the fit. Only used for visualisation

% Set global variables for figure plotting
set(0, 'DefaultAxesFontSize',18)
set(0, 'DefaultErrorBarLineWidth', 1);
set(0, 'DefaultLineLineWidth',2)

do_plot = 0;

data = readtable(data_file,'TreatAsEmpty',{'NA'});
%[pd,mu_data,~,data] = get_distribution_from_data(data,data_to_fit,ln_transform,gm_comp,do_plot);
len_vars = length(data_to_fit);
r = data{:,data_to_fit};
ln_ind = find(ln_transform);
r(:,ln_ind) = log(r(:,ln_ind));
r = r(~any(isnan(r'),1)',:); % Clear any rows with NaN, not absolutely necessary, but avoids a warning message
data=r;

data = (data - min(data,[],1))./(max(data,[],1) - min(data,[],1));

switch run_sims_flag
    case 1
        %%% setup:
        load(results_file)
        
        num_pps = size(store_sims,2); % number of plausible patients
        num_obs = 3; %%Nate edit%% %numel(mu_data);     % number of distributions to match (1d histograms)
        
        if length(vars_to_fit) == 1
            pps_obs(:,1) = cell2mat(arrayfun(@(x) log(store_sims(x).results{1}.(vars_to_fit{1})(1)), 1:num_pps, 'UniformOutput', false))';
        else
            pps_obs = [cell2mat(arrayfun(@(x) log(store_sims(x).results{1}.(vars_to_fit{1})(1)), 1:num_pps, 'UniformOutput', false))', ...
                cell2mat(arrayfun(@(y) cell2mat(arrayfun(@(x) store_sims(x).results{2}.(vars_to_fit{y})(1), ...
                1:num_pps, 'UniformOutput', false))',2:length(vars_to_fit),'UniformOutput',false))];
        end
        
        output_root = sprintf('%s',p.Results.output_root,num2str(num_pps),datestr(now, '_mmm_dd_yyyy'), datestr(now, '_hh_MM'));
        
        %%% generate test data:
%         %%% Data input check:
%         if (num_obs > 1)
%             mu1     = repmat(mu_data,num_pps,1); % for use in scaling VPs
%             sig_pps    = repmat(diag(sigma_data)',num_pps,1);
%             pps_obs_scal = (pps_obs-mu1)./sig_pps.^0.5;                % Normalize VPs
%             data_pdf   = mvnpdf(pps_obs_scal.*sig_pps.^0.5+mu1,mu_data,sigma_data); % Data density
%             
%         else
%             pps_obs_scal = (pps_obs - mu1)./sigma_data;
%             data_pdf = normpdf(pps_obs_scal.*sigma_data+mu1,mu_data,sigma_data);
%         end
        %%% rescale PP observations
        
       data_pdf = pdf(pd,pps_obs);
       pps_obs_scal = (pps_obs-min(pps_obs,[],1))./(max(pps_obs,[],1) - min(pps_obs,[],1));

%         data_pdf = mvksdensity(data,pps_obs_scal,'bandwidth',.05);

        %%% calculate the prob. of inclusion for each plausible patient:
        % Probability density (PDF) at the VP's location divided by the density of
        % plausible patients in the same region. The plausible patient density is
        % created by taking the volume of the sphere to the num_pts nearest
        % neighbors (default = 5).        

        vp_density  = mvksdensity(pps_obs_scal,pps_obs_scal,'bandwidth',.05)';
        p_incl      = (data_pdf)./vp_density';  % Uniform probability of inclusion
        
        %%% optimize the scaling factor for selection of VPs:
        runs    = 10; % Number of times to try the fitting
        fopt    = @(p)get_vp_select(p,p_incl,pps_obs_scal,runs,do_select,data);
        sf_max  = 1/max(p_incl); % Maximum scaling factor
        sf_lower = 0; % lower bound for the scaling factor on probability of inclusion
        sf_upper = log10(1000*sf_max); % upper bound for the scaling factor on the probability of inclusion
        
        options_sa = saoptimset('TolFun',1e-15,'ReannealInterval',50000,...
            'InitialTemperature',0.5,'MaxIter',1000,'Display','off',...
            'TemperatureFcn',@temperatureboltz,'AnnealingFcn',@annealingboltz,...
            'AcceptanceFcn',@acceptancesa);
        %   options=saoptimset('TolFun',1e-15,'ReannealInterval',50000,'InitialTemperature',1,'MaxIter',1000,'TemperatureFcn',@temperatureboltz,'AnnealingFcn', @annealingboltz,'AcceptanceFcn',@acceptancesa);
        
        if isempty(scaling_factor)
            k   = simulannealbnd(fopt,log10(sf_max),sf_lower,sf_upper,options_sa); % optimize the log10(scaling factor) value
            sf  = 10.^k(1); % scaling factor transformed back from log10
        else
            k = log10(scaling_factor);
            sf = scaling_factor;
        end
        
        %%% with the Optimal Scaling Factor, [Re]select VPs:
        % Iterate several times (default = 100) and select the best fit and mean
        % fits.
        num_tests = 100;
        hist_score_temp = zeros(num_tests,1);
        selected_temp   = zeros(num_pps,num_tests);
        for i=1:num_tests
            [hist_score_temp(i),selected_temp(:,i)] = get_vp_select(k,p_incl,pps_obs_scal,1,do_select,data);
        end
        [hist_score,I]  = min(hist_score_temp);
        select          = logical(selected_temp(:,I));
        hist_mu         = nanmean(hist_score_temp);
        hist_std        = nanstd(hist_score_temp);
        
        %         save(fullfile('SimOutput', output_root), 'store_sims', ...
        %             'p_pp', 'p_names', 'p_bnds', 'p_vary_ind', ...
        %             'init_pp', 'init_names', 'init_bnds', 'init_vary_ind', 'v_pp', 'v_names','v_bnds',...
        %             'dropout_pp','state_bnds', 'schedule_flag', 'screentime', 'starttime',...
        %             'TSTART','TSTOP','INTERVAL','mu','sigma', ...
        %             'select', 'hist_score', 'hist_mu', 'hist_std', 'p_incl', 'sf', '-v7.3');
        
        %verify_cohort(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,'data',data,'select',select);
        %paper_plots(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,'data',data);
        %paper_plots(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,'data',data,'select',select);
    case 2
        %load(results_file)
        %verify_cohort(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,'data',data,'select',select);
        %paper_plots(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,'data',data,'select',select);
end
end % function get_prevalence

% ************************************************************************
% Helper functions:

function vol = get_nsphere_vol(n,radius)
%% function get_nsphere_vol
% Calculates the volume of a n-dimensional sphere
% formula can be found on Wikipedia.

vol = pi^(n/2).*(1/gamma(n/2+1)).*radius.^n;

end

% ************************************************************************
function [hist_score,select,p_out] = ...
    get_vp_select(sf_log10, p_incld, pps, runs,do_select,data)
%% function get_vp_select(sf_log10, p_incld, pps, runs)
%   Selects virtual patients from the cohort of plausible patients
%   previously generated. And scores them for normality using a series of
%   K-S tests.
%
%   Inputs:
%       sf_log10 - 1x1 scalar, log10(scaling factor) for prob. of inclusion
%       p_incld - nx1 vector, probability of inclusion for each plausible
%           patient (pre-scaling-factor).
%       pps - nxm matrix, plausible patients for inclusion
%       runs - 1x1 scalar, number of times to attempt the fitting
%
%   Where n is the number of plausible patients amd m is the number of
%   observations being matched to data.
%
%   Outputs:
%       hist_score  - 1x1 scalar, sum of all univariate K-S tests across
%           all runs.
%       select      - nx1 boolean,
%       p_out       - 1x1 scalar, p-value for K-S tests
%

%%% setup
sf = 10.^sf_log10(1); % scale factor input is log10
data_dim = numel(pps(1,:));

% Pre-allocate output vectors:
p_val        = zeros(data_dim,1);
ksstat      = zeros(data_dim,1);
hist_score  = zeros(runs,1);
num_pps     = size(pps,1);

if ~do_select
    p_incld = ones(num_pps,1)/sf;
end

%%% perform Selection and Univariate K-S Tests:
for j = 1:runs
    r = rand(num_pps,1);
    select = r < p_incld*sf;
    num_vp = sum(select);
    if num_vp > 1
        for i = 1:data_dim
            %Test if dist. is close to the distribution of the data
            [~,p_val(i),ksstat(i)] = kstest2(pps(select,i),data(:,i));
            %             [h,p_val(i),ksstat(i)] = kstest(pps(select,i));
        end
        hist_score(j) = sum(ksstat);
        p_out = p_val;
    else
        hist_score(j) = 1*data_dim;
    end
end
hist_score = sum(hist_score)/runs;

end % function get_vp_select
