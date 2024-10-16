function power_analysis_weighted_samp_list=power_analysis(param_values,param_bounds,prob_inclusion,param_names,top_params,store_sims,layout,plot_only,power_analysis_weighted_samp_list)
    
    % This function performs a power analysis across varying trial sizes
    % for the top parameters identified in responder_rank_params. At each
    % sample size 1000 trials are simulated using sampling with replacment
    % from the plausible population. A rank-sum test is performed between
    % (non-)responders for each top parameter value. The fraction of
    % signification results across this set of trials is used to compute
    % the expected statistical power for the given sample size. This is
    % used to generate a power analysis plot (Fig. 5c)

    %Input
    % plausible_scaled_param -- 
    % prob_inclusion -- probability of inclusiuon computed from get_vpop,
    %                   a normalized ratio of probility of seeing given
    %                   plausible patient w.r.t. clinical data divided by
    %                   KNN density estimate of plausible patients in
    %                   output space. Sums to 1.
    % plausible_responder_stat
    % param_names -- string array of parameter names
    % top_params
    % layout -- layout handle for sub-plots in following functions

    %set observation and virtual population colors
    obspop_color = [0 0.4470 0.7410];
    vpop_color =[0.8500 0.3250 0.0980];
    %create color maps for observed and virtual populations
    color_map_obs = linspace(.7,.1,10)'.*repmat(obspop_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    color_map_vpop = linspace(.7,.1,10)'.*repmat(vpop_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);

    %get number of plausible patients (default to 10k for paper)
    num_plausibles = length(store_sims);
    %scale parameter values as % betwn bounds
    plausible_scaled_param = 100*((param_values - param_bounds(:,1))./(param_bounds(:,2) - param_bounds(:,1)))';
    %extract the BPC and responder status for the plausible patients
    plausible_BPC = cell2mat(arrayfun(@(x) store_sims(x).results{2}.best_dSLD(1), 1:num_plausibles, 'UniformOutput', false))';
    plausible_responder_stat = plausible_BPC<-30;

    %set list of trial samples sizes to evaluate (x-axis)
    sample_sizes = round(logspace(log10(20),log10(2000),10));
    %set number of trials to run at each sample size for power calc.
    num_trial_sims = 1000;
    
    %power_analysis_weighted_samp_list={};
    %create an array to store the fraction of significant trials (power) at
    %each sample size
    power = zeros(length(top_params),length(sample_sizes));
    %loop over the samples sizes
    for i=1:length(sample_sizes)
        %loop over the trials to be done at each sample size
        for j=1:num_trial_sims
            %sample a weighted set of indicies from the plausibles,
            %according to probability-of-inclusion, to yield desired sample
            %size, note: samples with replacment
            if(~plot_only)
                weighted_sample = subsamp(sample_sizes(i),prob_inclusion);
                power_analysis_weighted_samp_list{i,j}=weighted_sample;
            else
                weighted_sample=power_analysis_weighted_samp_list{i,j};
            end
            %subset the plausible parameter matrix, using above 
            weighted_scaled_param = repelem(plausible_scaled_param,weighted_sample,1);
            %convert weighted/subsetted parameter matrix above into a table
            %this allowed subsetting with parameter names (not essential
            %but method I used here due to prior analysis versions)
            weighted_scaled_param_table = array2table(weighted_scaled_param,'VariableNames',  param_names');
            %subset/replicate the responder lables of the original
            %plausible population for the new parameter sample above
            weighted_responder_stat = repelem(plausible_responder_stat,weighted_sample);
            %create a table of parameters for the responding patients
            values_responder = weighted_scaled_param_table(weighted_responder_stat==1,:);
            %create a table of parameters for the non-responding patients
            values_nonresponder = weighted_scaled_param_table(weighted_responder_stat==0,:);
            %loop over the top parameters 
            for p=1:length(top_params)
                %check if both responder types are present in the sample 
                % (at very small sample sizes (~20), sometimes no 
                % non-responders are selected)
                if sum(~weighted_responder_stat)>0 && sum(weighted_responder_stat)>0
                    %perform rank-sum test on given top parameter
                    p_val = ranksum(values_nonresponder{:,char(top_params(p))},...
                                    values_responder{:,char(top_params(p))});
                else
                    %if no non-responders were sampled, no way to achieve
                    %significance, automatically non-significant
                    p_val = 1e3;
                end
                %check for significance at 5% level
                if p_val<0.05
                    %if significant increment count of num. significant
                    %trials
                    power(p,i) = power(p,i)+1;
                end
            end
        end
    end
    %normalize count of num. significant trial by total number run at each
    %sample size
    power = 100*power/num_trial_sims;
    
    %preparr subplot for power analysis curves
    nexttile%(layout,10,[2,2])
    semilogx(sample_sizes',power(1,:)','-.','Color',vpop_color)
    hold on
    semilogx(sample_sizes',power(2,:)',':','Color',vpop_color)
    semilogx(sample_sizes',power(3,:)','--','Color',vpop_color)
    semilogx(sample_sizes',power(4,:)','-','Color',vpop_color)
    plot([min(sample_sizes) max(sample_sizes)],[80 80],'k:','LineWidth',1)
    xticks([20 50 100 250 500 1000 2000]);
    xlim([min(sample_sizes) max(sample_sizes)])
    l=legend(top_params,'Location','southeast');
    l.FontSize = 10;
    xlabel('Virtual Trial Sample Size')
    ylabel('Power (%)')
end