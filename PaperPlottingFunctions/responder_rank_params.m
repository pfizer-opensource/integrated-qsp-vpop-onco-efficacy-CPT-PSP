function [weight_scaled_param,weighted_responder_stat,...
          top_params,layout,responder_diff_fig] = responder_rank_params(param_values,param_bounds,...
                                                             param_names,prob_inclusion,store_sims,old_matlab)

    % Creates a very large virtual population by resampling plasuibles
    % according to probabiltiy of inclusion with replacment. Performs
    % rank-sum test on this large population over each parameter comparing
    % responders (>-30% BOR) vs non-responders. Plots top ranked parameters
    % in paired box plot to show quantiles of (non-)responders. (Fig. 5a)

    %Input
    % param_values -- matrix of plausible parameter values, col=patients &
    %                  rows=parameters, should be in MH iteration order
    % param_bounds -- table of bounds, rows are parameters, col1=min &
    %                 col2=max
    % param_names -- string array of parameter names
    % prob_inclusion -- probability of inclusiuon computed from get_vpop,
    %                   a normalized ratio of probility of seeing given
    %                   plausible patient w.r.t. clinical data divided by
    %                   KNN density estimate of plausible patients in
    %                   output space. Sums to 1.
    % store_sims -- data structure for storing simulation of the plausible
    %               population

    %Output
    % weight_scaled_param -- parameter matrix of the very large virtual
    %                        population which has been sampled with
    %                        replacment. Parameters are scaled as % between
    %                        min and max value.
    % weighted_responder_stat -- responders status for the very large
    %                            virtual population
    % layout -- layout handle for sub-plots in following functions
    % responder_diff_fig -- figure handle for parameter responder boxplot

    %set observation and virtual population colors
    respndr_color = [0.4940, 0.1840, 0.5560];
    noresp_color = [0.9290, 0.6940, 0.1250];
    
    %create color maps for observed and virtual populations
    color_map_resp = linspace(.7,.1,10)'.*repmat(respndr_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    color_map_norsp = linspace(.7,.1,10)'.*repmat(noresp_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);

    %get number of plausible patients (default to 10k for paper)
    num_plausibles = length(store_sims);
    %Number of Vpop Parameters (includes shell thickness and core-shell ratio)
    num_par=length(param_names);
    %scale parameter values as % betwn bounds
    plausible_scaled_param = 100*((param_values - param_bounds(:,1))./(param_bounds(:,2) - param_bounds(:,1)))';
    %select the number of virtual patients to plot
    %note: all plausible patients are used in the analysis, replicated
    %proportional to their probability of inclusion, number of virtual patients
    %determines the number of total replicates (should be >> plausible pop)
    num_virtual = 1e4;
    %convert prob inclusion to counts via num virtuals (round for whole num)
    weighted_plausible_counts = round(prob_inclusion*num_virtual);
    %re-sum to get exact number of total replicates over plausible pop
    num_tot_virtual = sum(weighted_plausible_counts);

    %extract the BPC and responder status for the plausible patients
    plausible_BPC = cell2mat(arrayfun(@(x) store_sims(x).results{2}.best_dSLD(1), 1:num_plausibles, 'UniformOutput', false))';
    plausible_responder_stat = plausible_BPC<-30;
    
    %replciate responder and param values according to weighted plausible counts
    weighted_responder_stat = repelem(plausible_responder_stat,weighted_plausible_counts);
    weight_scaled_param = repelem(plausible_scaled_param,weighted_plausible_counts,1);
    
    %create a long-form table of param names, param values, and responder stat
    %each replicated plausible patient takes up 12 rows, one per parameter
    %note: tables are needed for boxchart grouping, hence converting to that form
    Responder = repelem(weighted_responder_stat,num_par);
    Param = repmat(param_names,num_tot_virtual,1);
    Value = reshape(weight_scaled_param',prod(size(weight_scaled_param)),1);
    weighted_param_table = table(Responder,Param,Value);
    weighted_param_table.Param = categorical(weighted_param_table.Param,param_names);
    
    %create arrays for the rank-sum pvalues and median differences for each
    %parameter
    rank_sum_pvals=zeros(num_par,1);
    median_diffs=zeros(num_par,1);
    %loop over parameter names
    for i=1:length(param_names)
        %subset parameter values on responder status
        val_responder = weighted_param_table(weighted_param_table.Param==param_names(i)&weighted_param_table.Responder==1,:).Value;
        val_nonresponder = weighted_param_table(weighted_param_table.Param==param_names(i)&weighted_param_table.Responder==0,:).Value;
    
        %perform rank-sum test and record p-value
        rank_sum_pvals(i) = ranksum(val_nonresponder,val_responder);
        %compute median difference in parameter value and record
        median_diffs(i) = median(val_nonresponder)-median(val_responder);
    end
    
    %create a table of rank-sum p-values and median differences, one row per
    %parameter
    parameter_diffs_table = table(param_names,rank_sum_pvals,median_diffs);
    parameter_diffs_table.param_names = categorical(parameter_diffs_table.param_names,...
                                                    parameter_diffs_table.param_names);
    parameter_diffs_table = sortrows(parameter_diffs_table,'rank_sum_pvals','descend');
    
    %number of parameters (from lowest p-value up) to plot
    top_n = 4;
    %extract the top parameter names
    top_params = parameter_diffs_table.param_names(end-(top_n-1):end);
    %ensure weighted_param_table has same categorical names as parameter_diff_table
    weighted_param_table.Param = categorical(weighted_param_table.Param,parameter_diffs_table.param_names);
    %create a boolean array picking out only parameters of interest from
    %weighted_param_table array
    top_bool = ismember(weighted_param_table.Param,top_params);
    %subset weighted_param_table_top_subset selecting only top parameter values
    weighted_param_table_top_subset = weighted_param_table(top_bool,:);
    %removed unused categories for non-top parameters
    weighted_param_table_top_subset.Param = removecats(weighted_param_table_top_subset.Param);
    
    %uese box chart to plot parameter distributions by responder status
    responder_diff_fig = figure('Position', [1, 1, 1200, 500]);
    layout = tiledlayout(1,2);
    nexttile%(layout,1,[2,2])
    if(~old_matlab)
        %Note: Boxchart command does not work in older Matlab versions (i.e. Improve)
        boxchart(grp2idx(weighted_param_table_top_subset.Param)*2,weighted_param_table_top_subset.Value,'GroupByColor',weighted_param_table_top_subset.Responder,'Orientation','horizontal','MarkerStyle','.')
        colororder([noresp_color; respndr_color]); 
        yticks((1:top_n)*2);
        yticklabels(top_params)
        legend()
        xlabel('% Value Between Feasible Bounds')
        ylabel('Parameters')
        xlim([-10 110])
        ylim([-1.5,top_n*2+1])
        legend({'Non-responder','Responder'},'Location','south');
    else
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(1)&weighted_param_table_top_subset.Responder==0),'position',1,'orientation','horizontal','Colors',noresp_color,'BoxStyle',"filled",'Symbol','.')
        hold on
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(1)&weighted_param_table_top_subset.Responder==1),'position',1.25,'orientation','horizontal','Colors',respndr_color,'BoxStyle',"filled",'Symbol','.')
        
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(2)&weighted_param_table_top_subset.Responder==0),'position',2,'orientation','horizontal','Colors',noresp_color,'BoxStyle',"filled",'Symbol','.')
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(2)&weighted_param_table_top_subset.Responder==1),'position',2.25,'orientation','horizontal','Colors',respndr_color,'BoxStyle',"filled",'Symbol','.')
        
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(3)&weighted_param_table_top_subset.Responder==0),'position',3,'orientation','horizontal','Colors',noresp_color,'BoxStyle',"filled",'Symbol','.')
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(3)&weighted_param_table_top_subset.Responder==1),'position',3.25,'orientation','horizontal','Colors',respndr_color,'BoxStyle',"filled",'Symbol','.')
        
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(4)&weighted_param_table_top_subset.Responder==0),'position',4,'orientation','horizontal','Colors',noresp_color,'BoxStyle',"filled",'Symbol','.')
        boxplot(weighted_param_table_top_subset.Value(weighted_param_table_top_subset.Param==top_params(4)&weighted_param_table_top_subset.Responder==1),'position',4.25,'orientation','horizontal','Colors',respndr_color,'BoxStyle',"filled",'Symbol','.')
        
        ylim([-.5,4.5])
        yticks([1.125,2.25,3.125,4.125])
        yaxisproperties= get(gca, 'YAxis');
        yaxisproperties.TickLabelInterpreter = 'tex'; 
        yticklabels(top_params)
        
        xlabel('% Value Between Feasible Bounds')
        ylabel('Parameters')
        hLegend = legend(findall(gca,'Tag','Box'), {'Responder','Non-responder'},'Location','south');
    end

end