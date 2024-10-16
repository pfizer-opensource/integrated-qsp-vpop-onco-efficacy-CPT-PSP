function [mixing_fig,...
    autocorr_fig,...
    pairwise_fig]=mh_diagnostic_plots(param_values,param_bounds,param_names)

    %Generates MH diagnostic plots consisting of i) scaled parameter values
    %(% of min-max range) vs MH iteration ii) lagged autocorrelation of
    %parameter values iii) pairwise scatter for all parameters with
    %histograms of each parameter marginal. (Supplementary FigS2-4)

    %Inputs
    % param_values -- matrix of plausible parameter values, col=patients &
    %                  rows=parameters, should be in MH iteration order
    % param_bounds -- table of bounds, rows are parameters, col1=min &
    %                 col2=max
    % param_names -- string array of parameter names

    %Outputs
    % mixing_fig -- figure handle for the parameter vs iter fig (i)
    % autocorr_fig -- figure handle for the lagged autocorr fig (ii)
    % pairwise_fig -- figure handle for the pairwise scatter/hist (iii)

    %extract MH parameter values for plausibles and scale between bounds
    plausible_param_values_scaled = 100*((param_values - param_bounds(:,1))./(param_bounds(:,2) - param_bounds(:,1)))';
    %extract number of plausibles (10,000 default for paper)
    num_plausibles = size(plausible_param_values_scaled,1);
    
    %create figure for parameter mixing plots
    mixing_fig=figure('Position', [1, 1, 1000, 800]);
    t1=tiledlayout(4,3);
    %create figure for parameter autocorrelation plots
    autocorr_fig=figure('Position', [1, 1, 1000, 800]);
    t2=tiledlayout(4,3);
    %Loop over Vpop parameters
    for i=1:11
        %plot MH parameter value chain to show mixing
        nexttile(t1)
        plot(plausible_param_values_scaled(:,i))
        title(param_names(i))
        %compute and plot MH parameter value autocorrelation 
        ac=xcov(plausible_param_values_scaled(:,i),'normalized');
        nexttile(t2)
        plot(ac(num_plausibles:end))
        title(param_names(i))
    end
    %labels and titles
    title(t1,'MH Plausible Parameter Mixing')
    xlabel(t1,'MH Iteration')
    ylabel(t1,'Relative Parameter Value (%)')
    title(t2,'MH Plausible Parameter Autocorrelation')
    xlabel(t2,'Lag')
    ylabel(t2,'Autocorrelation')
    
    %create figure for pairwise and marginal parameter distributions for
    %plausible MH parameter values
    pairwise_fig=figure('Position', [1, 1, 1000, 1000]);
    t3=tiledlayout(1,1);
    nexttile(t3)
    [~,ax]=plotmatrix(plausible_param_values_scaled);
    for i=1:11
        ax(i,1).FontSize = 7; 
        ax(11,i).FontSize = 7; 
        ax(i,1).YLabel.String=param_names(i);
        ax(11,i).XLabel.String=param_names(i);
        ax(i,1).YLabel.FontSize=12;
        ax(11,i).XLabel.FontSize=12;
    end
    title(t3,'Plausible Pairwise and Marginal Parameter Distributions')
    xlabel(t3,{'','Relative Parameter Value (%)'})
    ylabel(t3,{'Relative Parameter Value (%)',''})
end