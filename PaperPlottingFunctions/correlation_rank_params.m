function correlation_rank_params(weight_scaled_param,weighted_responder_stat,param_names,layout)

    %Computes the pairwise spearman correlation between every parameter
    %pair and selects the most correlated (by abs) pairs to show pairwise
    %scatter plots that are colored by responder status. (Fig. 5b)

    %Input
    % weight_scaled_param -- parameter matrix of the very large virtual
    %                        population which has been sampled with
    %                        replacment. Parameters are scaled as % between
    %                        min and max value.
    % weighted_responder_stat -- responders status for the very large
    %                            virtual population
    % param_names -- string array of parameter names
    % layout -- layout handle for sub-plots in following functions

    %set observation and virtual population colors
    respndr_color = [0.4940, 0.1840, 0.5560];
    noresp_color = [0.9290, 0.6940, 0.1250];
    %create color maps for observed and virtual populations
    color_map_resp = linspace(.7,.1,10)'.*repmat(respndr_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    color_map_norsp = linspace(.7,.1,10)'.*repmat(noresp_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);

    %create wide table of each scaled parameter value, replicate weighted to
    %match probability of inclusion
    weighted_param_matrix = array2table(weight_scaled_param,'VariableNames',  param_names');
    %compute the absolute value of the Spearman correlation between each pair
    %of scaled parameters, we don't care about direction, only magnitude
    abs_corr = abs(corr(weighted_param_matrix{:,param_names},'Type','Spearman'));
    %zero the diagnoal elements of the correlation matrix (=1, hence abs_corr<1)
    %and zero the upper triangular elements of the correlation matrix
    %(these are replicates of the lower triangular elements which remain)
    abs_corr = tril(abs_corr.*(abs_corr<1));
    %generate matrix of parameter pair indices, used to look up top
    %correlated parameter pairs in above matrix
    mat_ind = [repmat([1:length(abs_corr)]',length(abs_corr),1) repelem([1:length(abs_corr)]',length(abs_corr),1)];
    %select the 4 largest correlation coef's (abs) from the linearized correlation matrix 
    [~,indx] = maxk(abs_corr(:),4);
    %use the index of the 4 largest coeffs to find the indicies of the
    %corresponding parameters pairs
    max_mat_ind = mat_ind(indx,:);
    %store the names of the top 4 most correlated parameter pairs
    cor_par_list = [param_names(max_mat_ind(:,1)) param_names(max_mat_ind(:,2))];

    %set up a incdices in layout to subplot for figure 5b   
    layout_indices = [3, 4, 7, 8];
    %loop over the top 4 parameters
    for i=1:4
        %select tile to plot on
        nexttile(layout,layout_indices(i))
        hold on
        %plot the responding patients in one color
        h1=scatter(weighted_param_matrix{weighted_responder_stat,cor_par_list(i,1)}, weighted_param_matrix{weighted_responder_stat,cor_par_list(i,2)},'filled','MarkerFaceColor',respndr_color,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.2);
        h1.SizeData = h1.SizeData/9;
        %plot the non-responding patients in another color
        h2=scatter(weighted_param_matrix{~weighted_responder_stat,cor_par_list(i,1)}, weighted_param_matrix{~weighted_responder_stat,cor_par_list(i,2)},'filled','MarkerFaceColor',noresp_color,'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.2);
        h2.SizeData = h2.SizeData/9;
        xlabel(param_names(max_mat_ind(i,1)))
        ylabel(param_names(max_mat_ind(i,2)))
        hold off
    end
    %create a legend
    hL1 = legend({'Non-responder','Responder'},'Location','bestoutside');

end