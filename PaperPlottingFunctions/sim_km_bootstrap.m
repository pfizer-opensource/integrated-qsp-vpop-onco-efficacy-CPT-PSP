function [km_uncertainty_fig,...
            km_bootstrap_selected_list] = sim_km_bootstrap(sample_size,...
                                                 num_replicate_trials,...
                                                 prob_inclusion,data,...
                                                 store_sims,plot_only,km_bootstrap_selected_list)

    %Generates a plot comparing bootstrap resampled uncertainty region for
    %KM estimate from the virtual population with the analytically derived
    %confidence region from the fit data (Figure 4b).

    %Input
    % sample_size -- sample size for bootrapped trial
    % num_replicate_trials -- number of replicated trials to use in
    %                         bootstrap computation of KM uncertainty
    % prob_inclusion -- probability of inclusiuon computed from get_vpop,
    %                   a normalized ratio of probility of seeing given
    %                   plausible patient w.r.t. clinical data divided by
    %                   KNN density estimate of plausible patients in
    %                   output space. Sums to 1.
    % data -- data structure containing clinical data
    % store_sims -- data structure for storing simulation of the plausible
    %               population

    %Output
    % km_uncertainty_fig -- figure handle for KM uncertainty plot

    %set observation and virtual population colors
    obspop_color = [0 0.4470 0.7410];
    vpop_color =[0.8500 0.3250 0.0980];
    
    %create color maps for observed and virtual populations
    color_map_obs = linspace(.7,.1,10)'.*repmat(obspop_color,10,1)+ ...
                    (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    color_map_vpop = linspace(.7,.1,10)'.*repmat(vpop_color,10,1)+...
                    (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    
    if(~plot_only)
        %extract BPC, DT and censoring labels for observed data
        best_change_data = sortrows([data.BESTPCHG'; 1:size(data.BESTPCHG)]',1,'descend');
        time_to_pfs_data = data.EEVALUMP(best_change_data(:,2));
        censoring_data = cell2mat(arrayfun(@(x)strcmp('NO PFS EVENT', data.EVENTNP_NEW{x}), 1:size(data,1),'UniformOutput', false))';
        censoring_data = censoring_data(best_change_data(:,2));
    end

    %km_bootstrap_selected_list={};

    %create arrays to store KM estamates from each bootstrap
    km_estimates =[];
    median_pfs = zeros(num_replicate_trials,1);
    %loop over number of replicate trials
    for i=1:num_replicate_trials
    
        %sample plausible indices with replacment, weighted by prob. inclusion
        if(~plot_only)
            selected = subsamp(sample_size,prob_inclusion);
            km_bootstrap_selected_list{i}=selected;
        else
            selected=km_bootstrap_selected_list{i};
        end
        %subset the stored simulations (replicate repeated plausibles, drop
        %those not selected)
        sample_store_sims = repelem(store_sims,1,selected);
        
        %retrieve the dropout times and censoring labels for the sample
        time_to_pfs_sim_t  = cell2mat(arrayfun(@(x) sample_store_sims(x).results{2}.time_to_pfs(1), 1:sample_size, 'UniformOutput', false))';
        censoring_sim_t  = cell2mat(arrayfun(@(x) sample_store_sims(x).results{2}.censor_label(1), 1:sample_size, 'UniformOutput', false))';
        
        %compute the KM estimate for the current sample
        [f_t,x_t]= ecdf(time_to_pfs_sim_t/7,'Function','survivor','Bounds','on','Censoring',censoring_sim_t);
        %store the x and y steps for the KM estimate of the given sample
        km_estimates(i).f_t=f_t;
        km_estimates(i).x_t=x_t;
        %find and store the median PFS for the given sample
        median_pfs(i) = km_estimates(i).x_t(cumsum(km_estimates(i).f_t<.5)==1);
        %stairs(x_t,f_t*100,'color',vpop_color,'LineWidth',2);
    end
    
    %vector of weeks to loop over
    week=0:62;
    %set perc such that 1-perc is the % for confidence intrerval over KM
    %survival probability
    perc=1;
    %arrayts to store median and upper/lower quantiles of bootstrapped KMs
    median_surv = [];
    upper = [];
    lower = [];
    %loop over weeks
    for w=week
        %array to store the interpolated survival probability from each KM curve
        surv_set = [];
        %loop over bootstrapped trial replicates
        for t=1:num_replicate_trials
            %interpolate the survival probability for each trial for the given week
            surv_p = interp1([0; km_estimates(t).x_t(2:end)],...
                    [1; km_estimates(t).f_t(2:end)],w,'previous');
            %save the survival probability for the given trial (on given week)
            surv_set = [surv_set surv_p];
        end
    
        %compute the median survival probability for the given week across
        %bootstrapped trials
        median_surv = [median_surv prctile(surv_set,50)];
        %compute the upper and lower quantiles for the confidence bounds
        upper = [upper prctile(surv_set,100-perc/2)];
        lower = [lower prctile(surv_set,perc/2)];
    end
    
    if(~plot_only)
        %plot Kaplan Meier plot computed on observed data
        [f,x,flo,fup] = ecdf(time_to_pfs_data/7,...
                            'Function','survivor',...
                            'Bounds','on','Censoring',censoring_data,...
                            'Alpha',0.01);
    end

    hold on
    %plot for the PFS/KM uncertainty plot
    km_uncertainty_fig=figure;
    %reshape lower and upper bounds for KM into a convex set of points for
    %plotting with fill() 
    x_stair = [week(1), repelem(week(2:end),2)];
    y_stair_low = [repelem(lower(1:end-1),2), lower(end)];
    y_stair_hgh = [repelem(upper(1:end-1),2), upper(end)];
    error_x = [x_stair, fliplr(x_stair)];
    error_y = [y_stair_low, fliplr(y_stair_hgh)]*100;
    fill(error_x,error_y,color_map_vpop(7,:),'LineStyle','none')
    
    hold on
    %add dotted lines to bounds of filled region for clearer visualization
    stairs(week,upper*100,':','color',vpop_color,'LineWidth',1.5,'HandleVisibility','off');
    stairs(week,lower*100,':','color',vpop_color,'LineWidth',1.5,'HandleVisibility','off');
    
    if(~plot_only)
        %format Greenwood confidence bounds from KM on original data into a convex
        %set of points for plotting with fill()
        x_stair = [x(2) ; repelem(x(3:end),2)];
        y_stair_low = [repelem(flo(2:end-1),2); flo(end)];
        y_stair_hgh = [repelem(fup(2:end-1),2); fup(end)];
        error_x = [x_stair; flipud(x_stair)];
        error_y = [y_stair_low; flipud(y_stair_hgh)]*100;
        fill(error_x,error_y,color_map_obs(7,:),'LineStyle','none','FaceAlpha',0.5,'HandleVisibility','off');%'Color',[map_obs(7,:) 0.2]
    end
    
    %plot interval for median PFS estimate from bootstrap resampled trials
    plot([prctile(median_pfs,perc/2) prctile(median_pfs,100-perc/2)],[50 50],'|-','Color',vpop_color*.85);
    
    if(~plot_only)
        stairs(x,f*100,'color',obspop_color,'LineWidth',2);
        xlabel('Time (weeks)'); ylabel('PFS Probability (%)')
        %compute observed median PFS from observed data
        median_pfs_data = x(cumsum(f<.5)==1);
        %plot a visual guide at 0.5 survival proability from 0 to median observed
        %PFS time
        %plot([0 median_pfs_data],[50,50],'--','LineWidth',1,'Color',obspop_color)
        %plot the median PFS with a vertical dashed line at the time
        plot([median_pfs_data median_pfs_data],[0,50],'--','LineWidth',1,'Color',obspop_color)
        %plot([0 70],[0 0],'k')
        yticks([0 25 50 75 100])
        xticks([15 30 45 60])
        xlim([0 60])
    end
    
    %add a label to the plot if desired
    %string = sprintf('Virtual Trial\n Median PFS Distribution');
    %text(median(median_pfs),-10,string,'FontSize',9,'HorizontalAlignment', 'center','VerticalAlignment','bottom','BackgroundColor',[1 1 1 0.7])
    
    %add legend
    if(~plot_only)
        legend({['99% Virtual Survival' newline 'Probability Interval by Week'],['99% Virtual Median' newline 'Survival Time Interval'],'Observed KM Estimate','Observed Median PFS'},'FontSize',7)
    else
        legend({['99% Virtual Survival' newline 'Probability Interval by Week'],['99% Virtual Median' newline 'Survival Time Interval']},'FontSize',7)
    end

end