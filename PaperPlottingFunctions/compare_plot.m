function compare_fig = compare_plot(store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTART,TSTOP,INTERVAL,pd,mu,sigma,plot_only,varargin)

    % This function plots a comparison between a selected virtual
    % population and the target clinical population used for selection.
    % These plots inlclude plots showing the overlap between the
    % virtual population and  the Gaussian mixture fit to the clinical data,
    % in each pair of dimensions of the score function. This function also 
    % comparisons via waterfall plots of BOR, median SLD time series and 
    % KM PFS estimates.
    
    % Input
    % store_sims -- data structure storing plausible patient simulations
    % v_pp -- parameter table for plausible population
    % v_names -- parameter names for plausible population
    % v_bnds -- parameter bounds for plausible population
    % schedule_flag -- pk index schedule flag, indicating which index to plot
    % screentime -- prefixed screening time for trial
    % TSTART -- start time
    % TSTOP -- stop time
    % INTERVAL -- time interval
    % pd -- distribution object used for scoring virtuals
    % mu -- means of gaussian mixture distribution object for scoring
    % sigma -- covariance of gaussian mixture distribution object for scoring
    % varargin-- variable length input argument list
    
    % Output
    % compare_fig -- figure handle for comparison figuer (Fig. 3)
    
    % initialise plotting related variables
    num_bins = 20;
    nrows = 2;
    ncols = 2;
    lw = 1;
    datacol = [1,0,0,0.9];
    simcol = 'b';
    
    p = inputParser;
    addRequired(p,'store_sims')
    addRequired(p,'v_pp')
    addRequired(p,'v_names')
    addRequired(p,'v_bnds')
    addRequired(p,'schedule_flag')
    addRequired(p,'screentime')
    addRequired(p,'TSTART')
    addRequired(p,'TSTOP')
    addRequired(p,'INTERVAL')
    addParameter(p, 'ln_transform',[1,0,0], @(x) isvector(x));
    if(~plot_only)
        addParameter(p,'data', readtable('./Data/synthetic_clinical_data.csv','TreatAsEmpty',{'NA'}))
    end
    addParameter(p,'select',true(size(store_sims,2),1));
    
    parse(p,store_sims,v_pp,v_names,v_bnds,schedule_flag,screentime,TSTOP,varargin{:})
    
    if(~plot_only)
        data = p.Results.data;
%         min_sld = min(data.SLD);
%         max_sld = max(data.SLD);
    end
    select      = p.Results.select;
    ln_transform = p.Results.ln_transform;
    
    store_sims_all = store_sims;
    
    store_sims = store_sims(select);
    

    %% obtain necessary model species from results
    uind = find(strcmp(schedule_flag, 'none'));
    tinds = find(~strcmp(schedule_flag, 'none'));
    num_treatments = length(tinds);
    
    num_vps = size(store_sims,2);
    
    if ~isempty(uind)
        baselineSLD = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.SLD_mm(1), 1:num_vps, 'UniformOutput', false));
        shell = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.shell(1), 1:num_vps, 'UniformOutput', false));
        prct_necrotic_core = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.prct_necrotic_core(1), 1:num_vps, 'UniformOutput', false));
        tum_doubling_time = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.tum_doubling_time(end), 1:num_vps, 'UniformOutput', false));
        
        prct_necrotic_core_ss = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.prct_necrotic_core_ss(1), 1:num_vps, 'UniformOutput', false));
        knecrotic = cell2mat(arrayfun(@(x) store_sims(x).results{uind}.knecrotic(1), 1:num_vps, 'UniformOutput', false));
    end
    
    dSLD_sim = cell(num_treatments,1);
    best_dSLD_sim = cell(num_treatments,1);
    time_to_best_sim = cell(num_treatments,1);
    time_to_pfs_sim = cell(num_treatments,1);
    
    if ~isempty(tinds)
        for i = 1:num_treatments
            dSLD_sim{i}         = cell2mat(arrayfun(@(x) store_sims(x).results{tinds(i)}.dSLD, 1:num_vps, 'UniformOutput', false))';
            best_dSLD_sim{i}    = cell2mat(arrayfun(@(x) store_sims(x).results{tinds(i)}.best_dSLD(1), 1:num_vps, 'UniformOutput', false))';
            time_to_best_sim{i} = cell2mat(arrayfun(@(x) store_sims(x).results{tinds(i)}.time_to_best(1), 1:num_vps, 'UniformOutput', false))';
            time_to_pfs_sim{i}  = cell2mat(arrayfun(@(x) store_sims(x).results{tinds(i)}.time_to_pfs(1), 1:num_vps, 'UniformOutput', false))';
            censoring_sim{i}  = cell2mat(arrayfun(@(x) store_sims(x).results{tinds(i)}.censor_label(1), 1:num_vps, 'UniformOutput', false))';
            
        end
    end
    
    
    %%%
    %% Check therapy response
    nrows = 2;
    ncols = 2;
    
    obspop_color = [0 0.4470 0.7410];
    vpop_color =[0.8500 0.3250 0.0980];
    map = linspace(.7,.1,10)'.*repmat(obspop_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    map_obs = linspace(.7,.1,10)'.*repmat(obspop_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);
    map_vpop = linspace(.7,.1,10)'.*repmat(vpop_color,10,1)+ (1-linspace(.7,.1,10)').*repmat([1 1 1],10,1);

    %load data
    if(~plot_only)
        time = cell2mat(arrayfun(@(x) str2double(cell2mat(regexp(data.Properties.VariableNames{x}, '\d*', 'match'))),...
            1:size(data,2),'UniformOutput', false));
        dSLD_data = data(:,{data.Properties.VariableNames{~isnan(time)}});
        dSLD_data = table2array(dSLD_data);
        num_vps_data = size(dSLD_data,1);
        
        time(isnan(time)) = [];
        
        best_change_data = sortrows([data.BESTPCHG'; 1:size(data.BESTPCHG)]',1,'descend');
        time_to_pfs_data = data.EEVALUMP(best_change_data(:,2));
        censoring_data = cell2mat(arrayfun(@(x)strcmp('NO PFS EVENT', data.EVENTNP_NEW{x}), 1:size(data,1),'UniformOutput', false))';
        censoring_data = censoring_data(best_change_data(:,2));
        best_change_data = best_change_data(:,1);
    end
    
    %%% Extract data and corresponding sim results
    if ~isempty(tinds)
        for t = 1:num_treatments
                        
            total_duration = TSTART:INTERVAL(tinds(t)):TSTOP(tinds(t));
            therapy_start_ind = find(total_duration == screentime(tinds(t)));
            time_sim = [0, (total_duration(therapy_start_ind+1:end)/7) - screentime(tinds(t))/7]';
            
            if(~plot_only)
                datainds_to_plot = 1:length(time);
                inds_to_plot = find(ismember(time_sim, time(datainds_to_plot)));
                max_time = max(time_sim(inds_to_plot));
            else
                datainds_to_plot = 1:13;
                inds_to_plot = find(ismember(time_sim, [0 6 12 18 24 30 36 42 48 54 60 66 72]));
                max_time = max(time_sim(inds_to_plot));
            end
            
            %CHANGED FOR COVERAGE PLOTS
            best_dSLD_sim_t = best_dSLD_sim{t};
            time_to_pfs_sim_t =  time_to_pfs_sim{t};
            censoring_sim_t = censoring_sim{t};
    
            best_dSLD_sim_t_sort = sortrows([best_dSLD_sim_t, [1:num_vps]'],1,'descend');
            time_to_pfs_sim_t_sort = time_to_pfs_sim_t(best_dSLD_sim_t_sort(:,2),:);

            ln_ind = find(ln_transform);
            alpha1 = pd.ComponentProportion;
            
            r=[baselineSLD',best_dSLD_sim_t(:,1),time_to_pfs_sim_t];
            r(:,ln_ind) = log(r(:,ln_ind));
            
            len_vars=3;
            total_plots = len_vars*(len_vars-1)/2;
            ncols = ceil(total_plots/2);
            [x,y] = meshgrid(1:len_vars, 1:len_vars);
            m = [x(:), y(:)]; n = m(m(:,2) < m(:,1),:);
            set(0, 'DefaultAxesFontSize',13)
            %figure()
    
            compare_fig = figure('Position', [1, 1, 1000, 800]);
            layout1 = tiledlayout(3,3);
%             layout2 = tiledlayout(layout1,3,1);
%             layout2.Layout.Tile = 1;
%             layout2.Layout.TileSpan = [3 1];
            data_names={'bSLD (mm)','BPC (%)','Dropout (weeks)'};
            %     figure(1)
            for i = 1:total_plots
                %subplot(3,12,[(i-1)*12+1:(i-1)*12+3])
                nexttile((i-1)*3+1)
                hold on;
                if i==1
                    [prob, X1Grid, X2Grid] = get_marginal(r,mu,sigma,alpha1,n(i,2),n(i,1));
                else
                    
                    [prob, X1Grid, X2Grid] = get_marginal(r,mu,sigma,alpha1,n(i,1),n(i,2));
                end
            
                lvls = max(max(prob))*logspace(0,-1,5);
                contourf(X1Grid, X2Grid, prob,lvls,'LineWidth',0.5);
    %             hold on
    %             contour(X1Grid, X2Grid, prob,lvls,'LineWidth',lw)
                colormap(map);
                %brighten(.5);
                
                tick_mm = [5,50,500];
                tick = log(tick_mm);
    
                tick_wks = [0,30,60];
                tick_days = tick_wks*7;
                
                box on;
                if i==1
                    scatter(r(:,n(i,2)),r(:,n(i,1)),70,vpop_color,'.');
                    xlabel(data_names{n(i,2)}); ylabel(data_names{n(i,1)});
                    if ln_transform(n(i,2))
                        set(gca,'XTick',tick,'XTickLabel',tick_mm);
                        xlim(log([5,500]))
                    elseif ln_transform(n(i,1))
                        set(gca,'YTick',tick,'YTickLabel',tick_mm);
                        ylim(log([5,500]))
                    end
                    if n(i,2)==3
                        set(gca,'XTick',tick_days,'XTickLabel',tick_wks);
                    elseif n(i,1)==3
                        set(gca,'YTick',tick_days,'YTickLabel',tick_wks);
                    end
                else
                    scatter(r(:,n(i,1)),r(:,n(i,2)),70,vpop_color,'.');
                    xlabel(data_names{n(i,1)}); ylabel(data_names{n(i,2)});
                    if ln_transform(n(i,2))
                        set(gca,'YTick',tick,'YTickLabel',tick_mm);
                        ylim(log([5,500]))
                    elseif ln_transform(n(i,1))
                        set(gca,'XTick',tick,'XTickLabel',tick_mm);
                        xlim(log([5,500]))
                    end
                    if n(i,2)==3
                        set(gca,'YTick',tick_days,'YTickLabel',tick_wks);
                    elseif n(i,1)==3
                        set(gca,'XTick',tick_days,'XTickLabel',tick_wks);
                    end
                end
    
                drawnow;
            end 
            
            %%%%%%%%%%%%%%%%% Plot results   %%%%%%%%%%%%%%%%%%
            
            %%%% WATERFALL PLOT
            %subplot(3,12,[6:12])
            nexttile(2,[1,2])
            hold on;
            if(~plot_only)
                s1 = stem(linspace(0,1,length(best_change_data))*100,best_change_data,'filled','Color',map_obs(5,:),'LineWidth',1,'MarkerSize',3,'Marker','none'); 
            end
            %alpha(s1,0.3);
            s2 = stem(linspace(0,1,num_vps)*100+0.2, best_dSLD_sim_t_sort(:,1),'filled','Color',map_vpop(5,:),'LineWidth',1,'MarkerSize',3,'Marker','none'); 
            %alpha(s2,0.3);
            if(~plot_only)
                scatter(linspace(0,1,length(best_change_data))*100,best_change_data,10,obspop_color,'filled','HandleVisibility','off','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
            end
            scatter(linspace(0,1,num_vps)*100+0.2, best_dSLD_sim_t_sort(:,1),10,vpop_color,'filled','HandleVisibility','off','MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0.7)
            xlim([0 100.03])
            xlabel 'Population Percentile (%)';  ylabel ({'Best % Change'})
            hold on; plot(linspace(0,100,5), ones(5,1)*-30,'--','LineWidth',1,'HandleVisibility','off','Color',"#77AC30");
            hold on; plot(linspace(0,100,5), ones(5,1)*20,'--','LineWidth',1,'HandleVisibility','off','Color',"#A2142F");
            set(gca,'LineWidth',lw,'box','on')
    
    
            %%%% SPIDER PLOT
            nexttile(5,[1,2])
            hold on;
            plot(time_sim(inds_to_plot), dSLD_sim{t}(:,inds_to_plot),'color',[vpop_color 0.1],'HandleVisibility','off');
            if(~plot_only)
                h = errorbar(time(datainds_to_plot),nanmean(dSLD_data(:,datainds_to_plot),1),nanstd(dSLD_data(:,datainds_to_plot),[],1),'color',obspop_color,'LineWidth',2);
            end
            errorbar(time_sim(inds_to_plot), nanmean(dSLD_sim{t}(:,inds_to_plot),1),nanstd(dSLD_sim{t}(:,inds_to_plot),[],1),'color',vpop_color,'LineWidth',2);
            set(gca,'ylim', [-100,100], 'xlim', [0, max_time],'LineWidth',lw,'box','on')
            hold on; plot(linspace(0,max_time,5), ones(5,1)*20,'--','LineWidth',1,'HandleVisibility','off','Color',"#A2142F");
            hold on; plot(linspace(0,max_time,5), ones(5,1)*-30,'--','LineWidth',1,'HandleVisibility','off','Color',"#77AC30");
            xlabel('Time (weeks)')
            ylabel({'% Change SLD'})
           
            
            %%%% PFS CURVE
            %subplot(3,12,[30,36])
            nexttile(8,[1,2])
            if(~plot_only)
                [f,x,flo,fup] = ecdf(time_to_pfs_data/7,'Function','survivor','Bounds','on','Censoring',censoring_data);
            end
            hold on
            [f_t,x_t,flo_t,fup_t]= ecdf(time_to_pfs_sim_t/7,'Function','survivor','Bounds','on','Censoring',censoring_sim_t);
            if(~plot_only)
                x_stair = [x(2) ; repelem(x(3:end),2)];
                y_stair_low = [repelem(flo(2:end-1),2); flo(end)];
                y_stair_hgh = [repelem(fup(2:end-1),2); fup(end)];
                error_x = [x_stair; flipud(x_stair)];
                error_y = [y_stair_low; flipud(y_stair_hgh)]*100;
                hold on;
                fill(error_x,error_y,map(7,:),'LineStyle','none','HandleVisibility','off')
                %hold on; ecdf(time_to_pfs_data/7,'Function','survivor','Bounds','on','Censoring',censoring_data);
                hold on; stairs(x,f*100,'color',obspop_color,'LineWidth',2);
            end
            hold on; stairs(x_t,f_t*100,'color',vpop_color,'LineWidth',2);
            if(~plot_only)
                x_cens_dat = unique(time_to_pfs_data(censoring_data==1)/7);
                x_cens_dat = x_cens_dat(~isnan(x_cens_dat));
                y_cens_dat = interp1(x(2:end),f(2:end)*100,x_cens_dat,'previous');
                plot(x_cens_dat,y_cens_dat,'+','Color',obspop_color*.8,'LineWidth',1.5)     
            end
    
            x_cens_sim = unique(time_to_pfs_sim_t(censoring_sim_t==1)/7);
            x_cens_sim = x_cens_sim(~isnan(x_cens_sim));
            y_cens_sim = interp1(x_t(2:end),f_t(2:end)*100,x_cens_sim,'previous');
            plot(x_cens_sim,y_cens_sim,'+','Color',vpop_color*.8,'LineWidth',1.5)
    
            set(gca,'LineWidth',lw,'box','on')
            xlabel('PFS (weeks)'); ylabel('Population %')

            hL = legend({'Observed Population','Virtual Population'},'Location','southoutside');
            %hL.Layout.Tile = 'South';

    
        end
        
    end



end

function [prob_marginal, X1Grid, X2Grid] = get_marginal(r,mu,sigma,alpha,i,j)
    if (i ~= j)
        [X1Grid, X2Grid] = meshgrid(linspace(min(r(:,i))-.1*range(r(:,i)),max(r(:,i)*1.1)+.1*range(r(:,i)),100),linspace(min(r(:,j))-.1*range(r(:,j)),max(r(:,j))+.1*range(r(:,j)),100));%sortrows(r(:,i)), sortrows(r(:,j))
        pd = gmdistribution(mu(:,[i,j]),sigma([i,j],[i,j],:),alpha);
        prob_marginal = reshape(arrayfun(@(x,y) pdf(pd, [x,y]), X1Grid(:), X2Grid(:)), size(X1Grid,1), size(X1Grid,2));
    else
        X1Grid = linspace(min(r(:,i))-.1*range(r(:,i)),max(r(:,i)*1.1)+.1*range(r(:,i)),100);
        pd = gmdistribution(mu(:,i),sigma(i,i,:),alpha);
        prob_marginal = arrayfun(@(x) pdf(pd, x), X1Grid(:));
        X2Grid = [];
    end
end