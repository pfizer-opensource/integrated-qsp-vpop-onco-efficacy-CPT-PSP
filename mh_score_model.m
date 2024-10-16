function [score,dropout_time_vp,censor_label_vp,constraints] = mh_score_model(results,state_bnds,schedule_index,screentime,vars_to_score,time_data,dropout_times_data,pd,mu,sigma,dropout_time_vp,only_constraints,censor_label_data)

% This function computes the MH score for a model simulation with respect
% to the Gaussian mixture fit to the clinical data. As part of computing
% the score a constraint check is performed, violations of which result in
% a score of zero.

% Input
% results -- results data structure from processing output of ODE simulation
% state_bnds -- bounds on the state variables
% schedule_index -- array of pk schedule index used in selection
% screentime -- prefix screening time
% vars_to_score -- variable names of outputs used in scoring (Guas. mixture) 
% time_data -- scanning times from data
% dropout_times_data -- list of possible dropout times from data
% pd -- probability distribution structure, fit to data for scoring
% mu -- means of Gaussian mixture
% sigma -- covariance of Gaussian mixture
% dropout_time_vp -- dropout time to score, if missing, score all times
% only_constraints -- flag indicating if function is being called for
%                     constraint check only
% censor_label_data -- censoring labels for each dropout time in data

% Output
% score -- score values for simulation (used in MH iteration accept/reject)
% dropout_time_vp -- if not passed as input, populated with posssible
%                    observed dropout times from data
% censor_label_vp -- censoring labels for each possible dropout time from data
% constraints -- logical check on subset of state bounds

uresults = results{1};
mdls = [uresults.shell; uresults.prct_necrotic_core; uresults.tum_doubling_time];
constraints = mdls(:)>=state_bnds(:,1) & mdls(:)<=state_bnds(:,2);
if only_constraints
    score = 0;
    censor_label_vp=0;
else
    pars = results{2}.pars;
    if ~isempty(dropout_time_vp)
        at = find(results{2}.t >= dropout_time_vp, 1, 'first');
        if isempty(at)
            score = eps;
        else
            tresults = structfun(@(x) x(1:at,:), rmfield(results{2},'pars'), 'UniformOutput', false);
            tresults = process_model_parameters_and_states(tresults.t, tresults.x, pars, schedule_index(2), screentime(2),time_data,[],[],dropout_time_vp);
            score = get_score(uresults,tresults,vars_to_score,constraints,pd,mu,sigma);
        end
    else
        % Erase tumour dynamics after PFS
        dropout_times_data(dropout_times_data > results{2}.time_to_pfs(1)) = [];
        censor_label_data(dropout_times_data > results{2}.time_to_pfs(1)) = [];
        %         tresults_all = results{2};
        tresults_t = results{2}.t;
        tresults_x = results{2}.x;
        si = schedule_index(2); st = screentime(2);
        parfor j = 1:length(dropout_times_data)
            tresults = process_model_parameters_and_states(tresults_t,tresults_x,pars,si,st,time_data,dropout_times_data(j),censor_label_data(j),dropout_times_data(j));
            score(j) = get_score(uresults,tresults,vars_to_score,constraints,pd,mu,sigma);
            dropout_time_vp(j) = dropout_times_data(j);
            censor_label_vp(j) = censor_label_data(j);
            best_dSLD(j) = tresults.best_dSLD;
            pfs_time(j) = tresults.time_to_pfs;
        end
        %%% troubleshooting - use this to plot PFS time vs BOR of all copies, 
        %%%                   colored according to their score
%         if any(score ~= eps)
%             surf([best_dSLD;best_dSLD], [pfs_time;pfs_time], ...
%                 [score;score], [score;score],...
%                 'EdgeColor','interp','LineWidth',2);
%             box on; view(2); colorbar; pause(1); drawnow;
%         end
    end
end
end % function mh_score_model

function score = get_score(uresults,tresults,vars,constraints,pd,mu,sigma)
if all(constraints)
    % Create the model observables:
    % In the M-H algorithm, higher is better for score, unlike other
    % algorithms:
    if ~isempty(pd)
        if length(vars) == 1
            vec = log(uresults.(vars{1}));
        else
            vec = [log(uresults.(vars{1})), cell2mat(arrayfun(@(x) tresults.(vars{x}), 2:length(vars),'UniformOutput',false))];
        end
        score = pdf(pd,vec);
    else
        if size(mu,2) == 1
            score = normpdf(log(uresults.(vars{1})),mu,sigma);
        else
            score = mvnpdf([log(uresults.(vars{1})), tresults.(vars{2}), tresults.(vars{3})],mu,sigma);
        end
    end
else
    score = eps; % Flag as a really, really small value
end
end
