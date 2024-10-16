function [schedule_index,cc_Drug,TSTOP,INTERVAL] = load_pk_table(pk_file,schedule_flag,TSTOP,screentime,starttime,verbose_fun)

%This function loads a pk table from pk_table.xlsx and formats it
%appropriately so that is can be used for pk interpolation in the clinical
%ODE function simulation. A single median pk is used for all subjects in
%the virtual population.

%Input
% pk_file -- name of the pk_table file to read table from
% schedule_flag -- string indicating which pk schedule to use
% TSTOP -- stop time for trial, overwritten by table length if table longer
% screentime -- prefixed time for screening (legacy)
% starttime -- prefixed start time interval (legacy)
% verbose_fun -- function handle for providing verbose feedback/logging

%Output
% schedule_index -- array indicating an numeric index per schedule loaded
% cc_Drug -- drug concentration and time table
% TSTOP -- updated stop time, set by length of pk table if table>TSTOP
% INTERVAL -- time interval

num_schedules = length(schedule_flag);
schedule_index = nan(num_schedules,1);
cc_Drug = cell(num_schedules,1);
INTERVAL = nan(num_schedules,1);

for i = 1:num_schedules
    switch schedule_flag{i}
        
        case 'none'
            schedule_index(i) = 0;
            cc_Drug{i}   =[];
            INTERVAL(i) = 1;
        case 'Drug'
            schedule_index(i) = 1;
            
            %%% A1:B2502 to get data till 500 days.
            %%% A1:B2487 to go upto 497 days only.
            cc_Drug_temp = readtable(pk_file,...
                'Range','A1:B5902','Sheet',1,...
                'ReadVariableNames',true); %B2502
            end_ind = find(cc_Drug_temp.time == TSTOP(i)-(screentime(i) + starttime(i)),1,'first');
            if ~isempty(end_ind)
                time_new   = cc_Drug_temp.time(1:end_ind);
                cc_new     = cc_Drug_temp.Cc_Drug(1:end_ind);
            else
                time_new   = cc_Drug_temp.time;
                cc_new     = cc_Drug_temp.Cc_Drug;
            end
            
            if (screentime(i) > 0 || starttime(i) > 0)
                time_new   = [0; screentime(i) + starttime(i) + time_new];
                cc_new     = [0; cc_new];
            end
            cc_Drug{i}   =[time_new, cc_new];
            TSTOP(i) = max(TSTOP(i),max(time_new));
            INTERVAL(i) = 7;  %% (wks converted to days)
            
            feval(verbose_fun,'Computing results for ALKi 100 mg QD')
    end
end
end