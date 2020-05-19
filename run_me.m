%% Celeris master Work flow for Coastal Model Test Bed
%    This script can operate in hindcast mode by inputing start and end 
%    dates or in Nowcast mode modeling the most recent day
%    
%   INPUTS:
%       homeDir - this is where the celeris model exe and all
%           supplimental files are located.  must also have the yaml_files
%           directory to make netCDF files
%   
% set home working directory 
homeDir = "C:\Users\acoll\Desktop\celerisDataGen";
% wall clock time (in seconds)  should be 30, will save the last 17 minutes
sim_time=20*60; 

water_level_change = 0;
for wc_in=18:30
    jj=140;
    for ii=1:59
        cd(homeDir)
        main(sim_time, jj, wc_in)
        fprintf('Success');
        jj=jj+1
    close all
    end
end
tocl
exit