function [z_setpoint, Vz_setpoint, csett] = set_Trajectory(z, Vz, csett)

% Author: Leonardo Bertelli
% Co-Author: Alessandro Del Duca
% Skyward Experimental Rocketry | ELC-SCS Dept | electronics@kywarder.eu
% email: leonardo.bertelli@skywarder.eu, alessandro.delduca@skywarder.eu
% Release date: 01/03/2021

if csett.iteration_flag == 1
    
    best_min   = inf;
    best_index = inf;

    for ind = 1:length(csett.data_trajectories)
       
        % Select a z trajectory and a Vz trajectory (to speed up select only the first values, not ALL)
        z_ref  = csett.data_trajectories(ind).Z_ref(1:150); 
        Vz_ref = csett.data_trajectories(ind).VZ_ref(1:150); 
        err    = (z_ref-z).^2 + (Vz_ref-Vz).^2; 

        % Find the nearest point to the current trajectory
        [min_value, index_min_value] = min( err ); 

        if (min_value < best_min)
             best_min = min_value;
             best_index = index_min_value;
             csett.chosen_trajectory = ind;  
             csett.Sbar = csett.data_trajectories(ind).Sbar*0.01/100;
        end

    end

    csett.index_min_value = best_index; % Save the actual index to speed up the research
    csett.iteration_flag  = 0; % Don't enter anymore the if condition
    
    % I select the reference altitude and the reference vertical velocity
    z_setpoint  =  csett.data_trajectories(csett.chosen_trajectory).Z_ref(csett.index_min_value);
    Vz_setpoint =  csett.data_trajectories(csett.chosen_trajectory).VZ_ref(csett.index_min_value);
    
    % Just for plot
    csett.starting_index = best_index;
    disp('trajectory chosen:')
    disp(csett.chosen_trajectory)

%% For the following iterations keep tracking the chosen trajectory
else

    % Select the z trajectory and the Vz trajectory 
    % To speed up the research, I reduce the vector at each iteration (add if-else for problem in index limits)
    z_ref  = csett.data_trajectories(csett.chosen_trajectory).Z_ref;  
    Vz_ref = csett.data_trajectories(csett.chosen_trajectory).VZ_ref;  

    % 1) Find the value of the altitude in z_reference nearer to z_misured 
    [~, current_index_min_value] = min(abs(z_ref(csett.index_min_value:end) - z));
    csett.index_min_value = csett.index_min_value + current_index_min_value -1; 
    
    % 2) Find the reference using Vz(z)
    z_setpoint  =  z_ref(csett.index_min_value);
    Vz_setpoint = Vz_ref(csett.index_min_value);
end  