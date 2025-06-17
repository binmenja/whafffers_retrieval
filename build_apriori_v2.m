clear all; close all; clc;

load('/home/binmenja/direct/field_campaigns/whafffers/retrievals/longterm_and_hourly_stats_era5_gault_feb.mat');

num_total_levels = size(p_levels, 1); % This should be 38 (surface + 37 pressure levels)
num_hourly_profiles = size(all_temp_hourly, 2); % Number of hourly observations

num_times_p1_lt_p2 = sum(all_p_hourly(1,:) < all_p_hourly(2,:));
percentage_p1_lt_p2 = (num_times_p1_lt_p2 / num_hourly_profiles) * 100;


Rd = 287.05; % Gas constant for dry air (J/kg/K)
g0 = 9.80665; % Standard gravity (m/s^2)
station_height_amsl = 129; % AMSL Gault station height in meters

all_gph_hourly = zeros(num_total_levels, num_hourly_profiles);

for j = 1:num_hourly_profiles

    current_p_profile = all_p_hourly(:, j); 
    current_t_profile = all_temp_hourly(:, j); 
    current_q_profile = all_q_hourly(:, j); %kg/kg

    % Calculate Tv
    current_tv_profile = current_t_profile .* (1 + 0.61 * current_q_profile);

    gph_for_this_profile = zeros(num_total_levels, 1);
    gph_for_this_profile(1) = station_height_amsl; % Surface level height (fixed)

    
    for i = 1:(num_total_levels - 1)
        P_i = current_p_profile(i);
        P_iplus1 = current_p_profile(i+1);

        
        if P_i <= P_iplus1
            
             warning('Pressure is not strictly decreasing from level %d to %d in profile %d. Adjusting P_i slightly.', i, i+1, j);
             P_i = P_iplus1 + 1e-6; % Make P_i infinitesimally larger than P_iplus1
        end

        % Mean virtual temperature in the layer
        Tv_mean_layer = (current_tv_profile(i) + current_tv_profile(i+1)) / 2;

        delta_gph = (Rd / g0) * Tv_mean_layer * log(P_i / P_iplus1);

        gph_for_this_profile(i+1) = gph_for_this_profile(i) + delta_gph;
    end
    
    all_gph_hourly(:, j) = gph_for_this_profile;
end


mean_gph_levels = mean(all_gph_hourly, 2); % Mean height for each of the 38 levels

data_priori.z = mean_gph_levels; % m, mean geopotential height for each of the 38 ERA5 levels
data_priori.t = temp_mean_longterm; % K, mean temperature on the 38 ERA5 levels
data_priori.q = q_mean_longterm .* 1000; % g/kg, mean specific humidity on the 38 ERA5 levels
data_priori.z_all = all_gph_hourly; % m, all hourly geopotential heights on their original 38 levels
data_priori.t_all = all_temp_hourly; % K, all hourly temperatures on their original 38 levels
data_priori.q_all = all_q_hourly .* 1000; % g/kg, all hourly specific humidity on their original 38 levels

save('era5_priori_data_v2.mat','data_priori');
disp('Apriori data saved to era5_priori_data_v2.mat');
disp(['Percentage of times surface pressure (p(1)) < first pressure level (p(2)): ', ...
      num2str(percentage_p1_lt_p2), '%']);