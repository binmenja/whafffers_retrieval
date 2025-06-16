function retrieval_results = retrieveFunction(utc_profile, path_modtran, v1, v2, resolution, fwhm, path_tape5_base, modroot, angle, emis, iLoc, ts, radiosonde_file, aeri_file, have_jacobian_ready, have_jacobian_iready, lambda_1st, is_q_log, dx_wv, dx_t, retrieval_type, priori_data_path, truth_profile_path)
    %   retrieval_results = retrieveFunction(...) executes an iterative retrieval
    %   process to derive atmospheric profiles (Temperature, Water Vapor, or both)
    %   from E-AERI observed radiance using a Levenberg-Marquardt optimization.
    %
    % Inputs:
    %   utc_profile: String, UTC time of radiosonde launch (e.g., '202502120753').
    %   path_modtran: String, path to MODTRAN executable.
    %   v1, v2: Numeric, spectral range lower and upper bounds (cm^-1).
    %   resolution: Numeric, spectral resolution (cm^-1).
    %   fwhm: Numeric, Full Width Half Maximum for spectral convolution.
    %   path_tape5_base: String, base path for MODTRAN TAPE5 files.
    %   modroot: String, root name for MODTRAN output files (e.g., 'whafffers_UTC_retrieval').
    %   angle: Numeric, viewing angle for MODTRAN (degrees).
    %   emis: Numeric, surface emissivity.
    %   iLoc: String, MODTRAN IESS parameter (location identifier).
    %   ts: Numeric, surface temperature (K).
    %   radiosonde_file: String, path to radiosonde NetCDF file.
    %   aeri_file: String, path to AERI NetCDF file.
    %   have_jacobian_ready: Logical (0 or 1), flag if Jacobians are pre-computed.
    %   have_jacobian_iready: Logical (0 or 1), flag if Jacobians are ready for the first iteration.
    %   lambda_1st: Numeric, initial Levenberg-Marquardt regularization parameter.
    %   is_q_log: Logical (0 or 1), flag if water vapor (q) is already in log scale.
    %   dx_wv: Numeric, water vapor uncertainty for convergence criterion (e.g., 0.1 for 10%).
    %   dx_t: Numeric, temperature uncertainty for convergence criterion (e.g., 0.5 for 0.5K).
    %   retrieval_type: Numeric (1: T only, 2: wv only, 3: T and wv).
    %   priori_data_path: String, path to a priori data .mat file (e.g., ERA5 data).
    %   truth_profile_path: String, path to truth profile .mat file (e.g., patched radiosonde).
    %
    % Outputs:
    %   retrieval_results: Struct containing all retrieval outputs and history.
    
    % Ensure output directories exist
    current_profile_dir = utc_profile;
    if ~exist(current_profile_dir, 'dir')
        mkdir(current_profile_dir);
    end
    if ~exist(path_tape5_base, 'dir')
        mkdir(path_tape5_base);
    end
    
    % Radiosonde data processing
    sonde_dt = ncread(radiosonde_file, 'DATETIME');
    ref_date_launch = datetime('2000-01-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC') + sonde_dt(1);
    
    % AERI data processing
    AERI_MOPD = ncread(aeri_file,'MOPD');
    AERI_fwhm = 1.2067/(2*AERI_MOPD);
    AERI_nesr_full = ncread(aeri_file, 'RADIANCE.SKY_NOISE');
    AERI_rad_full = ncread(aeri_file, 'RADIANCE.SKY_CLEANED');
    AERI_wnum = ncread(aeri_file, 'WAVENUMBER');
    AERI_dates_full_char = ncread(aeri_file, 'DATETIME');
    AERI_dates_full_str = string(AERI_dates_full_char.');
    
    years = double(extractBetween(AERI_dates_full_str, 1, 4));
    months = double(extractBetween(AERI_dates_full_str, 5, 6));
    days = double(extractBetween(AERI_dates_full_str, 7, 8));
    T_pos = strfind(AERI_dates_full_str, 'T');
    hour_aeri = zeros(size(AERI_dates_full_str));
    minute_aeri = zeros(size(AERI_dates_full_str));
    second_aeri = zeros(size(AERI_dates_full_str));
    
    for k = 1:length(AERI_dates_full_str)
        if ~isempty(T_pos{k})
            t_idx = T_pos{k};
            time_str = extractAfter(AERI_dates_full_str(k), t_idx);
            time_str = extractBefore(time_str, 'Z');
            hour_aeri(k) = str2double(extractBetween(time_str, 1, 2));
            minute_aeri(k) = str2double(extractBetween(time_str, 3, 4));
            second_aeri(k) = str2double(extractBetween(time_str, 5, 6));
        end
    end
    AERI_dates_full = datetime(years, months, days, hour_aeri', minute_aeri', second_aeri', 'TimeZone', 'UTC');
    
    time_diff = minutes(AERI_dates_full - ref_date_launch);
    start_idx = find(time_diff >= -2, 1, 'first');
    end_idx = find(time_diff <= 8, 1, 'last');
    if isempty(start_idx), start_idx = 1; end
    if isempty(end_idx), end_idx = length(AERI_dates_full); end
    
    AERI_rad = AERI_rad_full(:, start_idx:end_idx);
    AERI_nesr = AERI_nesr_full(:, start_idx:end_idx);
    
    % Retrieval type string
    switch retrieval_type
        case 1, variablename = 'T';
        case 2, variablename = 'wv';
        case 3, variablename = 'both';
    end
    
    % Read a priori data
    data_priori = importdata(priori_data_path);
    z_truth = (data_priori.z(:,1) + 12)./1000;
    t_prior = data_priori.t_all;
    q_prior = data_priori.q_all;
    
    x0_t = nanmean(t_prior,2);
    x0_q = nanmean(q_prior,2);
    nlev = length(z_truth);
    
    % Read truth data
    truth_profile_data = load(truth_profile_path);
    z_sonde = truth_profile_data.profile.z(:);
    p_sonde = truth_profile_data.profile.p(:);
    t_sonde = truth_profile_data.profile.t(:);
    q_sonde = truth_profile_data.profile.q(:);
    co2_sonde = truth_profile_data.profile.co2(:);
    ch4_sonde = truth_profile_data.profile.ch4(:);
    o3_sonde = truth_profile_data.profile.o3;
    
    t_truth = interp1(z_sonde,t_sonde,z_truth,'linear','extrap');
    q_truth = exp(interp1(z_sonde,log(q_sonde),z_truth,'linear','extrap'));
    q_truth(q_truth<0) = 1e-5;
    
    p = exp(interp1(z_sonde,log(p_sonde),z_truth,'linear','extrap'));
    co2 = exp(interp1(z_sonde,log(co2_sonde),z_truth,'linear','extrap'));
    ch4 = exp(interp1(z_sonde,log(ch4_sonde),z_truth,'linear','extrap'));
    o3 = exp(interp1(z_sonde,log(o3_sonde),z_truth,'linear','extrap'));
    
    % Prescribed vertically uniform GHGs
    co = zeros(1,nlev) + 0.1692;
    n2o = zeros(1,nlev) + 0.3171;
    
    if ~is_q_log
        q_prior = log(q_prior);
        q_truth = log(q_truth);
        x0_q = log(x0_q);
    end
    
    % A priori covariance Sa
    if strcmp(variablename, 'T')
        Sa = nancov(t_prior');
    elseif strcmp(variablename, 'wv')
        Sa = nancov(q_prior');
    elseif strcmp(variablename, 'both')
        combined_prior = [t_prior'; q_prior'];
        Sa = nancov(combined_prior);
    end
    
    % Initialize state vector x
    max_iterations = 20;
    if strcmp(variablename, 'both')
        x = NaN(2*nlev, max_iterations + 1);
    else
        x = NaN(nlev, max_iterations + 1);
    end
    
    if strcmp(variablename, 'T')
        xtrue = t_truth;
        xa = x0_t;
    elseif strcmp(variablename, 'wv')
        xtrue = q_truth;
        xa = x0_q;
    elseif strcmp(variablename, 'both')
        xtrue = [t_truth; q_truth];
        xa = [x0_t; x0_q];
    end
    
    x(:,1) = xa; % Initial guess
    
    % AERI measurement and Se (measurement error covariance)
    aeri_rad_mean = mean(AERI_rad,2,'omitnan');
    nesr_mean = mean(AERI_nesr,2,'omitnan');
    Se = diag(nesr_mean.^2);
    
    wv_begin = 600;
    wv_end = 1800;
    
    [~,index_aeri_begin] = find(abs(AERI_wnum' - wv_begin) == min(abs(AERI_wnum' - wv_begin)));
    [~,index_aeri_end] = find(abs(AERI_wnum' - wv_end) == min(abs(AERI_wnum' - wv_end)));
    AERI_wnum_adj = AERI_wnum(index_aeri_begin:index_aeri_end);
    
    measurement = aeri_rad_mean(index_aeri_begin:index_aeri_end);
    Se = Se(index_aeri_begin:index_aeri_end,index_aeri_begin:index_aeri_end);
    
    % Initialize cost function and related variables history
    J = NaN(1, max_iterations + 1);
    d = NaN(1, max_iterations);
    dy = NaN(1, max_iterations);
    lambda_history = NaN(1, max_iterations);
    F_output_history = NaN(length(AERI_wnum_adj), max_iterations);
    Spos_output_history = NaN(size(Sa,1), size(Sa,2), max_iterations);
    A_output_history = NaN(size(Sa,1), size(Sa,2), max_iterations);
    DFS_output_history = NaN(1, max_iterations);
    DFS_per_height_history = NaN(nlev, max_iterations);
    
    % Retrieval loop
    lambda = lambda_1st;
    F_prev_iter = []; % To store F from previous accepted iteration for dy calculation
    
    for i = 1:max_iterations
        disp(['ITERATION: ', num2str(i)]);
    
        % Prepare profile for forward model and Jacobian calculation
        profile_for_jac = struct();
        profile_for_jac.z = z_truth;
        profile_for_jac.p = p;
        profile_for_jac.co2 = co2;
        profile_for_jac.o3 = o3;
        profile_for_jac.ch4 = ch4;
        profile_for_jac.co = co;
    
        if strcmp(variablename, 'T')
            profile_for_jac.q = exp(q_truth);
            profile_for_jac.t = x(:,i);
        elseif strcmp(variablename, 'wv')
            profile_for_jac.q = exp(x(:,i));
            profile_for_jac.t = t_truth;
        elseif strcmp(variablename, 'both')
            profile_for_jac.q = exp(x(nlev+1:end,i));
            profile_for_jac.t = x(1:nlev,i);
        end
        
        % Jacobian calculation/loading
        jacobian_file_suffix = sprintf('x%d.mat', min(i-1,1)); % x0 for 1st iter, x1 for subsequent
        jacobian_path_t = fullfile(current_profile_dir, ['K_t_era5_', jacobian_file_suffix]);
        jacobian_path_q = fullfile(current_profile_dir, ['K_q_era5_', jacobian_file_suffix]);
    
        if ~have_jacobian_ready || (i == 1 && ~have_jacobian_iready) % Recompute if not ready or if it's the 1st iter and not using existing J for 1st
            compute_modtran_jacobian_temperature(utc_profile, profile_for_jac, jacobian_path_t);
            compute_modtran_jacobian_q(utc_profile, profile_for_jac, jacobian_path_q);
        end
        
        load(jacobian_path_t, 'jacobian_info');
        K_t = jacobian_info.jacobian .* 1e7; % Convert to RU/K
        sim_wnum = jacobian_info.wavenumbers;
    
        load(jacobian_path_q, 'jacobian_info');
        K_q = jacobian_info.jacobian .* 1e7; % Convert to RU/log(g/kg)
        
        % Handle potential NaNs in Jacobians by replacing with minimum absolute value
        K_t(isnan(K_t)) = min(abs(K_t(:)), [], 'all', 'omitnan');
        K_q(isnan(K_q)) = min(abs(K_q(:)), [], 'all', 'omitnan');
    
        % Forward simulation (F) for current state
        cloud_empty = struct('qi', [], 'ql', [], 'z', []);
        F_sim_output = run_single_simulation(path_modtran, path_tape5_base, ...
                profile_for_jac, modroot, resolution, fwhm, cloud_empty, ...
                v1, v2, angle, iLoc, emis, profile_for_jac.z(end), ts);
        F = F_sim_output.rad_plt .* 1e7; % Convert to RU
        
        % Adjust F and Jacobians to AERI resolution
        K_t_adj = zeros(size(AERI_wnum_adj, 1), size(K_t, 2));
        K_q_adj = zeros(size(AERI_wnum_adj, 1), size(K_q, 2));
        for il=1:size(K_t,2)
           K_t_adj(:,il) = band_conv_brb(sim_wnum, K_t(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
           K_q_adj(:,il) = band_conv_brb(sim_wnum, K_q(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
        end
        F = band_conv_brb(sim_wnum, F, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
        
        K_t = K_t_adj;
        K_q = K_q_adj;
        
        % Combine Jacobians based on retrieval type
        if strcmp(variablename, 'both')
            K = [K_t K_q];
        elseif strcmp(variablename,'T')
            K = K_t;
        else
            K = K_q;
        end
    
        F_output_history(:,i) = F;
    
        % Calculate cost function J for current iteration
        J(i) = (measurement - F)'*inv(Se)*(measurement - F) + (x(:,i) - xa)'*inv(Sa)*(x(:,i) - xa);
        
        J_candidate = J(i) + 1; % Initialize for while loop entry
    
        % Levenberg-Marquardt inner loop to find optimal lambda
        current_lambda = lambda; % Use a local lambda for this iteration's inner loop
        
        while J_candidate > J(i) && current_lambda < 1e30 % Loop until cost decreases or lambda too large
            % Check for NaNs and singular matrices
            if any(isnan(K(:))) || any(isnan(Se(:))) || any(isnan(Sa(:))) || any(isnan(measurement)) || any(isnan(F(:))) || any(isnan(x(:,i))) || any(isnan(xa))
                error('NaN values detected in inputs at iteration %d, sub-loop.', i);
            end
            if rcond(Se) < eps
                error('Se is singular at iteration %d, sub-loop. RCOND = %e', i, rcond(Se));
            end
    
            % Calculate Hessian approximation
            H = (1+current_lambda)*inv(Sa) + K'*inv(Se)*K;
            if rcond(H) < eps
                error('Hessian matrix H is singular at iteration %d, sub-loop. RCOND = %e', i, rcond(H));
            end
            
            % Calculate state vector update and candidate new state
            delta_x = H \ (K'*inv(Se)*(measurement - F) - inv(Sa)*(x(:,i)-xa));
            x_candidate = x(:,i) + delta_x;
    
            % Prepare profile for candidate forward simulation
            profile_for_F_candidate = struct('z', z_truth, 'p', p, 'co2', co2, 'ch4', ch4, 'o3', o3, 'co', co);
            if strcmp(variablename, 'T')
                profile_for_F_candidate.t = x_candidate;
                profile_for_F_candidate.q = exp(q_truth);
            elseif strcmp(variablename, 'wv')
                profile_for_F_candidate.t = t_truth;
                profile_for_F_candidate.q = exp(x_candidate);
            elseif strcmp(variablename, 'both')
                profile_for_F_candidate.t = x_candidate(1:nlev);
                profile_for_F_candidate.q = exp(x_candidate(nlev+1:end));
            end
    
            % Check for non-physical temperatures (if applicable for T retrieval)
            if any(profile_for_F_candidate.t <= 0) && (strcmp(variablename, 'T') || strcmp(variablename, 'both'))
                disp('Non-positive temperatures found in candidate profile. Increasing lambda.');
                J_candidate = J(i) + 1; % Force loop to continue with higher lambda
            else
                % Run forward simulation for candidate state
                F_candidate_sim_output = run_single_simulation(path_modtran, path_tape5_base, ...
                                        profile_for_F_candidate, modroot, resolution, fwhm, cloud_empty, ...
                                        v1, v2, angle, iLoc, emis, profile_for_F_candidate.z(end), ts);
                F_candidate = band_conv_brb(sim_wnum, F_candidate_sim_output.rad_plt .* 1e7, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
                
                % Calculate cost function for candidate state
                J_candidate = (measurement - F_candidate)'*inv(Se)*(measurement - F_candidate) + (x_candidate - xa)'*inv(Sa)*(x_candidate - xa);
            end
    
            % Adjust lambda
            if J_candidate > J(i)
                current_lambda = current_lambda * 10;
            else % Cost decreased
                x(:, i+1) = x_candidate; % Accept the new state
                F_accepted_this_iter = F_candidate; % Store F from this accepted step
                J(i+1) = J_candidate;
                if current_lambda >= 1
                    lambda_history(i) = current_lambda;
                    current_lambda = current_lambda / 10;
                else
                    lambda_history(i) = current_lambda;
                end
            end
        end % End of while loop for lambda adjustment
    
        % If cost didn't decrease in this outer iteration, break
        if J(i+1) > J(i)
            disp('Warning: Cost function did not decrease in this iteration. Breaking retrieval loop.');
            break;
        end
        
        F_output_history(:,i) = F_accepted_this_iter;
    
        % Calculate posterior covariance and averaging kernel
        H_final_iter = (1+lambda_history(i))*inv(Sa) + K'*inv(Se)*K;
        Spos_output_history(:,:,i) = inv(H_final_iter); % Posterior covariance
    
        A_classical = Spos_output_history(:,:,i) * K' * inv(Se) * K; % Averaging kernel
        DFS_classical = trace(A_classical);
        
        A_output_history(:,:,i) = A_classical;
        DFS_output_history(i) = DFS_classical;
        DFS_per_height_history(:,i) = diag(A_classical);
    
        % Calculate convergence metrics
        d(i) = norm(x(:,i+1) - x(:,i))^2; % L2 norm squared of state vector change
        if ~isempty(F_prev_iter)
            dy(i) = norm(F_accepted_this_iter - F_prev_iter)^2; % L2 norm squared of measurement vector change
        else
            dy(i) = 0; % No change for first iteration
        end
        F_prev_iter = F_accepted_this_iter; % Store for next iteration
    
        if retrieval_type == 1
            d_threshold = (zeros(1,size(x,1))+dx_t) * (pinv(Sa)+K'*pinv(Se)*K) * (zeros(size(x,1),1)+dx_t);
            %d_threshold(i) = (zeros(1,size(x,1))+dx_t) * pinv(Spos) * (zeros(size(x,1),1)+dx_t);
        elseif retrieval_type == 2
            %d_threshold(i) = (zeros(1,size(x,1))+log(1+dx_wv)) * pinv(Spos) * (zeros(size(x,1),1)+log(1+dx_wv));
            d_threshold = (zeros(1,size(x,1))+log(1+dx_wv)) * (pinv(Sa)+K'*pinv(Se)*K) * (zeros(size(x,1),1)+log(1+dx_wv));
        elseif retrieval_type == 3
            d_threshold = [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))] * (pinv(Sa)+K'*pinv(Se)*K) * [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))]';
            %d_threshold(i) = [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))] * pinv(Spos) * [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))]';
        end
        
        fprintf('Lambda: %.3e\n', lambda_history(i));
        fprintf('Cost function J(i): %.4f\n', J(i));
        fprintf('State step size d(i): %.4e\n', d(i));
        fprintf('Threshold d_threshold: %.4e\n', d_threshold);
        fprintf('DFS: %.2f\n', DFS_output_history(i));
        fprintf('Measurement step size dy(i): %.4e\n', dy(i));
    
        % Check for convergence
        if d(i) < d_threshold || i == max_iterations
            break;
        end
    end % End of main retrieval loop
    
    % Final results packaging
    retrieval_results = struct();
    retrieval_results.utc_profile = utc_profile;
    retrieval_results.retrieval_type = variablename;
    retrieval_results.iterations = i; % Actual number of performed iterations
    retrieval_results.x_retrieved = x(:, 1:i+1); % Final retrieved state and history
    retrieval_results.xa = xa;                  % A priori
    retrieval_results.xtrue = xtrue;            % Truth (from sonde)
    retrieval_results.Sa = Sa;                  % A priori covariance
    retrieval_results.Se = Se;                  % Measurement error covariance
    retrieval_results.AERI_wnum_adj = AERI_wnum_adj; % Adjusted AERI wavenumbers
    retrieval_results.measurement = measurement; % Mean AERI radiance measurement
    retrieval_results.Spos = Spos_output_history(:,:,i); % Final Posterior covariance
    retrieval_results.A = A_output_history(:,:,i);     % Final Averaging kernel
    retrieval_results.DFS = DFS_output_history(1:i);    % DFS history
    retrieval_results.DFS_per_height = DFS_per_height_history(:,1:i); % DFS per height history
    retrieval_results.J = J(1:i+1);             % Cost function history
    retrieval_results.d = d(1:i);               % State vector convergence history
    retrieval_results.dy = dy(1:i);             % Measurement space convergence history
    retrieval_results.lambda_history = lambda_history(1:i); % LM lambda history
    retrieval_results.F_output = F_output_history(:, 1:i); % Forward simulations history
    
    % Populate retrieved temperature and water vapor (physical units)
    if strcmp(variablename, 'both')
        retrieval_results.tx_retrieved = x(1:nlev, i+1);
        retrieval_results.qx_retrieved = exp(x(nlev+1:end, i+1));
    elseif strcmp(variablename, 'T')
        retrieval_results.tx_retrieved = x(:, i+1);
        retrieval_results.qx_retrieved = exp(q_truth); % q comes from truth if T-only retrieval
    elseif strcmp(variablename, 'wv')
        retrieval_results.tx_retrieved = t_truth; % T comes from truth if wv-only retrieval
        retrieval_results.qx_retrieved = exp(x(:, i+1));
    end
    
    % Save results to file
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    output_filename = sprintf('./%s/retrieval_results_%s_%s.mat', current_profile_dir, variablename, timestamp);
    save(output_filename, 'retrieval_results', '-v7.3');
    disp(['Saved retrieval results to: ', output_filename]);
    
    end
    