function [] = run_multiple_modtran_retrievals(datestring,n_recalc_jacob)

    % Created by Benjamin Riot--Bretecher (2025). 
    addpath('/home/binmenja/direct/matlab/mylib')
    addpath('/home/binmenja/direct/models/modtran6/bin/linux')
    addpath('/home/binmenja/direct/models/modtran6/bin/linux/matlib_modtran')
    addpath('/lustre03/project/6003571/binmenja/field_campaigns/whafffers/retrievals')

    if ~exist(datestring,'dir')
        system(strcat("mkdir -p ", datestring));
    end

    if nargin < 2
        n_recalc_jacob = 2; % default jacobian recalculation
    end

    % Retrieval settings
    have_jacobian_ready = 0; % whether the jacobian is already ready
    fprintf('Have jacobian ready: %d\n', have_jacobian_ready);
    lambda_1st = 10000;
    is_q_log = 0; % whether q is already in log scale; 0: no, 1: yes
    dx_wv = 0.1; % 10 percent uncertainty for wv
    dx_t = 0.5; % 0.5K uncertainty for T

    cloud.qi = []; cloud.ql = []; cloud.z = []; % clear-sky


    fprintf('Is q in log scale: %d\n', is_q_log);
    path_modtran = '/home/binmenja/direct/models/modtran6/bin/linux/'; v1 = 520; v2 = 1800; resolution = 0.1; fwhm = resolution*2;
    path_tape5 = ['/home/binmenja/direct/field_campaigns/whafffers/tape5/', datestring, '/'];
    if ~exist(path_tape5,'dir')
        system(["mkdir -p ", path_tape5]);
    end
    modroot = ['whafffers_', datestring, '_retrieval'];
    angle = 180; emis = 0; iLoc = '0'; ts = 0;

    % Extract datetime from profile name (datestring input)
    profile_dt = datetime(datestring, 'InputFormat', 'yyyyMMddHHmm', 'TimeZone', 'UTC');
    
    sonde_dir = '/home/binmenja/direct/field_campaigns/whafffers/radiosonde_data/';
    sonde_files = dir(fullfile(sonde_dir, 'balloon_sonde_mcgill_gault_*.nc'));
    sonde_match = '';

    for k = 1:length(sonde_files)
        fname = sonde_files(k).name;
        % Extract launch time portion from filename
        tokens = regexp(fname, '_(\d{8})T(\d{6})Z_', 'tokens');
        if ~isempty(tokens)
            launch_str = strcat(tokens{1}{1}, tokens{1}{2}(1:4));  % YYYYMMDDHHMM
            fprintf('Checking radiosonde file: %s (launch: %s)\n', fname, launch_str);
            if strcmp(launch_str, datestring)
                sonde_match = fullfile(sonde_dir, fname);
                break;
            end
        end
    end

    if isempty(sonde_match)
        error('No matching radiosonde file found for %s', datestring);
    else
        fprintf('Matched radiosonde file: %s\n', sonde_match);
    end

    radiosonde_file = sonde_match;
    
    % Extract launch datetime from radiosonde NetCDF
    sonde_dt = ncread(radiosonde_file, 'DATETIME'); % Format days since 2000-01-01 00:00:00 (MJD2K)
    ref_date_launch = datetime('2000-01-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC') + sonde_dt(1); % Convert to datetime
    disp(ref_date_launch)
    disp(['Radiosonde launch time: ', datestr(ref_date_launch)]);

     % === Find matching AERI file (by date only) ===
     aeri_dir = '/home/binmenja/direct/field_campaigns/whafffers/aeri_data/gault/';
     aeri_files = dir(fullfile(aeri_dir, 'atmospheric.physics_aeri_125_MCGILL_gault_*.nc'));
     aeri_match = '';
     
     for k = 1:length(aeri_files)
         fname = aeri_files(k).name;
         if contains(fname, datestr(profile_dt, 'yyyymmdd'))
             aeri_match = fullfile(aeri_dir, fname);
             break;
         end
     end
     
     if isempty(aeri_match)
         error('No matching AERI file found for %s', datestring);
     else
         fprintf('Matched AERI file: %s\n', aeri_match);
     end
     aeri_file = aeri_match;
    AERI_MOPD = ncread(aeri_file,'MOPD'); % AERI MOPD
    AERI_fwhm = 1.2067/(2*AERI_MOPD); % AERI FWHM
    AERI_cal_error_full = ncread(aeri_file, 'RADIANCE.SKY_ERROR');
    AERI_nesr_full = ncread(aeri_file, 'RADIANCE.SKY_NOISE'); % [RU]
    AERI_rad_full = ncread(aeri_file, 'RADIANCE.SKY_CLEANED'); % [RU]
    AERI_wnum = ncread(aeri_file, 'WAVENUMBER'); %[cm-1]
    AERI_dates_full = ncread(aeri_file, 'DATETIME'); % Format is YYYYMMDDThhmmssZ
    AERI_dates_full = string(AERI_dates_full.');
    years = double(extractBetween(AERI_dates_full, 1, 4));
    months = double(extractBetween(AERI_dates_full, 5, 6));
    days = double(extractBetween(AERI_dates_full, 7, 8));
    T_pos = strfind(AERI_dates_full, 'T');
    for i = 1:length(AERI_dates_full)
        if ~isempty(T_pos{i})
            t_idx = T_pos{i};
            time_str = extractAfter(AERI_dates_full(i), t_idx);
            time_str = extractBefore(time_str, 'Z'); %remove z
            hour_aeri(i) = str2double(extractBetween(time_str, 1, 2));
            minute_aeri(i) = str2double(extractBetween(time_str, 3, 4));
            second_aeri(i) = str2double(extractBetween(time_str, 5, 6));
        end
    end
    AERI_dates_full = datetime(years, months, days, hour_aeri', minute_aeri', second_aeri', 'TimeZone', 'UTC');
    time_diff_duration = AERI_dates_full - ref_date_launch; 
    time_diff = minutes(time_diff_duration); 
    start_idx = find(time_diff >= -2, 1, 'first'); 
    end_idx = find(time_diff <= 8, 1, 'last');    
    if isempty(start_idx), start_idx = 1; end
    if isempty(end_idx), end_idx = length(AERI_dates_full); end
    AERI_rad = AERI_rad_full(:, start_idx:end_idx);
    AERI_cal_error = AERI_cal_error_full(:, start_idx:end_idx); % 1-sigma absolute error
    AERI_nesr = AERI_nesr_full(:, start_idx:end_idx); % 1-sigma noise eq. spe. rad.

    output_dir = fullfile('output_retrievals_modtran/38levels', datestring);
    if ~exist(output_dir, 'dir'); mkdir(output_dir); end
    cd(output_dir);

    % retrieval type: 1(T only), 2(wv only), 3(T and wv)
    retrieval_type = 3; 
    switch retrieval_type
    case 1
        variablename = 'T';
    case 2
        variablename = 'wv';
    case 3
        variablename = 'both';
    end

    % read priori
    data_priori = importdata('/home/binmenja/direct/field_campaigns/whafffers/retrievals/era5_priori_data_v2.mat'); 
    z_truth = data_priori.z./1000; % km
    z_prior = data_priori.z./1000; % km
    t_prior = data_priori.t_all; % K
    q_prior = data_priori.q_all; % g/kg
    x0_t = nanmean(t_prior,2);
    x0_q = nanmean(q_prior,2); 
    nlev = length(z_truth);

    % read truth
    profile_file = fullfile('/home/binmenja/direct/field_campaigns/whafffers/200levels/', ...
        ['profile_' datestring 'UTC.mat']);
    if ~exist(profile_file, 'file')
    error('Profile file not found: %s', profile_file);
    end
    load(profile_file);      
    z_sonde = profile.z(:) + z_truth(1); % km; adjust to same base height!
    p_sonde = profile.p(:); % hPa
    t_sonde = profile.t(:); % K
    q_sonde = profile.q(:); % g/kg
    co2_sonde = profile.co2(:); % ppmv
    ch4_sonde = profile.ch4(:); % ppmv
    o3_sonde = profile.o3; % g/kg

    t_truth = interp1(z_sonde,t_sonde,z_truth,'linear','extrap');
    q_truth = exp(interp1(z_sonde,log(q_sonde),z_truth,'linear','extrap'));
    q_truth(q_truth<0) = 1e-5;

    p = exp(interp1(z_sonde,log(p_sonde),z_truth,'linear','extrap'));
    co2 = exp(interp1(z_sonde,log(co2_sonde),z_truth,'linear','extrap'));
    ch4 = exp(interp1(z_sonde,log(ch4_sonde),z_truth,'linear','extrap'));
    o3 = exp(interp1(z_sonde,log(o3_sonde),z_truth,'linear','extrap'));

    %% prescribed vertically uniform ghg, not used yet
    co = zeros(1,nlev) + 0.1692; % ppmv, vertically uniform
    n2o = zeros(1,nlev) + 0.3171; % ppmv, vertically uniform
    ccl4 = 1.0012e-4; % ppmv, vertically uniform
    f11 = 2.6009e-4; % ppmv, vertically uniform
    f12 = 5.4520e-4;% ppmv, vertically uniform

    % Check if z_truth, z_sonde and z_prior are the same
    if ~isequal(z_truth, z_sonde) || ~isequal(z_truth, z_prior)
        disp('z_truth, z_sonde and z_prior must be the same');
        fprintf('z_truth 1: %s\n', mat2str(z_truth(1)));
        fprintf('z_sonde 1: %s\n', mat2str(z_sonde(1)));
        fprintf('z_prior 1: %s\n', mat2str(z_prior(1)));
    end


    if ~is_q_log % enters if q is not log scale
        q_prior = log(q_prior);
        q_truth = log(q_truth);
        x0_q = log(x0_q);
        disp('q converted to log scale');    
    else
        disp('q already in log scale');
    end

    % Sa from ERA5
    if strcmp(variablename, 'T')
        Sa = nancov(t_prior'); % covariance matrix of T
    elseif  strcmp(variablename, 'wv')
        Sa = nancov(q_prior'); % covariance matrix of q
    elseif  strcmp(variablename, 'both')
        combined_prior = [t_prior', q_prior']; % combine T and q
        size(combined_prior)
        Sa = nancov(combined_prior); % covariance matrix of T and q
        size(Sa)
    end

    if strcmp(variablename, 'both')
        x = NaN(2*nlev, 0);
    else
        x = NaN(nlev, 0);
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

    x(:,1) = xa; % initial guess

    % aeri measurement and Se prescription
    aeri_rad_mean = mean(AERI_rad,2,'omitnan'); % RU
    nesr_mean = mean(AERI_nesr,2,'omitnan');
    Se = diag(nesr_mean.^2);

    wv_begin = 600% cm-1, lower bound of aeri wavenumber
    wv_end = 1800;  % cm-1, upper bound of aeri wavenumber
    
    [~,index_aeri_begin] = find(abs(AERI_wnum' - wv_begin) == min(abs(AERI_wnum' - wv_begin)));
    [~,index_aeri_end] = find(abs(AERI_wnum' - wv_end) == min(abs(AERI_wnum' - wv_end)));
    AERI_wnum_adj = AERI_wnum(index_aeri_begin:index_aeri_end); % cm-1

    measurement = aeri_rad_mean(index_aeri_begin:index_aeri_end);
    Se = Se(index_aeri_begin:index_aeri_end,index_aeri_begin:index_aeri_end); % RU^2 [chnl x chnl]

    %% initialize cost function
    clear J
    clear d

    profile_input.co2 = co2; % ppmv
    profile_input.o3 = o3; % g/kg
    profile_input.ch4 = ch4; % ppmv
    profile_input.co = co; % ppmv
    profile_input.z = z_truth; % km
    profile_input.p = p; % hPa
    ts = 0; % space temperature


    %% Retrieval
    lambda = lambda_1st; %for Levenberg-Marquardt regularization
    lambda_output = NaN(1, 20); 
    T = zeros(length(Sa),length(Se));
    for i = 1:20


        tic
        disp(['ITERATION: ', num2str(i)]);
        lambda_output(i) = lambda;
        if strcmp(variablename, 'T')
            qx = exp(q_truth); % This one is kept if q is not retrieved
            tx = x(:,i);
        elseif strcmp(variablename, 'wv')
            qx = exp(x(:,i));
            tx = t_truth; % This one is kept if T is not retrieved
        elseif strcmp(variablename, 'both')
            qx = exp(x(nlev+1:end,i));
            tx = x(1:nlev,i);
        end

        % Reassign profile_input for q and T
        profile_input.t = tx; % K
        profile_input.q = qx; % g/kg (not log anymore)
        pwd
        jacobian_folder = 'jacobians';
        jacobian_dir = fullfile('.', jacobian_folder);
        if ~exist(jacobian_dir, 'dir')
            mkdir(jacobian_dir);
            fprintf('Created directory: %s\n', fullfile(pwd, jacobian_folder));
        end

        jacobian_suffix = sprintf('x%d', i-1); % x0 for i=1, x1 for i=2, etc.
        jacobian_path_t = fullfile('./', datestring, ['K_t_era5_', jacobian_suffix, '.mat']); 
        jacobian_path_q = fullfile('./', datestring, ['K_q_era5_', jacobian_suffix, '.mat']);

        if have_jacobian_ready
            prev_suffix = sprintf('x%d', n_recalc_jacob-1);  % Always reuse latest jacobian
            jacobian_path_t = fullfile('./', datestring, ['K_t_era5_', prev_suffix, '.mat']);
            jacobian_path_q = fullfile('./', datestring, ['K_q_era5_', prev_suffix, '.mat']);
            if exist(jacobian_path_t, 'file') && exist(jacobian_path_q, 'file')
                load(jacobian_path_t); 
                K_t_ori = jacobian_info.jacobian .* 1e7;
                load(jacobian_path_q);
                K_q_ori = jacobian_info.jacobian .* 1e7;
                sim_wnum = jacobian_info.wavenumbers;
            else
                error('Jacobian files not found for iteration %d. Expected: %s and %s', ...
                    i, jacobian_path_t, jacobian_path_q);
            end
        elseif i <= n_recalc_jacob
            % Update input profile only from i >= 2 onward
            if i > 1
                profile_input.q = exp(x(nlev+1:end, i)); % g/kg
                profile_input.t = x(1:nlev, i);           % K
            end

            % Compute and save Jacobians
            compute_modtran_jacobian_temperature(datestring, profile_input, jacobian_path_t);
            compute_modtran_jacobian_q(datestring, profile_input, jacobian_path_q);

            % Load the newly computed jacobians
            load(jacobian_path_t); % loads variable: jacobian_info
            K_t_ori = jacobian_info.jacobian .* 1e7; % Units W/(cm^2*sr*cm^{-1})/K
            load(jacobian_path_q);
            K_q_ori = jacobian_info.jacobian .* 1e7; % Units W/(cm^2*sr*cm^{-1})/log(g/kg)
            sim_wnum = jacobian_info.wavenumbers;

            % After 2nd iteration, stop recomputation
            if i == n_recalc_jacob
                have_jacobian_ready = 1;
            end
        end

        % Convolve to AERI resolution
        for il = 1:size(K_t_ori,2)
            K_t(:,il) = band_conv_brb(sim_wnum, K_t_ori(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
            K_q(:,il) = band_conv_brb(sim_wnum, K_q_ori(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
        end

        % FORWARD SIMULATION
        if i == 1
            tic
            F = run_single_simulation(path_modtran, path_tape5, ...
                profile_input, modroot, resolution, fwhm, cloud, ...
                v1, v2, angle, iLoc, emis, profile.z(end), ts, true);
            F = F.rad_plt .* 1e7; % Convert to RU
        
            % Fix negative or NaN values
            bad_idx = isnan(F) | F < 0;
            if any(bad_idx)
                fprintf('[Iteration %d] Detected NaN or negative radiance values. Fixing...\n', i);
                F(bad_idx) = eps;
            end
        
            % Convolve to AERI resolution once
            F = band_conv_brb(sim_wnum, F, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
            
            F_new = F; % Cache forward model result
            toc
        else
            F = F_new; % Use cached result
        end
        
        cleanup_modtran_files(path_modtran, modroot); disp('Trying to clean up in current directory as well...');
        current_path = pwd; current_path = strcat(current_path, '/');cleanup_modtran_files(current_path, modroot);

        if strcmp(variablename, 'both')
            K = [K_t K_q]; % dimension: (nwv,nlev x 2)
        elseif strcmp(variablename,'T')
            K = K_t; % dimension: (nwv,nlev)
        else 
            K = K_q; % dimension: (nwv,nlev)
        end

        F_output(:,i) = F;

        J(i) = (measurement - F)'*inv(Se)*(measurement - F) + (x(:,i) - xa)'*inv(Sa)*(x(:,i) - xa);

        J(i+1) = J(i)+1;

        while J(i+1)>J(i)
            % disp(['Iteration ', num2str(i)]);
            % disp(['Size of K: ', mat2str(size(K))]);
            % disp(['Size of Se: ', mat2str(size(Se))]);
            % disp(['Size of Sa: ', mat2str(size(Sa))]);
            % disp(['Size of measurement: ', mat2str(size(measurement))]);
            % disp(['Size of F: ', mat2str(size(F))]);
            % disp(['Size of x(:,i): ', mat2str(size(x(:,i)))]);
            % disp(['Size of xa: ', mat2str(size(xa))]);
            
            x(:, i+1) = x(:,i) + inv((1+lambda)*inv(Sa) + K'*inv(Se)*K)*(K'*inv(Se)*(measurement - F)-inv(Sa)*(x(:,i)-xa));

            if strcmp(variablename, 'T')
                tx = x(:,i+1);
            elseif strcmp(variablename, 'wv')
                qx = exp(x(:,i+1));
            elseif strcmp(variablename, 'both')
                qx = exp(x(nlev+1:end,i+1)); % back to g/kg
                tx = x(1:nlev,i+1);
            end
            if any(tx <= 0)
                disp('Non-positive temperatures found:');
                disp(find(tx <= 0));
                J(i+1) = NaN;
            else
                profile_input.t = tx; % K
                profile_input.q = qx; % g/kg
                tic
                F_new = run_single_simulation(path_modtran, path_tape5, ...
                profile_input, modroot, resolution, fwhm, cloud, ...
                v1, v2, angle, iLoc, emis, profile.z(end), ts, true);
                F_new = F_new.rad_plt .* 1e7;
                toc
                F_new = band_conv_brb(sim_wnum, F_new, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');

                % Fix negative or NaN values
                bad_idx = isnan(F) | F < 0;
                if any(bad_idx)
                    fprintf('[Iteration %d] Detected NaN or negative radiance values. Fixing...\n', i);
                    F(bad_idx) = eps;
                end
                
                J(i+1) = (measurement - F_new)'*inv(Se)*(measurement - F_new) + (x(:,i+1) - xa)'*inv(Sa)*(x(:,i+1) - xa);

                fprintf('Forward simulation done for iteration %d.\n', i+1);

                
                disp('Trying to clean up in current directory...');
                cleanup_modtran_files(path_modtran, modroot);
                disp('Trying to clean up in current directory as well...');
                current_path = pwd;
                current_path = strcat(current_path, '/');
                cleanup_modtran_files(current_path, modroot);
            
            
            end
            if isnan(J(i+1))
                J(i+1) = J(i)+1;
            end

            if J(i+1)>J(i) 
                if lambda<1e30
                    lambda = lambda*10;
                else
                    break
                end
            elseif J(i+1)<J(i)
                if lambda>=1
                    lambda_output(i) = lambda;
                    lambda = lambda/10;
                else
                    lambda_output(i) = lambda;
                end
            end
            
        end

        % Lambda method, lambda not zero
        Spos = inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K)*((lambda_output(i)+1).^2*inv(Sa)+K'*inv(Se)*K)*inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K);
        Spos_output(:,:,i) = Spos;
        A = inv(K'*inv(Se)*K+(lambda_output(i)+1).*inv(Sa))*K'*inv(Se)*K; 
        A_output(:,:,i) = A;
        DFS_output(i) = trace(A); % Degrees of freedom for signal
        DFS_per_height(:,i) = diag(A); % Degrees of freedom per height


        % Fair version, lambda = 0
        G = (K' * inv(Se) * K + inv(Sa)) \ (K' * inv(Se)); 
        Spos_CR = inv(K'*inv(Se)*K+inv(Sa));
        Spos_output_CR(:,:,i) = Spos_CR;
        A_CR = Spos_CR * K' * inv(Se) * K;
        A_output_CR(:,:,i) = A_CR;
        DFS_output_CR(i) = trace(A_CR); % Degrees of freedom for signal (Rodgers full form)
        DFS_per_height_CR(:,i) = diag(A_CR); % Degrees of freedom per height (Rodgers full form)
        DFS_T_output_CR(i) = trace(A_CR(1:nlev, 1:nlev)); % DFS for T (Rodgers full form)
        DFS_Q_output_CR(i) = trace(A_CR(nlev+1:end, nlev+1:end)); % DFS for q (Rodgers full form)

        %dx^2 threshold
        if retrieval_type == 1
            d_threshold = (zeros(1,size(x,1))+dx_t) * (pinv(Sa)+K'*pinv(Se)*K) * (zeros(size(x,1),1)+dx_t);
        elseif retrieval_type == 2
            d_threshold = (zeros(1,size(x,1))+log(1+dx_wv)) * (pinv(Sa)+K'*pinv(Se)*K) * (zeros(size(x,1),1)+log(1+dx_wv));
        elseif retrieval_type == 3
            d_threshold = [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))] * (pinv(Sa)+K'*pinv(Se)*K) * [zeros(1,size(x,1)./2)+dx_t (zeros(1,size(x,1)./2)+log(1+dx_wv))]';
        end

        
        d(i) = (x(:,i)-x(:,i+1))' * (pinv(Sa)+K'*pinv(Se)*K) * (x(:,i)-x(:,i+1));
        Sy = Se * inv(K*Sa*K'+Se) * Se;
        dy(i) = (F-F_new)'*inv(Sy)*(F-F_new);
        drad(i,:) = measurement - F_new; % RU

        fprintf('Lambda: %.3e\n', lambda_output(i));
        fprintf('Cost function J(i): %.4f\n', J(i));
        fprintf('Step size d(i): %.4e\n', d(i));
        fprintf('Threshold d_threshold: %.4e\n', d_threshold);
        fprintf('d value: %.4e\n', d(i));
        fprintf('DFS (standard): %.2f\n', DFS_output(i));
        fprintf('DFS for T (lambda equals zero): %.2f\n', DFS_T_output_CR(i));
        fprintf('DFS for q (lambda equals zero): %.2f\n', DFS_Q_output_CR(i));
        fprintf('Measurement mismatch dy(i): %.4e\n', dy(i));

        if  d(i) < min(d_threshold,length(Sa)./20)
            break
        end
        toc
    end

    retrieval_results = struct();

    retrieval_results.utc_profile = datestring;
    retrieval_results.drad = drad; % RU
    retrieval_results.retrieval_type = variablename;
    retrieval_results.iterations = i;
    retrieval_results.x_retrieved = x(:, 1:i+1); % Final retrieved state (log for q)
    retrieval_results.xa = xa;                 % A priori
    retrieval_results.xtrue = xtrue;           % Truth (from sonde), log for q
    retrieval_results.Sa = Sa;                 % A priori covariance
    retrieval_results.K = K;                   % Jacobian matrix
    retrieval_results.K_t = K_t;               % Jacobian for T
    retrieval_results.K_q = K_q;               % Jacobian for q
    retrieval_results.Se = Se;                 % Measurement error covariance
    retrieval_results.G = G; % Gain matrix
    retrieval_results.Spos = Spos_output(:,:,i);   % Posterior covariance (classic)
    retrieval_results.Spos_CR = Spos_output_CR(:,:,i); % Posterior covariance (Rodgers full form)
    retrieval_results.A = A_output(:,:,i);     % Averaging kernel (classic)
    retrieval_results.A_CR = A_output_CR(:,:,i); % Averaging kernel (Rodgers full form)
    retrieval_results.DFS = DFS_output(1:i);     % Degrees of freedom for signal (classic)
    retrieval_results.DFS_CR = DFS_output_CR(1:i); % DFS (lambda = 0 version)
    retrieval_results.DFS_per_height = DFS_per_height(:,1:i); % DFS per height
    retrieval_results.DFS_per_height_CR = DFS_per_height_CR(:,1:i); % DFS per height (lambda = 0 version)
    retrieval_results.DFS_T = DFS_T_output_CR(1:i); % DFS for T
    retrieval_results.DFS_Q = DFS_Q_output_CR(1:i); % DFS for q
    retrieval_results.J = J(1:i+1);            % Cost function history
    retrieval_results.d = d(1:i);              % State vector convergence history
    retrieval_results.d_threshold = d_threshold; % Threshold for convergence
    retrieval_results.dy = dy(1:i);            % Measurement space convergence history
    retrieval_results.lambda_history = lambda_output(1:i); % LM lambda history
    retrieval_results.F_output = F_output(:, 1:i); % Forward simulations
    retrieval_results.z = z_truth; % km, height levels
    retrieval_results.tfinal =  x(1:38, 1:i+1);
    retrieval_results.qfinal =  x(39:end, 1:i+1);

    if strcmp(variablename, 'both')
        retrieval_results.tx_retrieved = x(1:nlev, i+1);
        retrieval_results.qx_retrieved = exp(x(nlev+1:end, i+1)); % g/kg
    elseif strcmp(variablename, 'T')
        retrieval_results.tx_retrieved = x(:, i+1);
    elseif strcmp(variablename, 'wv')
        retrieval_results.qx_retrieved = exp(x(:, i+1));
    end
    retrieval_results.AERI_wnum = AERI_wnum_adj; % cm^{-1}
    retrieval_results.measurement = measurement; % RU

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    disp('Current directory:');
    disp(pwd);
    disp('Saving retrieval results...');
    disp(['Output directory: ', output_dir]);
    disp(['Output filename: retrieval_results_modtran_', variablename, '_', timestamp, '.mat']);
    output_filename = sprintf('retrieval_results_modtran_%s_%s.mat', variablename, timestamp);    
    save(output_filename, 'retrieval_results', '-v7.3');
    disp(['Saved retrieval results to: ', output_filename]);
    addpath('/home/binmenja/direct/field_campaigns/whafffers/retrievals')
    plot_retrieval_results(output_filename)
end