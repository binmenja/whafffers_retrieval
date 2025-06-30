close all; clear all; clc;
% Created by Benjamin Riot--Bretecher (2025). 
addpath('/home/binmenja/direct/matlab/mylib')
addpath('/home/binmenja/direct/models/modtran6/bin/linux')
addpath('/home/binmenja/direct/models/modtran6/bin/linux/matlib_modtran')
addpath('/home/binmenja/direct/models/LBLRTM/scripts_brb')

utc_profile = '202502120753'; % UTC time of the radiosonde launch
if ~exist(utc_profile,'dir')
    system(strcat("mkdir -p ", utc_profile));
end

% Retrieval settings
have_jacobian_ready = 0; % whether the jacobian is already ready
fprintf('Have jacobian ready: %d\n', have_jacobian_ready);
lambda_1st = 10000;
is_q_log = 0; % whether q is already in log scale; 0: no, 1: yes
dx_wv = 0.1; % 10 percent uncertainty for wv
dx_t = 0.5; % 0.5K uncertainty for T
altitude_toa=100;
cloud.qi = []; cloud.ql = []; cloud.z = []; % clear-sky


fprintf('Is q in log scale: %d\n', is_q_log);
path_modtran = '/home/binmenja/direct/models/modtran6/bin/linux/'; v1 = 0.2;
v2 = 2008; resolution = 0.1; fwhm = resolution*2; zenith = 0; atmos_id = '5';
path_tape5 = ['/home/binmenja/direct/field_campaigns/whafffers/tape5/', utc_profile, '/'];
if ~exist(path_tape5,'dir')
    system(["mkdir -p ", path_tape5]);
end
modroot = ['whafffers_', utc_profile, '_retrieval'];
angle = 180; emis = 0; iLoc = '0'; ts = 0;

radiosonde_file = '/home/binmenja/direct/field_campaigns/whafffers/radiosonde_data/balloon_sonde_mcgill_gault_20250212T075314Z_20250212T090502Z_001.nc';
sonde_dt = ncread(radiosonde_file, 'DATETIME'); % Format days since 2000-01-01 00:00:00 (MJD2K)
ref_date_launch = datetime('2000-01-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC') + sonde_dt(1); % Convert to datetime
disp(ref_date_launch)

% AERI
aeri_file = '/home/binmenja/direct/field_campaigns/whafffers/aeri_data/gault/atmospheric.physics_aeri_125_MCGILL_gault_20250212T000424Z_20250212T235249Z_001.nc';
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
start_idx = find(time_diff >= -1, 1, 'first'); % start (~ 1 minute before launch)
end_idx = find(time_diff <= 10, 1, 'last');    % end (10 minutes after launch)
if isempty(start_idx), start_idx = 1; end
if isempty(end_idx), end_idx = length(AERI_dates_full); end
AERI_rad = AERI_rad_full(:, start_idx:end_idx);
AERI_cal_error = AERI_cal_error_full(:, start_idx:end_idx); % 1-sigma absolute error
AERI_nesr = AERI_nesr_full(:, start_idx:end_idx); % 1-sigma noise eq. spe. rad.

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
q_prior = convert_gkg_to_ppmv(q_prior, 'H2O'); % convert to ppmv
x0_t = nanmean(t_prior,2);
x0_q = nanmean(q_prior,2); 
nlev = length(z_truth);

% read truth
load('/home/binmenja/direct/field_campaigns/whafffers/profile_202502120753UTC.mat'); % Already patched with ERA5 data
z_sonde = profile.z(:) + z_truth(1); % km; adjust to same base height!
p_sonde = profile.p(:); % hPa
t_sonde = profile.t(:); % K
q_sonde = profile.q(:); % g/kg
q_sonde = convert_gkg_to_ppmv(q_sonde, 'H2O'); % convert to ppmv
co2_sonde = profile.co2(:); % ppmv
ch4_sonde = profile.ch4(:); % ppmv
o3_sonde = profile.o3; % g/kg
o3_sonde = convert_gkg_to_ppmv(o3_sonde, 'O3'); % convert to ppmv

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
profile_input.o3 = o3; % ppmv
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

    if ~have_jacobian_ready
        if i==1
            disp('Computing jacobian for the first iteration...');
            tic
            [wv_aj,K_t,K_q] = lblrtm_aeri_aj_clearsky(v1,v2,profile_input.z,profile_input.p,profile_input.t,profile_input.t(1),profile_input.q,...
                profile_input.co2,profile_input.o3,n2o,co,profile_input.ch4,ccl4,f11,f12,zenith,atmos_id,altitude_toa); % W/cm^2/sr/cm^-1/K or log(mixing ratio)
            K_t = K_t .* 1e7; % convert the unit to RU/K
            K_q = K_q .* 1e7; % convert the unit to RU/log(mixing ratio (ppmv))
            K_t(isnan(K_t)) = min(abs(K_t(:)));  
            K_q(isnan(K_q)) = min(abs(K_q(:)));  
            save K_t_era5_x0.mat K_t;
            save K_q_era5_x0.mat K_q;
            save wv_aj_era5_x0.mat wv_aj;
            toc
        elseif i==2
            tic
            disp('entering the second iteration, recomputing jacobian');
            [wv_aj,K_t,K_q] = lblrtm_aeri_aj_clearsky(v1,v2,profile_input.z,profile_input.p,profile_input.t,profile_input.t(1),profile_input.q,...
                profile_input.co2,profile_input.o3,n2o,co,profile_input.ch4,ccl4,f11,f12,zenith,atmos_id,altitude_toa); % W/cm^2/sr/cm^-1/K or log(mixing ratio)
            K_t = K_t .* 1e7; % convert the unit to RU/K
            K_q = K_q .* 1e7; % convert the unit to RU/log(mixing ratio (ppmv))
            K_t(isnan(K_t)) = min(abs(K_t(:)));  
            K_q(isnan(K_q)) = min(abs(K_q(:)));  
            save K_t_era5_x1.mat K_t;
            save K_q_era5_x1.mat K_q;
            save wv_aj_era5_x1.mat wv_aj;
            toc
        else
            disp('Using jacobian from the second iteration. No need to recompute or reload.');
        end
        
    else 
        if i==1
            wv_aj = importdata('wv_aj_era5_x0.mat');
            K_t = importdata('K_t_era5_x0.mat'); % RU/K
            K_q = importdata('K_q_era5_x0.mat'); % RU/log(mixing ratio (ppmv))
        elseif i>=2
            wv_aj = importdata('wv_aj_era5_x1.mat');
            K_t = importdata('K_t_era5_x1.mat'); % RU/K
            K_q = importdata('K_q_era5_x1.mat'); % RU/log(mixing ratio (ppmv))
        else
            disp('Using jacobian from the second iteration. No need to recompute or reload.');
        end
    end

    tic
    if i==1
        disp('Computing forward simulation for the first iteration...');
        tic
        [~, F] = lblrtm_aeri_clearsky(v1,v2,profile_input.z,profile_input.p,profile_input.t,profile_input.t(1),profile_input.q,profile_input.co2,...
            profile_input.o3,n2o,co,profile_input.ch4,ccl4,f11,f12,zenith,atmos_id,altitude_toa);
        toc
    else
        F = F_new; % use the previous forward simulation result, which is the same
    end
    toc
    [~,index_aj_begin] = find(abs(wv_aj - wv_begin) == min(abs(wv_aj - wv_begin)));
    [~,index_aj_end] = find(abs(wv_aj - wv_end) == min(abs(wv_aj - wv_end)));

    if( (index_aj_end-index_aj_begin+1)~=(index_aeri_end-index_aeri_begin+1) ...
        | abs(nanmean(wv_aj(index_aj_begin:index_aj_end) - AERI_wnum_adj(index_aeri_begin:index_aeri_end)')) > 0.01 )
        disp('aj wv not consistent');
        return;
    end


    if strcmp(variablename, 'both')
        K = [K_t(index_aj_begin:index_aj_end,:) K_q(index_aj_begin:index_aj_end,:)]; % dimension: (nwv,nlev x 2)
    elseif strcmp(variablename,'T')
        K = K_t(index_aj_begin:index_aj_end,:); % dimension: (nwv,nlev)
    else 
        K = K_q(index_aj_begin:index_aj_end,:); % dimension: (nwv,nlev)
    end

    F = (F(index_aj_begin:index_aj_end))' .* 1e7; % convert radiance unit to RU

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
            disp('Computing forward simulation for the next iteration...');
            tic
            [~, F_new] = lblrtm_aeri_clearsky(v1,v2,z_truth,p,tx,tx(1),qx,co2,o3,n2o,co,ch4,ccl4,f11,f12,zenith,atmos_id,altitude_toa);
            toc
            F_new = (F_new(index_aj_begin:index_aj_end))' .* 1e7;
            
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

retrieval_results.utc_profile = utc_profile;
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
output_filename = sprintf('./%s/retrieval_results_lblrtm_%s_%s.mat', utc_profile, variablename, timestamp);
save(output_filename, 'retrieval_results', '-v7.3');
disp(['Saved retrieval results to: ', output_filename]);
