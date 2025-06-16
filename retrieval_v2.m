close all; clear all; clc;
% Created by Benjamin Riot--Bretecher (2025). Based code from Lei Liu.
addpath('/home/binmenja/direct/matlab/mylib')
addpath('/home/binmenja/direct/models/modtran6/bin/linux')
utc_profile = '202502120753'; % UTC time of the radiosonde launch
if ~exist(utc_profile,'dir')
    system(strcat("mkdir -p ", utc_profile));
end

path_modtran = '/home/binmenja/direct/models/modtran6/bin/linux/';
v1 = 520;
v2 = 1800;
resolution = 0.1; % cm-1
fwhm = resolution*2;
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

have_jacobian_ready = 0; % whether the jacobian is already ready
have_jacobian_iready = 1; %whether jacobian ready for first iteration 
fprintf('Have jacobian ready: %d\n', have_jacobian_ready);
fprintf('Have jacobian for first iteration: %d\n', have_jacobian_iready);
lambda_1st = 10000;
is_q_log = 0; % whether q is already in log scale; 0: no, 1: yes
fprintf('Is q in log scale: %d\n', is_q_log);

dx_wv = 0.1; % 10 precent uncertainty for wv
dx_t = 0.5; % 0.5K uncertainty for T
cf_id = 0; % 1: combine dx and cf convergence criteria, 0: only use dx convergence critieria
cf_grad = 1e-3; % abs(J(i+1) - J(i)) < cf_grad * J(i)

aeri_file = '/home/binmenja/direct/field_campaigns/whafffers/aeri_data/gault/atmospheric.physics_aeri_125_MCGILL_gault_20250212T000424Z_20250212T235249Z_001.nc';

% Convolute total error to AERI resolution
AERI_MOPD = ncread(aeri_file,'MOPD'); % AERI MOPD
AERI_fwhm = 1.2067/(2*AERI_MOPD); % AERI FWHM
AERI_cal_error_full = ncread(aeri_file, 'RADIANCE.SKY_ERROR');
AERI_nesr_full = ncread(aeri_file, 'RADIANCE.SKY_NOISE'); % [RU]
AERI_rad_full = ncread(aeri_file, 'RADIANCE.SKY_CLEANED'); % [RU]
AERI_wnum = ncread(aeri_file, 'WAVENUMBER'); %[cm-1]
AERI_dates_full = ncread(aeri_file, 'DATETIME'); % Format is YYYYMMDDThhmmssZ
% Convert the char var AERI_dates to string
AERI_dates_full = string(AERI_dates_full.');

% Custom approach for ISO 8601 format
years = double(extractBetween(AERI_dates_full, 1, 4));
months = double(extractBetween(AERI_dates_full, 5, 6));
days = double(extractBetween(AERI_dates_full, 7, 8));
% disp(days)
% Find position of 'T' to locate the time part
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
% Find indices for -2 to +8 minutes around the launch time
time_diff_duration = AERI_dates_full - ref_date_launch; % Get duration object
time_diff = minutes(time_diff_duration); % Convert duration to numeric minutes
start_idx = find(time_diff >= -2, 1, 'first'); % First index where time_diff >= -2 minutes
end_idx = find(time_diff <= 8, 1, 'last');    % Last index where time_diff <= +8 minutes
% Ensure indices are within valid bounds
if isempty(start_idx), start_idx = 1; end
if isempty(end_idx), end_idx = length(AERI_dates_full); end

% Extract the timestamp around the radiosonde launch time
AERI_rad = AERI_rad_full(:, start_idx:end_idx);
AERI_cal_error = AERI_cal_error_full(:, start_idx:end_idx); % 1-sigma absolute error
AERI_nesr = AERI_nesr_full(:, start_idx:end_idx); % 1-sigma noise eq. spe. rad.

% retrieval type: 1(T only), 2(wv only), 3(T and wv)
retrieval_type = 3; %1/2/3 
switch retrieval_type
case 1
    variablename = 'T';
case 2
    variablename = 'wv';
case 3
    variablename = 'both';
end

% read priori
data_priori = importdata('/home/binmenja/direct/field_campaigns/whafffers/retrievals/era5_priori_data.mat'); % ERA5 data, not yet available
z_truth = (data_priori.z(:,1) + 12)./1000; % m, sfc = 115m + 12m = 127m
z_prior = (data_priori.z(:,1) + 12)./1000;
t_prior = data_priori.t_all; % K
q_prior = data_priori.q_all; % g/kg

x0_t = nanmean(t_prior,2);
x0_q = nanmean(q_prior,2); 

nlev = length(z_truth);

% read truth
load('/home/binmenja/direct/field_campaigns/whafffers/profile_202502120753UTC.mat'); % Already patched with ERA5 data
z_sonde = profile.z(:); % km
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

if ~is_q_log % enters if q is not log scale
    q_prior = log(q_prior);
    q_truth = log(q_truth);
    x0_q = log(x0_q);
    qStr = 'log(q)'; tStr = 'T';
    
else
    disp('q already in log scale');
    qStr = 'q'; tStr = 'T';
end

% Sa from ERA5
if strcmp(variablename, 'T')
    Sa = nancov(t_prior'); % covariance matrix of T
elseif  strcmp(variablename, 'wv')
    Sa = nancov(q_prior'); % covariance matrix of q
elseif  strcmp(variablename, 'both')
    disp('look here:')
    combined_prior = [t_prior', q_prior']; % combine T and q
    size(combined_prior)
    Sa = nancov(combined_prior); % covariance matrix of T and q
    size(Sa)
    disp('---')
end

% increase Sa
Sa = Sa;

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

% Se: only diag part

% aeri measurement and Se prescription
aeri_rad_mean = mean(AERI_rad,2,'omitnan'); % RU
nesr_mean = mean(AERI_nesr,2,'omitnan');
% nesr_mean(nesr_mean == 0) = eps; 
Se = diag(nesr_mean.^2);
% Se

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


%% retrieval
lambda = lambda_1st; %for Levenberg-Marquardt regularization
lambda_output = NaN(1, 20); 
T = zeros(length(Sa),length(Se));
for i = 1:20
    tic
     disp(['ITERATION: ', num2str(i)]);
    
    if strcmp(variablename, 'T')
        qx = exp(q_truth);
        tx = t_truth;
        tx = x(:,i);
    elseif strcmp(variablename, 'wv')
        qx = exp(q_truth);
        qx = exp(x(:,i));
        tx = t_truth;   
    elseif strcmp(variablename, 'both')
        qx = exp(q_truth);
        qx = exp(x(nlev+1:end,i));
        tx = t_truth;
        tx = x(1:nlev,i);
    end

   

    if have_jacobian_ready
        if i==1
            load(strcat('./',utc_profile,'/K_t_era5_x0.mat')) % Unit is W/(cm^2*sr*cm^{-1})/K
            K_t = jacobian_info.jacobian .* 1e7; % convert the unit to RU/K
            load(strcat('./',utc_profile,'/K_q_era5_x0.mat'))% Unit is W/(cm^2*sr*cm^{-1})/log(g/kg)
            K_q = jacobian_info.jacobian .* 1e7; % convert to RU/log(g/kg)
        elseif i>=2
            load(strcat('./',utc_profile,'/K_t_era5_x1.mat')) % Unit is W/(cm^2*sr*cm^{-1})/K
            K_t = jacobian_info.jacobian .* 1e7; % convert the unit to RU/K
            load(strcat('./',utc_profile,'/K_q_era5_x1.mat'))% Unit is W/(cm^2*sr*cm^{-1})/log(g/kg)
            K_q = jacobian_info.jacobian .* 1e7; % convert to RU/log(g/kg)
        end
    else 
        if i==1
            profile.z = z_truth; % km
            profile.p = p; % hPa
            profile.t = tx; % K
            profile.q = qx; % g/kg
            profile.co2 = co2; % ppmv
            profile.o3 = o3; % g/kg
            profile.ch4 = ch4; % ppmv
            profile.co = co; % ppmv
            if ~have_jacobian_iready
                compute_modtran_jacobian_temperature(utc_profile,profile,strcat('./',utc_profile,'/K_t_era5_x0.mat'));
            end
            load(strcat('./',utc_profile,'/K_t_era5_x0.mat')) % Unit is W/(cm^2*sr*cm^{-1})/K
            K_t = jacobian_info.jacobian .* 1e7; % convert the unit to RU/K
            if ~have_jacobian_iready
                compute_modtran_jacobian_q(utc_profile,profile,strcat('./',utc_profile,'/K_q_era5_x0.mat'));
            end
            load(strcat('./',utc_profile,'/K_q_era5_x0.mat'))% Unit is W/(cm^2*sr*cm^{-1})/log(g/kg)
            K_q = jacobian_info.jacobian .* 1e7; % convert to RU/log(g/kg)
            sim_wnum = jacobian_info.wavenumbers;
            K_t(isnan(K_t)) = min(abs(K_t(:)));
            K_q(isnan(K_q)) = min(abs(K_q(:)));
        elseif i==2
            profile.q = exp(x(nlev+1:end, i)); % updated q
            profile.t = x(1:nlev, i); % updated T
            compute_modtran_jacobian_temperature(utc_profile,profile,strcat('./',utc_profile,'/K_t_era5_x1.mat'));
            load(strcat('./',utc_profile,'/K_t_era5_x1.mat')) % Unit is W/(cm^2*sr*cm^{-1})/K
            K_t = jacobian_info.jacobian .* 1e7; % convert the unit to RU/K
            compute_modtran_jacobian_q(utc_profile,profile,strcat('./',utc_profile,'/K_q_era5_x1.mat'));
            load(strcat('./',utc_profile,'/K_q_era5_x1.mat'))% Unit is W/(cm^2*sr*cm^{-1})/log(g/kg)
            K_q = jacobian_info.jacobian .* 1e7; % convert to RU/log(g/kg)
            sim_wnum = jacobian_info.wavenumbers;
	        K_t(isnan(K_t)) = min(abs(K_t(:)));
            K_q(isnan(K_q)) = min(abs(K_q(:)));
        else
            load(strcat('./',utc_profile,'/K_t_era5_x1.mat')) % Unit is W/(cm^2*sr*cm^{-1})/K
            K_t = jacobian_info.jacobian .* 1e7; % convert the unit to RU/K
            load(strcat('./',utc_profile,'/K_q_era5_x1.mat'))% Unit is W/(cm^2*sr*cm^{-1})/log(g/kg)
            K_q = jacobian_info.jacobian .* 1e7; % convert to RU/log(g/kg)
            sim_wnum = jacobian_info.wavenumbers;
	        K_t(isnan(K_t)) = min(abs(K_t(:)));
            K_q(isnan(K_q)) = min(abs(K_q(:)));
        end
    end
    cloud.qi = []; cloud.ql = []; cloud.z = [];
    profile.t = tx;
    profile.q = qx;
    ts = 0; % space temp
    %F = run_single_simulation(path_modtran, path_tape5, ...
    %    profile, modroot, resolution, AERI_fwhm, cloud, v1, v2, ...
    %    angle, iLoc, emis, profile.z(end));
    tic
    F = run_single_simulation(path_modtran, path_tape5, ...
            profile, modroot, resolution, fwhm, cloud, ...
            v1, v2, angle, iLoc, emis, profile.z(end), ts);
    F = F.rad_plt .* 1e7; % convert the unit to RU
    toc
    fprintf('Forward simulation done.\n');

    % Adjust to AERI resolution
    for il=1:size(K_t,2)
       K_t_adj(:,il) = band_conv_brb(sim_wnum, K_t(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
       K_q_adj(:,il) = band_conv_brb(sim_wnum, K_q(:,il), AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
    end
    F = band_conv_brb(sim_wnum, F, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
    size(F)
    K_t = K_t_adj;
    K_q = K_q_adj;
    

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
        disp(['Iteration ', num2str(i)]);
        disp(['Size of K: ', mat2str(size(K))]);
        disp(['Size of Se: ', mat2str(size(Se))]);
        disp(['Size of Sa: ', mat2str(size(Sa))]);
        disp(['Size of measurement: ', mat2str(size(measurement))]);
        disp(['Size of F: ', mat2str(size(F))]);
        disp(['Size of x(:,i): ', mat2str(size(x(:,i)))]);
        disp(['Size of xa: ', mat2str(size(xa))]);
        if any(isnan(K(:))) || any(isnan(Se(:))) || any(isnan(Sa(:))) || any(isnan(measurement)) || any(isnan(F)) || any(isnan(x(:,i))) || any(isnan(xa))
            error('NaN values detected in inputs at iteration %d', i);
        end
        if rcond(Se) < eps
            error('Se is singular at iteration %d. RCOND = %e', i, rcond(Se));
        end
        H = (1+lambda)*pinv(Sa) + K'*pinv(Se)*K;
        if rcond(H) < eps
            error('Hessian matrix is singular at iteration %d. RCOND = %e', i, rcond(H));
        end
        x(:, i+1) = gather(x(:,i) + inv((1+lambda)*inv(Sa) + K'*inv(Se)*K)*(K'*inv(Se)*(measurement - F)-inv(Sa)*(x(:,i)-xa)));


        if strcmp(variablename, 'T')
            tx = x(:,i+1);
        elseif strcmp(variablename, 'wv')
            qx = exp(x(:,i+1));
        elseif strcmp(variablename, 'both')
            qx = exp(x(nlev+1:end,i+1));
            tx = x(1:nlev,i+1);
        end
        if any(tx <= 0)
            disp('Non-positive temperatures found:');
            disp(find(tx <= 0));
            J(i+1) = NaN;
        else
            cloud.qi = []; cloud.ql = []; cloud.z = [];
            profile.t = tx;
            profile.q = qx;
	    ts = 0 ;
            %F_new = run_single_simulation(path_modtran, path_tape5, profile, modroot, resolution, fwhm, cloud, v1, v2, ...
            %angle, iLoc, emis, profile.z(end));
            tic
            F_new = run_single_simulation(path_modtran, path_tape5, ...
            profile, modroot, resolution, fwhm, cloud, ...
            v1, v2, angle, iLoc, emis, profile.z(end), ts);
            F_new = F_new.rad_plt .* 1e7;
            toc
            fprintf('Forward simulation done for iteration %d.\n', i+1);
            F_new = band_conv_brb(sim_wnum, F_new, AERI_wnum_adj, AERI_fwhm, AERI_MOPD, 'Sinc');
            J(i+1) = (measurement - F_new)'*inv(Se)*(measurement - F_new) + (x(:,i+1) - xa)'*inv(Sa)*(x(:,i+1) - xa);
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
    if isnan(lambda_output(i))
        
        fprintf('Warning: lambda_output(%d) was not assigned (cost function did not decrease). I must check the values of the matrices.\n', i);
        
    end
    Spos = inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K)*((lambda_output(i)+1).^2*inv(Sa)+K'*inv(Se)*K)*inv((lambda_output(i)+1)*inv(Sa)+K'*inv(Se)*K);
    Spos_output(:,:,i) = gather(Spos);

    Mnew = pinv(K'*inv(Se)*K+inv(Sa)+lambda_output(i)*inv(Sa));
    Gnew = Mnew * K' * inv(Se);
    T = Gnew + ((eye(length(Sa))-Gnew*K-Mnew*inv(Sa)))*T;

    A_output_CR(:,:,i) = T*K;
    DFS_output_CR(i) = trace(A_output_CR(:,:,i));
    Spos_output_CR(:,:,i) = T*Se*T'; 

    %dx2_threshold
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

    A = gather(inv(K'*inv(Se)*K+(lambda_output(i)+1).*inv(Sa)))*gather(K'*inv(Se)*K); 
    DFS = trace(A); 

    A_output(:,:,i) = A; 
    DFS_output(i) = DFS;
    DFS_per_height(:,i) = diag(A); % BRB
    d(i) = (x(:,i)-x(:,i+1))' * (pinv(Sa)+K'*pinv(Se)*K) * (x(:,i)-x(:,i+1));
    
    Sy = Se * inv(K*Sa*K'+Se) * Se;

    dy(i) = (F-F_new)'*inv(Sy)*(F-F_new);

    fprintf('Lambda: %.3e\n', lambda_output(i));
    fprintf('Cost function J(i): %.4f\n', J(i));
    fprintf('Step size d(i): %.4e\n', d(i));
    fprintf('Threshold d_threshold: %.4e\n', d_threshold);
    fprintf('DFS (standard): %.2f\n', DFS_output(i));
    fprintf('DFS (constrained retrieval): %.2f\n', DFS_output_CR(i));
    fprintf('Degrees of freedom per height: %.2f\n', DFS_per_height(:,i));
    fprintf('Measurement mismatch dy(i): %.4e\n', dy(i));

    if  d(i) < min(d_threshold,length(Sa)./20)
        break
    end
    toc
end

retrieval_results = struct();

retrieval_results.utc_profile = utc_profile;
retrieval_results.retrieval_type = variablename;
retrieval_results.iterations = i;
retrieval_results.x_retrieved = x(:, 1:i+1); % Final retrieved state
retrieval_results.xa = xa;                 % A priori
retrieval_results.xtrue = xtrue;           % Truth (from sonde)
retrieval_results.Sa = Sa;                 % A priori covariance
retrieval_results.Spos = Spos_output(:,:,i);   % Posterior covariance (classic)
retrieval_results.Spos_CR = Spos_output_CR(:,:,i); % Posterior covariance (Rodgers full form)
retrieval_results.A = A_output(:,:,i);     % Averaging kernel (classic)
retrieval_results.A_CR = A_output_CR(:,:,i); % Averaging kernel (Rodgers full form)
retrieval_results.DFS = DFS_output(1:i);     % Degrees of freedom for signal (classic)
retrieval_results.DFS_CR = DFS_output_CR(1:i); % DFS (Rodgers)
retrieval_results.DFS_per_height = DFS_per_height(:,1:i); % New: Save DFS per height
retrieval_results.J = J(1:i+1);            % Cost function history
retrieval_results.d = d(1:i);              % State vector convergence history
retrieval_results.dy = dy(1:i);            % Measurement space convergence history
retrieval_results.lambda_history = lambda_output(1:i); % LM lambda history
retrieval_results.F_output = F_output(:, 1:i); % Forward simulations

if strcmp(variablename, 'both')
    retrieval_results.tx_retrieved = x(1:nlev, i+1);
    retrieval_results.qx_retrieved = exp(x(nlev+1:end, i+1));
elseif strcmp(variablename, 'T')
    retrieval_results.tx_retrieved = x(:, i+1);
elseif strcmp(variablename, 'wv')
    retrieval_results.qx_retrieved = exp(x(:, i+1));
end

timestamp = datestr(now,'yyyymmdd_HHMMSS');
output_filename = sprintf('./%s/retrieval_results_%s_%s.mat', utc_profile, variablename, timestamp);
save(output_filename, 'retrieval_results', '-v7.3');
disp(['Saved retrieval results to: ', output_filename]);
