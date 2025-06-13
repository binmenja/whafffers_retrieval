function [] = compute_modtran_jacobian_temperature(date_str, profile, saveFileName)
    % Perturbs T at each level, computes (ΔRadiance / ΔT), and saves results

    perturbation_k = 1;  % K
    if nargin < 3
        error('Must provide date_str, profile, and saveFileName');
    end

    addpath('/home/binmenja/direct/models/modtran6/bin/linux/matlib_modtran');
    path_modtran = '/home/binmenja/direct/models/modtran6/bin/linux/';
    output_base = ['/home/binmenja/direct/field_campaigns/whafffers/retrievals/', date_str, '/'];
    path_tape5 = ['/home/binmenja/direct/field_campaigns/whafffers/tape5/', date_str, '/'];
    system(['mkdir -p ', path_tape5]);

    resolution = 0.1; fwhm = resolution * 2;
    v1 = 520; v2 = 1800;
    angle = 180; emis = 0; iLoc = '0';
    zt = profile.z(end); ts = 0;
    cloud.qi = []; cloud.ql = []; cloud.z = [];

    % Check if baseline exists
    baseline_file = fullfile(output_base, 'baseline_spectra_aeri.mat');
    %if exist(baseline_file, 'file')
    %    load(baseline_file, 'baseline_spectra');
    %disp('Loaded baseline spectra.');
    %else
        disp('Running baseline simulation...');
        baseline_modroot = ['baseline_', date_str];
        baseline_spectra = run_single_simulation(path_modtran, path_tape5, ...
            profile, baseline_modroot, resolution, fwhm, cloud, ...
            v1, v2, angle, iLoc, emis, zt, ts);
        save(baseline_file, 'baseline_spectra', '-v7.3');
        try
            cleanup_modtran_files(path_modtran, baseline_modroot);
        catch
            warning('Could not clean up baseline MODTRAN files.');
        end
    %end

    n_levels = length(profile.z);
    jacobian_matrix = zeros(length(baseline_spectra.rad_plt), n_levels);
    perturbed_spectra = cell(n_levels, 1);
    actual_perturbation = zeros(n_levels, 1);

    for level = 1:n_levels
        perturbed_profile = profile;
        perturbed_profile.t(level) = perturbed_profile.t(level) + perturbation_k;

        level_modroot = sprintf('jac_temppert_%s_level_%d', date_str, level);
        disp(['Running simulation for level ', num2str(level)]);

        perturbed_result = run_single_simulation(path_modtran, path_tape5, ...
            perturbed_profile, level_modroot, resolution, fwhm, cloud, ...
            v1, v2, angle, iLoc, emis, zt, ts);

        delta_rad = perturbed_result.rad_plt - baseline_spectra.rad_plt;
        jacobian_matrix(:, level) = delta_rad / perturbation_k;

        perturbed_spectra{level} = perturbed_result;
        actual_perturbation(level) = perturbation_k;

        try
            current_path = pwd;
            current_path = strcat(current_path, '/');
            cleanup_modtran_files(path_modtran, level_modroot);
        catch
            warning(['Could not delete files in linux for ', level_modroot]);
            try
                cleanup_modtran_files(current_path, level_modroot);
            catch
                warning(['Could not delete files for ', level_modroot, ' in current directory']);
            end
        end
    end

    % Save everything
    jacobian_info.jacobian = jacobian_matrix;
    jacobian_info.baseline_spectra = baseline_spectra;
    jacobian_info.perturbed_spectra = perturbed_spectra;
    jacobian_info.date = date_str;
    jacobian_info.wavenumbers = baseline_spectra.wavnum_plt;
    jacobian_info.z_levels = profile.z;
    jacobian_info.method = 'temperature';
    jacobian_info.perturbation_k = perturbation_k;
    jacobian_info.actual_perturbations = actual_perturbation;

    save(saveFileName, 'jacobian_info', '-v7.3');
    disp(['Saved temperature Jacobian to: ', saveFileName]);
end
