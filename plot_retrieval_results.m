function [] = plot_retrieval_results(filename)

    load(filename)

    [~, base, ~] = fileparts(filename);
    outdir = fullfile('figures_png', base);
    if ~exist(outdir, 'dir'); mkdir(outdir); end

    fnt = 'Montserrat'; fs = 40;
    nlev = length(retrieval_results.z); % Number of vertical levels

    if ~isfield(retrieval_results, 'Gnew') && isfield(retrieval_results, 'Sa') && ...
            isfield(retrieval_results, 'K') && isfield(retrieval_results, 'Se')
        K = retrieval_results.K;
        Se_inv = inv(retrieval_results.Se);
        Sa_inv = inv(retrieval_results.Sa);
        retrieval_results.Gnew = (K' * Se_inv * K + Sa_inv) \ (K' * Se_inv);
        fprintf('Computed Gnew (Rodgers-style Gain matrix).\n');
    end

    % Helper for saving
    savefigpng = @(name) saveas(gcf, fullfile(outdir, [name, '_',retrieval_results.utc_profile,'.png']));

    % === Temperature Retrieval ===
    if isfield(retrieval_results, 'tx_retrieved')
        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        plot(retrieval_results.tx_retrieved, retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
        plot(retrieval_results.xtrue(1:nlev), retrieval_results.z, 'k--', 'LineWidth', 2);
        plot(retrieval_results.xa(1:nlev), retrieval_results.z, 'b:', 'LineWidth', 2);
        legend('Retrieved T', 'Truth T', 'Prior T', 'FontSize', fs, 'FontName', fnt, 'Location','best');
        xlabel('Temperature [K]'); ylabel('Height [km]');
        title('Temperature Retrieval'); set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('temperature_retrieval');
    end

    % === Water Vapor Retrieval ===
    if isfield(retrieval_results, 'qx_retrieved')
        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        qx = retrieval_results.qx_retrieved;
        q_truth = exp(retrieval_results.xtrue(nlev+1:end));
        q_prior = exp(retrieval_results.xa(nlev+1:end));
        plot(qx, retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
        plot(q_truth, retrieval_results.z, 'k--', 'LineWidth', 2);
        plot(q_prior, retrieval_results.z, 'b:', 'LineWidth', 2);
        legend('Retrieved q', 'Truth q', 'Prior q', 'FontSize', fs, 'FontName', fnt, 'Location','best');
        xlabel('Specific Humidity [ppmv]'); ylabel('Height [km]');
        title('Water Vapor Retrieval'); set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('wv_retrieval');
    end

    % === Temperature Differences ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    plot(retrieval_results.tx_retrieved - retrieval_results.xtrue(1:nlev), retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
    plot(retrieval_results.xa(1:nlev) - retrieval_results.xtrue(1:nlev), retrieval_results.z, 'b--', 'LineWidth', 3);
    legend('Retrieved - Truth', 'Prior - Truth', 'FontSize', fs, 'FontName', fnt);
    xlabel('T [K]'); ylabel('Height [km]');
    title('Temperature Differences'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('temperature_diff');

    % === Water Vapor Differences ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    plot(qx - q_truth, retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
    plot(q_prior - q_truth, retrieval_results.z, 'b--', 'LineWidth', 3);
    legend('Retrieved - Truth', 'Prior - Truth', 'FontSize', fs, 'FontName', fnt);
    xlabel('q [ppmv]'); ylabel('Height [km]');
    title('Water Vapor Differences'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('wv_diff');

    % === Residuals ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    imagesc(retrieval_results.AERI_wnum,1:size(retrieval_results.drad,1),retrieval_results.drad);
    xlabel('Wavenumber [$cm^{-1}$]','Interpreter','latex'); ylabel('Iteration');
    title('Residuals: Meas - Sim [RU]'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('residuals');

    % === Cost Function ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    plot(1:length(retrieval_results.J), retrieval_results.J, '-o', 'LineWidth', 3);
    xlabel('Iteration'); ylabel('Cost Function J');
    title('OEM Cost Function'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('cost_function');

    % === DFS Breakdown ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    plot(1:length(retrieval_results.DFS_CR), retrieval_results.DFS_CR, '-o', 'LineWidth', 3); hold on;
    plot(1:length(retrieval_results.DFS_T), retrieval_results.DFS_T, '--x', 'LineWidth', 3);
    plot(1:length(retrieval_results.DFS_Q), retrieval_results.DFS_Q, ':s', 'LineWidth', 3);
    legend('Total DFS', 'DFS for T', 'DFS for q', 'FontSize', fs, 'FontName', fnt);
    xlabel('Iteration'); ylabel('DFS');
    title('Degrees of Freedom for Signal'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('dfs_breakdown');

    % === DFS per Layer ===
    figure('Color','w', 'Position', [100, 100, 1440, 1440]);
    barh(retrieval_results.DFS_per_height_CR(:,end));
    xlabel('DFS per layer'); ylabel('Layer index');
    title('DFS per Layer (Final Iteration)'); set(gca, 'FontSize', fs, 'FontName', fnt);
    savefigpng('dfs_per_layer');

    % === Averaging Kernel ===
    if isfield(retrieval_results, 'Gnew')
        A = retrieval_results.Gnew * retrieval_results.K;
        AK_diag_T = diag(A(1:nlev,1:nlev));
        AK_diag_q = diag(A(nlev+1:end,nlev+1:end));
        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        plot(AK_diag_T, retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
        plot(AK_diag_q, retrieval_results.z, 'b--', 'LineWidth', 3);
        xlabel('AK Diagonal'); ylabel('Height [km]');
        title('Averaging Kernel Diagonal'); legend('T', 'q');
        set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('ak_diag');
    end

    % === Jacobian Norms ===
    if isfield(retrieval_results, 'K_t') && isfield(retrieval_results, 'K_q')
        Kt_norm = vecnorm(retrieval_results.K_t, 2, 1);
        Kq_norm = vecnorm(retrieval_results.K_q, 2, 1);
        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        plot(Kt_norm, retrieval_results.z, 'r-', 'LineWidth', 3); hold on;
        plot(Kq_norm, retrieval_results.z, 'b--', 'LineWidth', 3);
        xlabel('Jacobian Norm [RU/K or RU/log(q)]'); ylabel('Height [km]');
        legend('|K_t|', '|K_q|'); title('Jacobian Norms');
        set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('jacobian_norms');
    end

    % === Final chi² ===
    if isfield(retrieval_results, 'drad') && isfield(retrieval_results, 'Se')
        r = retrieval_results.drad(end,:)';
        chi2 = r' / retrieval_results.Se * r;
        fprintf('Final iteration chi² = %.2f (expected ~%d)\n', chi2, size(retrieval_results.Se,1));
    end

    % === SVD of K ===
    if isfield(retrieval_results, 'K')
        svals = svd(retrieval_results.K);
        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        semilogy(svals, 'ko-', 'LineWidth', 2);
        xlabel('Index'); ylabel('Singular value');
        title('Singular Values of Jacobian');
        set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('svd');

        figure('Color','w', 'Position', [100, 100, 1440, 1440]);
        cumvar = cumsum(svals.^2) / sum(svals.^2);
        plot(cumvar, 'k-', 'LineWidth', 3); grid on;
        xlabel('Index'); ylabel('Cumulative variance explained');
        title('Cumulative Information Content of K');
        set(gca, 'FontSize', fs, 'FontName', fnt);
        savefigpng('svd_cumulative');
    end
end
