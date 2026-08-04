function rdsim4(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,tv_start,tv_end,tv_points,workers,dt,sigma_eta)

    % Sweep tau_v and tau_c using fixed-instantaneous-variance active noise.
    % Rows of each output matrix correspond to tau_v; columns correspond to
    % tau_c. This workflow uses Y(t=0) = -6.
    snapshot_dt = 0.1;
    noise_mode = "active_fixed_variance";
    Y_initial = -6;

    if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && ...
            isfinite(dt) && dt > 0)
        error('dt must be a positive finite scalar.');
    end
    if ~(isnumeric(physical_runtime) && isreal(physical_runtime) && ...
            isscalar(physical_runtime) && isfinite(physical_runtime) && ...
            physical_runtime > 0)
        error('physical_runtime must be a positive finite scalar.');
    end
    if ~(isnumeric(sigma_eta) && isreal(sigma_eta) && ...
            isscalar(sigma_eta) && isfinite(sigma_eta) && sigma_eta >= 0)
        error('sigma_eta must be a nonnegative finite scalar.');
    end
    if ~(isnumeric(tcorr_start) && isreal(tcorr_start) && ...
            isscalar(tcorr_start) && isfinite(tcorr_start) && tcorr_start > 0 && ...
            isnumeric(tcorr_end) && isreal(tcorr_end) && ...
            isscalar(tcorr_end) && isfinite(tcorr_end) && ...
            tcorr_end >= tcorr_start)
        error(['tcorr_start and tcorr_end must be positive finite scalars, ' ...
            'with tcorr_end >= tcorr_start.']);
    end
    if ~(isnumeric(tv_start) && isreal(tv_start) && ...
            isscalar(tv_start) && isfinite(tv_start) && tv_start > 0 && ...
            isnumeric(tv_end) && isreal(tv_end) && isscalar(tv_end) && ...
            isfinite(tv_end) && tv_end >= tv_start)
        error(['tv_start and tv_end must be positive finite scalars, ' ...
            'with tv_end >= tv_start.']);
    end
    validate_sweep_points(tcorr_points, 'tcorr_points');
    validate_sweep_points(tv_points, 'tv_points');

    num_steps = round(physical_runtime/dt);
    steps_per_snapshot = round(snapshot_dt/dt);
    num_snapshots = round(physical_runtime/snapshot_dt);
    time_tolerance = 100*eps(max([physical_runtime,snapshot_dt,dt]));

    if abs(num_steps*dt-physical_runtime) > time_tolerance
        error('physical_runtime must be an integer multiple of dt.');
    end
    if abs(steps_per_snapshot*dt-snapshot_dt) > time_tolerance
        error('snapshot_dt = 0.1 must be an integer multiple of dt.');
    end
    if abs(num_snapshots*snapshot_dt-physical_runtime) > time_tolerance
        error('physical_runtime must be an integer multiple of snapshot_dt = 0.1.');
    end

    %===================
    % Initializing workers for parallelization
    %===================

    pc = parcluster('local');
    mkdir(strcat('/scratch/',getenv('USER'),'/slurmJobs/',getenv('slurm_ID')));
    pc.JobStorageLocation = strcat('/scratch/', getenv('USER'), ...
        '/slurmJobs/',getenv('slurm_ID'));

    pool = gcp("nocreate");
    if isempty(pool)
        parpool(pc,workers);
    end

    pdmat = -1*ones(tv_points, tcorr_points);
    prpmat = -1*ones(tv_points, tcorr_points);
    pltmat = -1*ones(tv_points, tcorr_points);
    tv_values = linspace(tv_start, tv_end, tv_points);
    tcorr_values = linspace(tcorr_start, tcorr_end, tcorr_points);
    tv_step = tv_values(2) - tv_values(1);
    tcorr_step = tcorr_values(2) - tcorr_values(1);

    filename_parameters = strcat("noise_", noise_mode, ...
        "_sigmaeta_", filename_number(sigma_eta), ...
        "_Yinitial_", filename_number(Y_initial), ...
        "_tv_", filename_number(tv_start), ...
        "_to_", filename_number(tv_end), ...
        "_step_", filename_number(tv_step), ...
        "_tcorr_", filename_number(tcorr_start), ...
        "_to_", filename_number(tcorr_end), ...
        "_step_", filename_number(tcorr_step), ...
        "_dt_", filename_number(dt), ...
        "_runtime_", filename_number(physical_runtime), ...
        "_snapshotdt_", filename_number(snapshot_dt), ".csv");
    pdmat_filename = strcat("pdmat4_", filename_parameters);
    prpmat_filename = strcat("prpmat4_", filename_parameters);
    pltmat_filename = strcat("pltmat4_", filename_parameters);

    for tv_count = 1:tv_points
        for tcorr_count = 1:tcorr_points
            tv = tv_values(tv_count);
            tcorr = tcorr_values(tcorr_count);
            [pdmat(tv_count,tcorr_count), ...
                prpmat(tv_count,tcorr_count), ...
                pltmat(tv_count,tcorr_count)] = gen_pd_point_par2( ...
                    points,physical_runtime,nruns,tcorr,tv,dt,sigma_eta, ...
                    snapshot_dt,noise_mode,Y_initial);
            writematrix(pdmat, pdmat_filename);
            writematrix(prpmat, prpmat_filename);
            writematrix(pltmat, pltmat_filename);
        end
    end
end

function validate_sweep_points(value, argument_name)
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && value >= 2 && value == floor(value))
        error('%s must be an integer greater than or equal to 2.', argument_name);
    end
end

function text = filename_number(value)
    text = sprintf('%.12g', value);
    text = strrep(text, '-', 'm');
    text = strrep(text, '.', 'p');
    text = strrep(text, '+', '');
end
