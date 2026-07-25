function rdsim3(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,A_start,A_end,A_points,workers,dt,tv)

    % Sweep fixed-strength active-noise amplitude A and correlation time
    % tau_c at a fixed tau_v. At tau_c = 0, gen_pd_point_par2 uses the
    % corresponding Gaussian white-noise limit. Rows of each output matrix
    % correspond to A; columns correspond to tau_c.
    snapshot_dt = 0.1;

    if ~(isnumeric(dt) && isreal(dt) && isscalar(dt) && ...
            ~isnan(dt) && ~isinf(dt) && dt > 0)
        error('dt must be a positive finite scalar.');
    end
    if ~(isnumeric(physical_runtime) && isreal(physical_runtime) && ...
            isscalar(physical_runtime) && ~isnan(physical_runtime) && ...
            ~isinf(physical_runtime) && physical_runtime > 0)
        error('physical_runtime must be a positive finite scalar.');
    end
    if ~(isnumeric(tv) && isreal(tv) && isscalar(tv) && ...
            ~isnan(tv) && ~isinf(tv) && tv > 0)
        error('tv must be a positive finite scalar.');
    end
    if ~(isnumeric(tcorr_start) && isreal(tcorr_start) && ...
            isscalar(tcorr_start) && isfinite(tcorr_start) && tcorr_start >= 0 && ...
            isnumeric(tcorr_end) && isreal(tcorr_end) && ...
            isscalar(tcorr_end) && isfinite(tcorr_end) && ...
            tcorr_end >= tcorr_start)
        error(['tcorr_start and tcorr_end must be nonnegative finite scalars, ' ...
            'with tcorr_end >= tcorr_start.']);
    end
    if ~(isnumeric(A_start) && isreal(A_start) && isscalar(A_start) && ...
            isfinite(A_start) && A_start >= 0 && ...
            isnumeric(A_end) && isreal(A_end) && isscalar(A_end) && ...
            isfinite(A_end) && A_end >= A_start)
        error(['A_start and A_end must be nonnegative finite scalars, ' ...
            'with A_end >= A_start.']);
    end
    validate_sweep_points(tcorr_points, 'tcorr_points');
    validate_sweep_points(A_points, 'A_points');

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

    % create a local cluster object
    pc = parcluster('local');

    % explicitly set JobStorageLocation to job-specific temp directory
    mkdir(strcat('/scratch/',getenv('USER'),'/slurmJobs/',getenv('slurm_ID')));
    pc.JobStorageLocation = strcat('/scratch/', getenv('USER'),'/slurmJobs/',getenv('slurm_ID'));

    pool = gcp("nocreate");
    if isempty(pool)
      parpool(pc,workers);
    end

    pdmat = -1*ones(A_points, tcorr_points);
    prpmat = -1*ones(A_points, tcorr_points);
    pltmat = -1*ones(A_points, tcorr_points);
    A_values = linspace(A_start, A_end, A_points);
    tcorr_values = linspace(tcorr_start, tcorr_end, tcorr_points);
    A_step = sweep_step(A_values);
    tcorr_step = sweep_step(tcorr_values);
    filename_parameters = strcat("A_", filename_number(A_start), ...
        "_to_", filename_number(A_end), ...
        "_step_", filename_number(A_step), ...
        "_tcorr_", filename_number(tcorr_start), ...
        "_to_", filename_number(tcorr_end), ...
        "_step_", filename_number(tcorr_step), ...
        "_tv_", filename_number(tv), ...
        "_dt_", filename_number(dt), ...
        "_runtime_", filename_number(physical_runtime), ...
        "_snapshotdt_", filename_number(snapshot_dt), ".csv");
    pdmat_filename = strcat("pdmat3_", filename_parameters);
    prpmat_filename = strcat("prpmat3_", filename_parameters);
    pltmat_filename = strcat("pltmat3_", filename_parameters);

    for A_count = 1:A_points
        for tcorr_count = 1:tcorr_points
            A = A_values(A_count);
            tcorr = tcorr_values(tcorr_count);
            [pdmat(A_count,tcorr_count), ...
                prpmat(A_count,tcorr_count), ...
                pltmat(A_count,tcorr_count)] = gen_pd_point_par2( ...
                    points,physical_runtime,nruns,tcorr,tv,dt,A,snapshot_dt);
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

function step = sweep_step(values)
    step = values(2) - values(1);
end

function text = filename_number(value)
    text = sprintf('%.12g', value);
    text = strrep(text, '-', 'm');
    text = strrep(text, '.', 'p');
    text = strrep(text, '+', '');
end
