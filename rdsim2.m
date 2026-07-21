function rdsim2(points,physical_runtime,nruns,tcorr_start,tcorr_end,tcorr_points,tv_start,tv_end,tv_points,workers,dt,A)

    % physical_runtime is simulation time, not the number of integration steps.
    % Dynamics are integrated with dt and saved for patch analysis every 0.1.
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
    if ~(isnumeric(A) && isreal(A) && isscalar(A) && ...
            ~isnan(A) && ~isinf(A) && A >= 0)
        error('A must be a nonnegative finite scalar.');
    end

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

    pdmat = -1*ones(tv_points, tcorr_points);
    prpmat = -1*ones(tv_points, tcorr_points);
    pltmat = -1*ones(tv_points, tcorr_points);
    tv_step = (tv_end - tv_start)/(tv_points-1);
    tcorr_step = (tcorr_end - tcorr_start)/(tcorr_points-1);
    filename_parameters = strcat("tv_", filename_number(tv_start), ...
        "_to_", filename_number(tv_end), ...
        "_step_", filename_number(tv_step), ...
        "_tcorr_", filename_number(tcorr_start), ...
        "_to_", filename_number(tcorr_end), ...
        "_step_", filename_number(tcorr_step), ...
        "_dt_", filename_number(dt), ...
        "_runtime_", filename_number(physical_runtime), ...
        "_A_", filename_number(A), ...
        "_snapshotdt_", filename_number(snapshot_dt), ".csv");
    pdmat_filename = strcat("pdmat2_", filename_parameters);
    prpmat_filename = strcat("prpmat2_", filename_parameters);
    pltmat_filename = strcat("pltmat2_", filename_parameters);

    for tv_count = 1:tv_points
        for tcorr_count = 1:tcorr_points
            tv = tv_start + tv_step * (tv_count - 1);
            tcorr = tcorr_start + tcorr_step * (tcorr_count - 1);
            [pdmat(tv_count,tcorr_count), ...
                prpmat(tv_count,tcorr_count), ...
                pltmat(tv_count,tcorr_count)] = gen_pd_point_par2( ...
                    points,physical_runtime,nruns,tcorr,tv,dt,A,snapshot_dt);
            writematrix(pdmat, pdmat_filename);
            writematrix(prpmat, prpmat_filename);
            writematrix(pltmat, pltmat_filename);
        end
    end
end

function text = filename_number(value)
    text = sprintf('%.12g', value);
    text = strrep(text, '-', 'm');
    text = strrep(text, '.', 'p');
    text = strrep(text, '+', '');
end
