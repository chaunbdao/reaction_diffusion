function results = compare_gen_pd_point(points,runtime,runs,t_corr,t_v,workers)
%COMPARE_GEN_PD_POINT Time serial and parallel gen_pd_point implementations.
%   Run with no inputs to compare 10 runs using default simulation settings:
%       results = compare_gen_pd_point()
%
%   Or override parameters:
%       results = compare_gen_pd_point(1000,5000,10,5,50,4)

    if nargin < 1 || isempty(points)
        points = 1000;
    end
    if nargin < 2 || isempty(runtime)
        runtime = 5000;
    end
    if nargin < 3 || isempty(runs)
        runs = 10;
    end
    if nargin < 4 || isempty(t_corr)
        t_corr = 0.5;
    end
    if nargin < 5 || isempty(t_v)
        t_v = 25;
    end
    if nargin < 6 || isempty(workers)
        workers = 4;
    end

    fprintf("Benchmarking with points=%d, runtime=%d, runs=%d, t_corr=%g, t_v=%g\n", ...
        points, runtime, runs, t_corr, t_v);

    rng(1);
    serial_timer = tic;
    serial_value = gen_pd_point(points,runtime,runs,t_corr,t_v);
    serial_time = toc(serial_timer);

    has_parallel_toolbox = license("test","Distrib_Computing_Toolbox") && exist("gcp","file") == 2;
    if ~has_parallel_toolbox
        parallel_value = NaN;
        parallel_time = NaN;
        speedup = NaN;

        results = table( ...
            serial_value, ...
            parallel_value, ...
            serial_time, ...
            parallel_time, ...
            speedup, ...
            workers);

        fprintf("Serial:   value=%g, time=%.2f s\n", serial_value, serial_time);
        fprintf("Parallel: unavailable; Parallel Computing Toolbox is required for gcp/parpool/parfor.\n");
        return
    end

    pool = gcp("nocreate");
    if isempty(pool) || pool.NumWorkers ~= workers
        if ~isempty(pool)
            delete(pool);
        end
        parpool("local", workers);
    end

    rng(1);
    parallel_timer = tic;
    parallel_value = gen_pd_point_par(points,runtime,runs,t_corr,t_v);
    parallel_time = toc(parallel_timer);

    speedup = serial_time / parallel_time;

    results = table( ...
        serial_value, ...
        parallel_value, ...
        serial_time, ...
        parallel_time, ...
        speedup, ...
        workers);

    fprintf("Serial:   value=%g, time=%.2f s\n", serial_value, serial_time);
    fprintf("Parallel: value=%g, time=%.2f s, workers=%d\n", parallel_value, parallel_time, workers);
    fprintf("Speedup:  %.2fx\n", speedup);
end
