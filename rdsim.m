function rdsim(points,runtime,nruns,tcorr_start,tcorr_end,tcorr_points,tv_start,tv_end,tv_points,workers)
    
    %===================
    % Initializing workers for parallelization
    %===================

    % create a local cluster object
    pc = parcluster('local');

    % explicitly set JobStorageLocation to job-specific temp directory
    mkdir(strcat('/scratch/',getenv('USER'),'/slurmJobs/',getenv('slurm_ID')));
    pc.JobStorageLocation = strcat('/scratch/', getenv('USER'),'/slurmJobs/', getenv('slurm_ID'));

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
        "_step_", filename_number(tcorr_step), ".csv");
    pdmat_filename = strcat("pdmat_", filename_parameters);
    prpmat_filename = strcat("prpmat_", filename_parameters);
    pltmat_filename = strcat("pltmat_", filename_parameters);

    for tv_count = 1:tv_points
        for tcorr_count = 1:tcorr_points
            tv = tv_start + tv_step * (tv_count - 1);
            tcorr = tcorr_start + tcorr_step * (tcorr_count - 1);
            [pdmat(tv_count,tcorr_count), ...
                prpmat(tv_count,tcorr_count), ...
                pltmat(tv_count,tcorr_count)] = gen_pd_point_par(points,runtime,nruns,tcorr,tv);
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
