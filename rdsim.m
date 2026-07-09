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
      parpool(pc,workers)
    end

    pdmat = -1*ones(tv_points, tcorr_points);
    tv_step = (tv_end - tv_start)/(tv_points-1);
    tcorr_step = (tcorr_end - tcorr_start)/(tcorr_points-1);
    for tv_count = 1:tv_points
        for tcorr_count = 1:tcorr_points
            tv = tv_start + tv_step * (tv_count - 1);
            tcorr = tcorr_start + tcorr_step * (tcorr_count - 1);
            pdmat(tv_count,tcorr_count) = gen_pd_point_par(points,runtime,nruns,tcorr,tv);
        end 
    end 
    writematrix(pdmat, "pdmat.csv")
end