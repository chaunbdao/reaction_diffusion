function [pd_val, patch_run_percent, avg_patch_lifetime] = gen_pd_point_par2(points,physical_runtime,runs,t_corr,t_v,dt,A,snapshot_dt)

    %=====================================
    % Important Parameters
    %=====================================

    corr_time = t_corr/dt;
    num_steps = round(physical_runtime/dt);
    steps_per_snapshot = round(snapshot_dt/dt);
    num_snapshots = round(physical_runtime/snapshot_dt);


    % Choose one noise process:
    %   'active_fixed_strength'
    %   'active_fixed_variance'
    %   'gaussian_white'

    noise_mode = ['active_fixed_strength'];

    if (t_corr == 0)
        noise_mode = ['gaussian_white'];
    end

    noise_amplitude = A; % A for fixed-strength active noise
    eta_std_fixed_variance = 2;
    sigma_active = NaN;
    sigma_white = NaN;

    if strcmp(noise_mode,'active_fixed_strength')
        % C(Delta t) = (A/tauc) exp(-|Delta t|/tauc)
        sigma_active = sqrt(2*noise_amplitude)/t_corr;
    elseif strcmp(noise_mode,'active_fixed_variance')
        % C(Delta t) = eta_std^2 exp(-|Delta t|/tauc)
        eta_std = eta_std_fixed_variance;
        sigma_active = eta_std*sqrt(2/t_corr);
    elseif strcmp(noise_mode,'gaussian_white')
        % White-noise limit of C(Delta t) = (A/tauc) exp(-|Delta t|/tauc)
        sigma_white = sqrt(2*noise_amplitude);
    else
        error('Unknown noise_mode: %s', noise_mode);
    end


    %===================================
    % Creating the Graph Laplacian
    %=====================================
    % might consider putting more effort into optimizing this

    e = ones(points,1);
    glaplacian = spdiags([e 0*e e],-1:1,points,points);
    glaplacian(1,points)=1;
    glaplacian(points,1)=1;

    %===================================
    % Precomputing the timestep matrices
    %===================================

    DX = 1;
    DY = 5;
    gamma = 1/t_v;

    betavar = 0.7*gamma;
    alphavar = 0.5*gamma;
    epsilon = 1;
    a = 0.1*sqrt(epsilon);

    mu1 = DX*dt/a/a/2;
    mu2 = DY*dt/a/a/2;
    identity = speye(points);

    Bmat1 = (1-2*mu1)*identity + mu1*glaplacian;
    Bmat2 = (1-2*mu2-alphavar*dt/2)*identity + mu2*glaplacian;
    Amat1 = (1+2*mu1)*identity - mu1*glaplacian;
    Amat2 = (1+2*mu2+alphavar*dt/2)*identity - mu2*glaplacian;
    C = (dt/2)*identity;

    BigA = [Amat1 sparse(points,points); -C*gamma Amat2];

    % decomposition objects cannot be serialized reliably to process
    % workers. Build one from the precomputed matrix on each worker, then
    % reuse it for every timestep and every run assigned to that worker.
    BigA_solver_constant = parallel.pool.Constant( ...
        @() decomposition(BigA,'lu'));

    %=====================================
    % Simulation Loop
    %=====================================

    total_patch_spatial_size = 0;
    total_patch_lifetime = 0;
    num_patches = 0;
    runs_with_patches = 0;
    parfor run_count = 1:runs
        BigA_solver = BigA_solver_constant.Value;

        % Some initial conditions for system
        xmat = zeros(points,num_snapshots);
        X = ones(points,1)*-1.0;
        Y = ones(points,1)*-0.4;
        %Y = ones(points,1)*(0.7 - 1)/0.6;
        %X(1)=1;
        eta = zeros(points,1);
        snapshot_count = 0;
        for w = 1:num_steps
            [Xnew,Ynew,eta_new] = rd_step_active( ...
                X,Y,eta,dt,points,noise_mode,corr_time,sigma_active,sigma_white, ...
                gamma,betavar,epsilon,Bmat1,Bmat2,BigA_solver);
            eta=eta_new;

            X = Xnew;
            Y = Ynew;
            if mod(w,steps_per_snapshot) == 0
                snapshot_count = snapshot_count + 1;
                xmat(:,snapshot_count) = X;
            end
        end

        xmask = xmat>0.5;
        CC = bwconncomp(xmask,8);
        CC2 = CC2periodic(CC,[1,0]);

        stats = periodic_component_stats(CC2, points, run_count);

        stats = stats( ...
            stats.SpaceTimeArea >= 100 & ...
            ~stats.TouchesInitialTime & ...
            ~stats.TouchesFinalTime,:);

        patches_this_run = height(stats);
        total_patch_spatial_size = total_patch_spatial_size + sum(stats.PeriodicHeightPixels);
        total_patch_lifetime = total_patch_lifetime + sum(stats.DurationPixels) * snapshot_dt;
        num_patches = num_patches + patches_this_run;
        if patches_this_run > 0
            runs_with_patches = runs_with_patches + 1;
        end
    end

    if num_patches > 0
        pd_val = total_patch_spatial_size / num_patches;
        avg_patch_lifetime = total_patch_lifetime / num_patches;
    else
        pd_val = 0;
        avg_patch_lifetime = 0;
    end

    patch_run_percent = 100 * runs_with_patches / runs;

end

function stats_table = periodic_component_stats(CC2, points, run_count)
        region_stats = regionprops("table", CC2, "BoundingBox", "Area");
        num_objects = CC2.NumObjects;

        Run = zeros(num_objects,1);
        Component = zeros(num_objects,1);
        PeriodicHeightPixels = zeros(num_objects,1);
        DurationPixels = zeros(num_objects,1);
        SpaceTimeArea = zeros(num_objects,1);
        RegionPropsHeightPixels = zeros(num_objects,1);
        RegionPropsDurationPixels = zeros(num_objects,1);
        RegionPropsArea = zeros(num_objects,1);
        RegionPropsXMin = zeros(num_objects,1);
        RegionPropsYMin = zeros(num_objects,1);
        HeightDifference = zeros(num_objects,1);
        DurationDifference = zeros(num_objects,1);
        AreaDifference = zeros(num_objects,1);
        TouchesInitialTime = false(num_objects,1);
        TouchesFinalTime = false(num_objects,1);

        for k = 1:num_objects
            [rows,cols] = ind2sub(CC2.ImageSize, CC2.PixelIdxList{k});

            Run(k) = run_count;
            Component(k) = k;
            PeriodicHeightPixels(k) = circular_height(unique(rows), points);
            DurationPixels(k) = max(cols) - min(cols) + 1;
            SpaceTimeArea(k) = numel(rows);

            RegionPropsXMin(k) = region_stats.BoundingBox(k,1);
            RegionPropsYMin(k) = region_stats.BoundingBox(k,2);
            RegionPropsDurationPixels(k) = region_stats.BoundingBox(k,3);
            RegionPropsHeightPixels(k) = region_stats.BoundingBox(k,4);
            RegionPropsArea(k) = region_stats.Area(k);
            HeightDifference(k) = RegionPropsHeightPixels(k) - PeriodicHeightPixels(k);
            DurationDifference(k) = RegionPropsDurationPixels(k) - DurationPixels(k);
            AreaDifference(k) = RegionPropsArea(k) - SpaceTimeArea(k);
            TouchesInitialTime(k) = any(cols == 1);
            TouchesFinalTime(k) = any(cols == CC2.ImageSize(2));
        end

        stats_table = table( ...
            Run, ...
            Component, ...
            PeriodicHeightPixels, ...
            DurationPixels, ...
            SpaceTimeArea, ...
            RegionPropsHeightPixels, ...
            RegionPropsDurationPixels, ...
            RegionPropsArea, ...
            RegionPropsXMin, ...
            RegionPropsYMin, ...
            HeightDifference, ...
            DurationDifference, ...
            AreaDifference, ...
            TouchesInitialTime, ...
            TouchesFinalTime);
end

function height_pixels = circular_height(rows, points)
        rows = sort(rows(:));

        if isempty(rows)
            height_pixels = 0;
        elseif numel(rows) == 1
            height_pixels = 1;
        else
            gaps = diff(rows);
            wrap_gap = rows(1) + points - rows(end);
            height_pixels = points - max([gaps; wrap_gap]) + 1;
        end
end

function [Xnew, Ynew, eta_new] = rd_step_active( ...
    X,Y,eta,dt,points,noise_mode,corr_time,sigma_active,sigma_white, ...
    gamma,betavar,epsilon,Bmat1,Bmat2,BigA_solver)

        phi = (1/epsilon)*(1.-(X.*X)).*(X-Y);
        phinew = phi;

        b1 = Bmat1*X+dt/2*(phi+phinew);
        b2 = dt*betavar + Bmat2*Y + gamma*X*dt/2;
        XYnew = BigA_solver\[b1;b2];
        Xnew = XYnew(1:points);
        Ynew = XYnew(points+1:end);
        phinewer = (1/epsilon)*(1.-(Xnew.*Xnew)).*(Xnew-Ynew);

        %fixed point iteration, up to machine precision
        iterations = 0;
        while (norm(phinewer-phinew)>(1*10^-10) && iterations <100)
            b1 = Bmat1*X+dt/2*(phi+phinewer);
            XYnew = BigA_solver\[b1;b2];
            Xnew = XYnew(1:points);
            Ynew = XYnew(points+1:end);
            phinew = phinewer;
            phinewer = (1/epsilon)*(1.-(Xnew.*Xnew)).*(Xnew-Ynew);
            iterations = iterations + 1; 
        end
        if strcmp(noise_mode,'gaussian_white')
            eta_new = normrnd(0,sigma_white/sqrt(dt),points,1);
            Xnew = Xnew + eta_new*dt;
        elseif strcmp(noise_mode,'active_fixed_strength') || strcmp(noise_mode,'active_fixed_variance')
            eta_new = (1-1/corr_time)*eta + normrnd(0,sigma_active*sqrt(dt),points,1);
            Xnew = Xnew + eta_new*dt;
        else
            error('Unknown noise_mode: %s', noise_mode);
        end
        %Ynew = Ynew + normrnd(0,0.5*sqrt(dt),points,1);
end
