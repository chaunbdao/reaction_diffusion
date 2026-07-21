function periodic_test_optimized(points,runtime,t_corr,plot_each_step)

    % runtime is the total physical simulation time, not the number of
    % integration steps. For example, runtime = 500 simulates t = 0..500.
    % plot_each_step controls live animation. Snapshots are always saved
    % every snapshot_dt, and periodic_test_optimized always performs exactly one run.

    if ~(isscalar(plot_each_step) && ...
            (islogical(plot_each_step) || ...
            (isnumeric(plot_each_step) && ismember(plot_each_step,[0,1]))))
        error('plot_each_step must be a logical scalar (true or false).');
    end
    plot_each_step = logical(plot_each_step);

    %=====================================
    % Important Parameters
    %=====================================

    dt = 0.001;
    snapshot_dt = 0.1;
    corr_time = t_corr/dt;

    num_steps = round(runtime/dt);
    steps_per_snapshot = round(snapshot_dt/dt);
    num_snapshots = round(runtime/snapshot_dt);

    time_tolerance = 100*eps(max([runtime,snapshot_dt,dt]));
    if abs(num_steps*dt-runtime) > time_tolerance
        error('runtime must be an integer multiple of dt.');
    end
    if abs(steps_per_snapshot*dt-snapshot_dt) > time_tolerance
        error('snapshot_dt must be an integer multiple of dt.');
    end
    if abs(num_snapshots*snapshot_dt-runtime) > time_tolerance
        error('runtime must be an integer multiple of snapshot_dt.');
    end

    snapshot_times = (1:num_snapshots)*snapshot_dt;


    % Choose one noise process:
    %   'active_fixed_strength'
    %   'active_fixed_variance'
    %   'gaussian_white'

    noise_mode = ['active_fixed_strength'];

    if (t_corr == 0)
        noise_mode = ['gaussian_white'];
    end

    noise_amplitude = 0.3; % A for fixed-strength active noise
    eta_std_fixed_variance = 2;

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

    e = ones(points,1);
    glaplacian = spdiags([e 0*e e],-1:1,points,points);
    glaplacian(1,points)=1;
    glaplacian(points,1)=1;

    %===================================
    % Precomputing the timestep matrices
    %=====================================

    DX = 1;
    DY = 5;
    %gamma = 5;
    t_v = 10;
    %t_v = 6;
    gamma = 1/t_v;

    betavar = 0.7*gamma;
    alphavar = 0.5*gamma;
    epsilon = 1;
    %epsilon = 0.01;
    a = 0.1*sqrt(epsilon);
    %epsilon = (a*10)^2;

    mu1 = DX*dt/a/a/2;
    mu2 = DY*dt/a/a/2;
    identity = speye(points);

    Amat1 = (1+2*mu1)*identity - mu1*glaplacian;
    Bmat1 = (1-2*mu1)*identity + mu1*glaplacian;
    Amat2 = (1+2*mu2+alphavar*dt/2)*identity - mu2*glaplacian;
    Bmat2 = (1-2*mu2-alphavar*dt/2)*identity + mu2*glaplacian;
    C = (dt/2)*identity;

    BigA = [Amat1 sparse(points,points); -C*gamma Amat2];
    BigA_solver = decomposition(BigA,'lu');

    %=====================================
    % Simulation Loop
    %=====================================

    % Some initial conditions for the single diagnostic run
    X = ones(points,1)*-1.0;
    Y = ones(points,1)*-0.6;
    %Y = ones(points,1)*(0.7 - 1)/0.6;
    %X(1)=1;
    eta = zeros(points,1);

    xvals = linspace(0,2*3.14,points);
    xmat = zeros(points,num_snapshots);
    ymat = zeros(points,num_snapshots);

    if plot_each_step
        figure
        hX = plot(xvals,X);
        hold on
        hY = plot(xvals,Y);
        hold off
        ylim([-2 2])
        drawnow
    end

    snapshot_count = 0;
    for w = 1:num_steps
        [Xnew,Ynew,eta_new] = rd_step_active(X,Y,eta);
        eta=eta_new;
        X = Xnew;
        Y = Ynew;

        if plot_each_step
            hX.YData = X;
            hY.YData = Y;
            drawnow
        end

        if mod(w,steps_per_snapshot) == 0
            snapshot_count = snapshot_count + 1;
            meanX = mean(X);
            meanY = mean(Y);
            stdX = std(X);
            stdY = std(Y);

            [w*dt meanX stdX meanY stdY max(X)]
            xmat(:,snapshot_count) = X;
            ymat(:,snapshot_count) = Y;
        end
    end

    figure
    [tvals,dvals] = meshgrid(snapshot_times,xvals);
    s = surf(tvals,dvals,xmat);
    s.EdgeColor = 'none';
    xlabel('time')
    ylabel('space')
    writematrix(xmat,"periodic_test_xmat.csv")
    writematrix(ymat,"periodic_test_ymat.csv")

    xmask = xmat>0.5;
    CC = bwconncomp(xmask,8);
    CC2 = CC2periodic(CC,[1,0]);

    comparison_stats = periodic_component_stats(CC2, points, 1);
    comparison_stats = comparison_stats( ...
        comparison_stats.SpaceTimeArea >= 100 & ...
        ~comparison_stats.TouchesInitialTime & ...
        ~comparison_stats.TouchesFinalTime,:);
    periodic_test_summary = periodic_summary_table(comparison_stats)
    writetable(comparison_stats, "periodic_test_data.csv")
    writetable(periodic_test_summary, "periodic_test_summary.csv")

    function summary_table = periodic_summary_table(stats_table)
        summary_table = stats_table(:,{ ...
            'Component', ...
            'PeriodicHeightPixels', ...
            'RegionPropsHeightPixels', ...
            'HeightDifference', ...
            'DurationPixels', ...
            'DurationTime', ...
            'RegionPropsDurationPixels', ...
            'RegionPropsDurationTime', ...
            'DurationDifference', ...
            'SpaceTimeArea', ...
            'RegionPropsArea', ...
            'AreaDifference', ...
            'TouchesInitialTime', ...
            'TouchesFinalTime'});
    end

    function stats_table = periodic_component_stats(CC2, points, run_count)
        region_stats = regionprops("table", CC2, "BoundingBox", "Area");
        num_objects = CC2.NumObjects;

        Run = zeros(num_objects,1);
        Component = zeros(num_objects,1);
        PeriodicHeightPixels = zeros(num_objects,1);
        DurationPixels = zeros(num_objects,1);
        DurationTime = zeros(num_objects,1);
        SpaceTimeArea = zeros(num_objects,1);
        RegionPropsHeightPixels = zeros(num_objects,1);
        RegionPropsDurationPixels = zeros(num_objects,1);
        RegionPropsDurationTime = zeros(num_objects,1);
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
            DurationTime(k) = DurationPixels(k)*snapshot_dt;
            SpaceTimeArea(k) = numel(rows);

            RegionPropsXMin(k) = region_stats.BoundingBox(k,1);
            RegionPropsYMin(k) = region_stats.BoundingBox(k,2);
            RegionPropsDurationPixels(k) = region_stats.BoundingBox(k,3);
            RegionPropsDurationTime(k) = RegionPropsDurationPixels(k)*snapshot_dt;
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
            DurationTime, ...
            SpaceTimeArea, ...
            RegionPropsHeightPixels, ...
            RegionPropsDurationPixels, ...
            RegionPropsDurationTime, ...
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

    function [Xnew, Ynew, eta_new] = rd_step_active(X,Y,eta)

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

    %function [patterns] = analyze(X)

    %end

end
