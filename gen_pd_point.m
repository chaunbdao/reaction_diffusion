function pd_val = gen_pd_point(points,runtime,runs,t_corr,t_v)

    %=====================================
    % Important Parameters
    %=====================================

    dt = 0.1;
    corr_time = t_corr/dt;


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
    % might consider putting more effort into optimizing this

    e = ones(points,1);
    glaplacian = spdiags([e 0*e e],-1:1,points,points);
    glaplacian(1,points)=1;
    glaplacian(points,1)=1;

    %=====================================
    % Simulation Loop
    %=====================================

    total_patch_spatial_size = 0;
    num_patches = 0;
    for run_count = 1:runs
        % Some initial conditions for system
        xmat = [];
        X = ones(points,1)*-1.0;
        Y = ones(points,1)*-0.4;
        %Y = ones(points,1)*(0.7 - 1)/0.6;
        %X(1)=1;
        eta = zeros(points,1);
        for w = 1:runtime
            [Xnew,Ynew,eta_new] = rd_step_active(X,Y,eta);
            eta=eta_new;

            X = Xnew;
            Y = Ynew;
            xmat = [xmat X];
            [run_count w]
        end

        xmask = xmat>0.5;
        CC = bwconncomp(xmask,8);
        CC2 = CC2periodic(CC,[1,0]);

        stats = periodic_component_stats(CC2, points, run_count);

        stats = stats( ...
            stats.SpaceTimeArea >= 100 & ...
            ~stats.TouchesInitialTime & ...
            ~stats.TouchesFinalTime,:);

        total_patch_spatial_size = total_patch_spatial_size + sum(stats.PeriodicHeightPixels);
        num_patches = num_patches + height(stats);
    end

    if num_patches > 0
        pd_val = total_patch_spatial_size / num_patches;
    else
        pd_val = 0;
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

function [Xnew, Ynew, eta_new] = rd_step_active(X,Y,eta)

        %Constants relevant to the equation
        %==================================
        DX = 1;
        DY = 5;
        %gamma = 5;
        gamma = 1/t_v;

        betavar = 0.7*gamma;
        alphavar = 0.5*gamma;
        epsilon = 1;
        %epsilon = 0.01;
        a = 0.1*sqrt(epsilon);
        %epsilon = (a*10)^2;
        %===================================


        phi = (1/epsilon)*(1.-(X.*X)).*(X-Y);
        phinew = phi;

        mu1 = DX*dt/a/a/2;
        mu2 = DY*dt/a/a/2;
        Amat1 = -1*mu1*glaplacian;
        Amat2 = -1*mu2*glaplacian;
        Bmat1 = mu1*glaplacian;
        Bmat2 = mu2*glaplacian;

        for m = 1:points
            Amat1(m,m) = (1+2*mu1);
            Bmat1(m,m) = (1-2*mu1);
            Amat2(m,m) = (1+2*mu2) + alphavar*dt/2;
            Bmat2(m,m) = (1-2*mu2) - alphavar*dt/2;
        end

        C = spdiags((dt/2)*ones(points,1),0,points,points);

        BigA = [Amat1 zeros(points); -C*gamma Amat2];
        b1 = Bmat1*X+dt/2*(phi+phinew);
        b2 = dt*betavar + Bmat2*Y + gamma*X*dt/2;
        XYnew = BigA\[b1;b2];
        Xnew = XYnew(1:points);
        Ynew = XYnew(points+1:end);
        phinewer = (1/epsilon)*(1.-(Xnew.*Xnew)).*(Xnew-Ynew);

        %fixed point iteration, up to machine precision
        while norm(phinewer-phinew)>(1*10^-15)
            b1 = Bmat1*X+dt/2*(phi+phinew);
            XYnew = BigA\[b1;b2];
            Xnew = XYnew(1:points);
            Ynew = XYnew(points+1:end);
            phinew = phinewer;
            phinewer = (1/epsilon)*(1.-(Xnew.*Xnew)).*(Xnew-Ynew);
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
