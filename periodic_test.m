function periodic_test(points,runtime,t_corr,runs)

    %=====================================
    % Important Parameters
    %=====================================

    dt = 0.1;
    corr_time = t_corr/dt;


    %Turn on only one of the noise processes below, choosed fixed
    %integrated strength of fixed instantaneous variance

    %=====
    % Noise configutation for fixed integrated strength
    % C(Delta t) = (A/tauc) exp(-|Delta t|/tauc)

    noise_amplitude = 1.25;
    sigma_active = sqrt(2*noise_amplitude)/t_corr;
    %=====

    %=====
    % Noise configutation for fixed instantaneous variance
    % C(Delta t) = eta_std^2 exp(-|Delta t|/tauc)

    %eta_std = 0.5;
    %sigma_active = eta_std*sqrt(2/t_corr);
    %=====


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

    csv = '.csv';
    runtxt = 'runs';
    pointtxt = 'points';
    tau='tau';
    time='time';
    simdat = 'simdat1d_';
    filename = strcat(simdat,num2str(points),pointtxt,num2str(runtime),time,num2str(corr_time),tau,num2str(runs),runtxt,csv);
    datamat=zeros(runtime,runs);
    for run_count = 1:runs
        % Some initial conditions for system
        X = ones(points,1)*-1.0;
        Y = ones(points,1)*-0.6;
        %X(1)=1;
        eta = zeros(points,1);

        if (runs==1)
            figure
            xvals = linspace(0,2*3.14,points);
            xmat = [];
            ymat = [];
        end

        for w = 1:runtime
            [Xnew,Ynew,eta_new] = rd_step_active(X,Y,eta);
            eta=eta_new;

            meanX = mean(Xnew);
            meanY = mean(Ynew);
            stdX = std(Xnew);
            stdY = std(Ynew);

            [run_count w meanX stdX meanY stdY max(Xnew)]
            X = Xnew;
            Y = Ynew;
            datamat(w,run_count)=meanX;
            if (runs==1)
                plot(xvals,X)
                hold on
                plot(xvals,Y)
                hold off
                ylim([-2 2])
                pause(0.05)
                xmat = [xmat X];
                ymat = [ymat Y];
            end
        end
    end
    if (runs>1)
        writematrix(datamat,filename)
    else
        [tvals,dvals] = meshgrid(1:runtime,xvals);
        s = surf(tvals,dvals,xmat);
        s.EdgeColor = 'none';
        writematrix(xmat,"periodic_test_xmat.csv")
        writematrix(ymat,"periodic_test_ymat.csv")

        xmat = xmat>0.5;
        CC = bwconncomp(xmat,8);
        CC2 = CC2periodic(CC,[1,0]);

        comparison_stats = periodic_component_stats(CC2, points, run_count);
        comparison_stats = comparison_stats( ...
            comparison_stats.SpaceTimeArea >= 100 & ...
            ~comparison_stats.TouchesInitialTime & ...
            ~comparison_stats.TouchesFinalTime,:);
        periodic_test_summary = periodic_summary_table(comparison_stats)
        writetable(comparison_stats, "periodic_test_data.csv")
        writetable(periodic_test_summary, "periodic_test_summary.csv")
    end

    function summary_table = periodic_summary_table(stats_table)
        summary_table = stats_table(:,{ ...
            'Component', ...
            'PeriodicHeightPixels', ...
            'RegionPropsHeightPixels', ...
            'HeightDifference', ...
            'DurationPixels', ...
            'RegionPropsDurationPixels', ...
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
        DX =1;
        DY = 5;
        %gamma = 5;
        t_v = 100;
        %t_v = 6;
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
        eta_new = (1-1/corr_time)*eta + normrnd(0,sigma_active*sqrt(dt),points,1);
        Xnew = Xnew + eta_new*dt;
        %Ynew = Ynew + normrnd(0,0.5*sqrt(dt),points,1);
    end

    %function [patterns] = analyze(X)

    %end

end
