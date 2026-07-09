function generate_data_1d(points,runtime,t_corr,runs)

    %=====================================
    % Important Parameters 
    %=====================================

    dt = 0.1;
    corr_time = t_corr/dt;


    %Turn on only one of the noise processes below, choosed fixed
    %integrated strength of fixed instantaneous variance

    %=====
    % Noise configutation for fixed integrated strength 
    % C(Δt) = (A/τc) exp(-|Δt|/τc)
    
    %noise_amplitude = 1.25;
    %sigma_active = sqrt(2*noise_amplitude)/t_corr;
    %=====

    %=====
    % Noise configutation for fixed instantaneous variance
    % C(Δt) = eta_std² exp(-|Δt|/τc)
    
    eta_std = 0.5;
    sigma_active = eta_std*sqrt(2/t_corr);
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
        writematrix(xmat,"xmat.csv")
        writematrix(ymat,"ymat.csv")

        xmat = xmat>0.5;
        CC = bwconncomp(xmat,8);
        CC2 = CC2periodic(CC,[1,0]);
        stats = regionprops("table", CC2, "BoundingBox", "Area");
        if isempty(stats) || isempty(stats.BoundingBox)
            newstats = stats
        else
            idx = stats.Area>100 & stats.BoundingBox(:,4)<1000;
            newstats = stats(idx,:)
        end
    end

    function [Xnew, Ynew, eta_new] = rd_step_active(X,Y,eta)
		
		%Constants relevant to the equation
        %==================================
		DX =1;
		DY = 5;
        %gamma = 5;
        t_v = 25;
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
