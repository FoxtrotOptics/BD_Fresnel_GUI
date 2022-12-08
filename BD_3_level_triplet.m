%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   BD_Fresnel_GUI.m                
%
%   Initial version (1.0):    10/01/15 - Manuel Ferdinandus
%   Latest revision (1.6):    03/09/17 - Manuel Ferdinandus
%
%   GUI wrapper for the beam deflection signal calculator
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT

function output = BD_3_level_triplet(alpha_2_e,I_0_e,lambda_e,n_g_e,tau_e,w_0_e,...
    n_g_p,L,N0,X0,...
    sigma_S_01,sigma_S_12,sigma_T_12,...
    tau_ISC,tau_S_10,tau_S_21,tau_TS,tau_T_21,...
    t_in,tau_d,x,y,diagnostics)

    if diagnostics
    
        tic
        disp([datestr(now),': BD_3_level_triplet - Calculating excited state effects.'])

    end
        
    %% physical constants

    h = 6.626070e-34; % (J/s);
    c = 299792458; % (m/s);

    %% calculated parameters

    v_e = c/n_g_e; % group veloctity of excitation (m/s)
    v_p = c/n_g_p; % group velocity of probe (m/s)

    E_photon = (h*c)/lambda_e; % excitation photon energy (J)

    %% propagation parameters
    
    pop_names = {'N_S_0', 'N_S_1', 'N_S_2', 'N_T_1', 'N_T_2'};

    indata = zeros(1,length(pop_names));
    indata(1) = N0; % the initial conditions for this position

    z_steps = 50; % steps through the sample
    
    irr_pts = 10; % radial points
    
    % ODE tolerance parameters
    relTol = 1e-3;
    absTol = 1e-3;
    absTol_arr = absTol*ones(length(indata),1); % absolute error tolerarance array

    %% propagate the excitation

    % define array of z values
    z = linspace(0,L,z_steps);
    delta_z = z(2) - z(1); % step size in space

    % define the time range
    t_min = min(tau_d);
    t_max = max(tau_d);
    t = [min([-3*tau_e,t_min]),t_max]; % range in time

    rad_fac = linspace(1,0,irr_pts); % radial irradiance factor

    delay_points = length(tau_d); % delay points

    I_p = ones(delay_points,z_steps,irr_pts);

    options = odeset('RelTol',relTol,'AbsTol',absTol_arr); % options for the ODE

    I_0_e_norm = I_0_e*repmat(rad_fac,[length(z),1]); % array of initial irradiance for transmission calculation
    I_0_e_arr = I_0_e_norm; % array of stored irradiances
    
    for l = 1:irr_pts

        for i = 1:length(z) - 1 % loop in space through the sample

            % the differential equation to solve
            diff_eq = @(t, y) (excite_prop(t,y,z(i),...
            I_0_e_arr(i,l),v_e,v_p,tau_e,E_photon,...
            sigma_S_01,sigma_S_12,sigma_T_12,...
            tau_ISC,tau_S_10,tau_S_21,tau_TS,tau_T_21,...
            alpha_2_e));

            [t,Y] = ode15s(diff_eq,t,indata,options); % run the solver to step in time

            if i*l == 1 % set up the output array on the first iteration

                Y_out = zeros(length(t),length(z),length(indata),irr_pts); % define the output array
                [~,zero_ind] =  min(abs(t)); % find the zero position element

            end

            Y_out(:,i,:,l) = Y; % save the populations

            % update the irradiance for this step
            I_0_e_arr(i + 1,l) = I_0_e_arr(i,l) -(sigma_S_01*Y(zero_ind,1) + sigma_S_12*Y(zero_ind,2) + sigma_T_12*Y(zero_ind,4) +...
                alpha_2_e*I_0_e_arr(i,l))*I_0_e_arr(i,l)*delta_z;

        end

    end

    %interpolate the time to the delay axis
    [T_grid,z_grid,N_grid,r_grid] = ndgrid(t,z,1:length(indata),1:irr_pts); % z, t, 2-d grid
    [TAU,Z,N_grid_out,r_grid_out] = ndgrid(tau_d,z,1:length(indata),1:irr_pts); % z, t, 2-d grid to interpolate to
    Y_interp = interpn(T_grid,z_grid,N_grid,r_grid,Y_out,TAU,Z,N_grid_out,r_grid_out);

    %% propagate the probe through the sample

    for j = 1:length(z) - 1 % propagate the probe through the sample

        I_p(:,j + 1,:) = squeeze(I_p(:,j,:)) - ((squeeze(sigma_S_01*Y_interp(:,j,1,:) + sigma_S_12*Y_interp(:,j,2,:) + ...
            sigma_T_12*Y_interp(:,j,4,:)) + squeeze(alpha_2_e*I_p(:,j,:))) .* squeeze(I_p(:,j,:)))*delta_z;

    end

    I_p_end = I_p(:,end,:); % load the transmission at the end of the sample
    
    %% interpolate to excitation irradiance over the probe extent
    
    F = griddedInterpolant({tau_d,flip(rad_fac)},flip(squeeze(I_p_end),2));

    [X,Y,tau_d_3D] = meshgrid(x,y,tau_d);
    I_e_grid = exp(-((X - w_0_e*X0).^2 + Y.^2)/w_0_e^2);
    T_grid = zeros(length(x),length(y),length(tau_d));
    
    for m = 1:length(x)
        
        T_grid(m,:,:) = F(tau_d_3D(m,:,:),I_e_grid(m,:,:));
        
    end
    
    %% make plots of distributions if diagnostics requested

    if diagnostics % if diagnostic plots requested

        for n = 1:length(indata) % population plots

            figure(n)
            surf(TAU(:,1:end - 1,n,1)/tau_e,Z(:,1:end - 1,n,1)/L,squeeze(Y_interp(:,1:end - 1,n,1))/N0)
            axis tight
            xlabel('T (tau_e)');
            ylabel('Z (L)');
            zlabel([cell2mat(pop_names(n)),'/N_t']);
            title(['Population of ',cell2mat(pop_names(n)),' vs. time']);

        end

        figure(n + 1) % excitation irradiance plot
        plot(z/L,I_0_e_arr ./ I_0_e_norm)
        axis tight
        xlabel('Z (L)');
        ylabel('T (unitless)');
        title('Transmission of excitation vs. position');
        fac_str = strsplit(num2str(rad_fac,2),' ');
        I_ePlotlegend = legend(fac_str,'Location','southwest','Box','off');
        title(I_ePlotlegend,'I_e[r]/I_e_0');
                
        figure(n + 2) % probe transmission plot
        plot(tau_d/tau_e,squeeze(I_p_end))
        axis tight
        xlabel('T_d (tau_e)');
        ylabel('T (unitless)');
        title('Transmission of probe vs. delay');
        I_pPlotlegend = legend(fac_str,'Location','southwest','Box','off');
        title(I_pPlotlegend,'I_e[r]/I_e_0');

        drawnow
        
    end
    
    %% output

    % replicate the array (X,Y,T_d) for the time axis (X,Y,T,T_d)
    output = repmat(T_grid,[1,1,1,length(t_in)]);
    output = permute(output,[1,2,4,3]);
    
    if diagnostics

        BD_3_level_triplet = toc;
        disp([datestr(now),': BD_3_level_triplet - Calculation time: ',num2str(BD_3_level_triplet),' s.'])

    end
    
end