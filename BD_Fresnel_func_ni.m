%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   BD_Fresnel_Temporal_Dynamics.m          
%   
%   Initial version (1.0):    07/23/13 - Matthew Reichert
%   Latest revision (1.5):    03/09/17 - Manuel Ferdinandus
%
%   Calculcates cross-correlation signal of a Beam Deflection pump-probe 
%   experiment with NLR and 2PA.  Includes GVM. Uses Fresnel Propagation to
%   get to detector.  Calculations are symmetric in Y and T_d, thus only
%   simulate half those sp35116625aces to save time and memory.  Can't do the sme
%   for X or T
% 
%   This is an initial attempt at including time dynamics, where we at come
%   non-instatantaneou response function chi_n (see Beam Deflection
%   Derivation - Time Dynamics.docx).  
%
%   Expressions:    _p = pump
%                   _n = non-instantaneous component
%                   c  = collision induced
%                   v  = vibrational
%                   l  = librational
%                   d  = diffusive reorientation
%                   
%   Assumptions:    Zero GVD
%                   Linear absorption in the pump beam
%   
%   Notes:  Based on "Two-photon spectroscopy and analysis with a 
%           white-light continuum probe" by Negres et. al. (2002)
%           The [] after the description gives the units.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT

function[output] = BD_Fresnel_func_ni(params)

tic

% parse input parameters
lambda_p = params.lambda_p*(1e-9); % pump wavelength [m]
tau_p = params.tau_p*(1e-15); % pump temporal pulse width (FWHM) [s]
w_0_p = params.w_0p*(1e-6); % pump spot size (HW1/e^2 of I) [m]

lambda = params.lambda*(1e-9); % probe wavelength [m]
tau = params.tau*(1e-15); % probe temporal pulse width (HW1/e of I) [s]
w_0 = params.w_0*(1e-6); % probe spot size (HW1/e^2 or I) [m]

d = params.d*(1e-2); % distance from Sample to Detector [c]
X0 = params.x_0; % pump-probe x-axis displacement [w_p]

L = params.l*(1e-3); % thickness [m]

alpha = params.alpha*100; % absorption coefficient of probe
alpha_2d = params.alpha_2d*(1e-11); % effective degenerate two photon absorption of pump

sigma_p = params.sigma_p; % linear absorption parameter of pump

I_0_p = params.I_0_p; % peak irradiance of pump [W/m^2]

% group velocity calculation
ng_p = params.ng_p; % group index at pump
ng = params.ng; % group index at probe

rho = params.rho; % GVM Parameter

X_max = params.X_max; % maximum spatial extent in X
X_num = params.X_num; % number of points in X

T_min = params.T_min; % min range in T [tau_p]
T_max = params.T_max; % max range in T [tau_p]
T_num = params.T_num; % number of points in T [tau_p]

T_d_min = params.T_d_min; % minimum of Normalized Delay Vector [tau_p]
T_d_max = params.T_d_max; % maximum of Normalized Delay Vector [tau_p]
T_d_num = params.T_d_num; % number of points in T_d

%% CONSTANTS

c = 299792458; % speed of light [m/s]

%% CALCULATED PARAMETERS

% determine if we need to consider fabry-perot effects in the reflectivity
Delta_x_p = c*tau_p;
Delta_x_ratio = L/Delta_x_p;

% calculate absorption and reflection of pump
if Delta_x_ratio < 1 % if the sample length is longer than the pulse length
   
    disp([datestr(now),': BD_Frensel_func - Pulse length is longer than sample. Ensure you have corrected for multilayer effects.'])
    
end

sigma = alpha*L/2; % linear absorption parameter of probe

% normalize dimensions
t = tau/tau_p; % normalized Probe Pulse Width
W_0_p = w_0_p/w_0; % normalized Pump Radius
D = d/w_0; % normalized Distance from Sample to Detector
LAMBDA = lambda/w_0; % normalized Wavelength
K_0 = 2*pi/LAMBDA; % normalized Wavenumber


dX = abs(2*X_max/X_num); % differential unit of normalized spatial vector in X-direction [w_0]
dT = 2*abs(T_max + abs(rho))/T_num; % differential unit of normalized temporal vector [tau_p]
dT_d = abs(T_d_max - T_d_min)/T_d_num; % differential unit of normalized delay vector [tau_p]

%% INITIALIZATION

X = -X_max:dX:X_max; % normalized Spatial Vector in X-direction [w]
Y = -X_max:dX:X_max; % normalized Spatial Vector in Y-direction [w]

% adjust the temporal range for the beam walk through
if rho > 0 % this is valid only for rho /= 0
    
    T_L = T_min + heaviside(-rho)*rho;
    T_U = T_max + heaviside(rho)*rho;

else % if rho is zero no adjustment is required
    
    T_L = T_min;
    T_U = T_max;
    
end

T = T_L:dT:T_U; % normalized temporal vector [tau_p]
T_d = T_d_min:dT_d:T_d_max; % normalized temporal delay vector [tau_p]
[Y,X,T,T_d] = ndgrid(Y,X,T,T_d); % normalized temporal matrix [tau_p] (t,t_d)

% account for the saturable absorption
a_p_0 = exp(-((X + X0*W_0_p) .^2 + Y.^2) / W_0_p^2 ); % excitation spatial profile
sigma_p2 = I_0_p*a_p_0*alpha_2d*L/2; % effective degenerate two photon abosrption parameter of pump

% integration range should cover entire time pulse is in sample
T_fac = 5;

T1_max = T_d_max + 2*T_fac + abs(rho);
T1_min = -T1_max;
T1_num = 2^14 - 1;

dT1 = abs(T1_max - T1_min)/T1_num;
T1 = T1_min:dT1:T1_max;

% the array of responses for each mechanism
n_responses = length(params.mech);
R_n = zeros(n_responses,length(T1));
I_RESP = 0*T;
a_ESA = I_RESP + 1;
NI_resp = R_n;

%% REPSONSE FUNCTIONS

rho_min = 1e-4;
for i = 1:n_responses

    param = params.mech(i).param;

    switch cell2mat(params.mech(i).name)

        case 'None (both)'
        
        case 'Electronic (both)' % electronic

            alpha_2e = param(1).value*(1e-11); % ND-2PA Coefficient [m/W]
            n_2e = param(2).value*(1e-19); % ND-NLR [m^2/W]

            % calculate the response due to the electronic mechanism
            if rho < rho_min && sigma_p == 0
                
                I_RESP = I_0_p*(2*n_2e + 1i*alpha_2e)*exp(-(T + T_d).^2).*(-1 + sigma_p2);
            
            elseif rho < rho_min && sigma_p ~= 0

                A1 = exp(-2*sigma_p - (T + T_d).^2);
                A2 = sigma_p - sigma_p2 + sigma_p*coth(sigma_p);
                A3 = sinh(sigma_p)^2;
                
                I_RESP = I_0_p*(2*n_2e + 1i*alpha_2e)*(A1 .* A2 .* A3)/(2*sigma_p^2);
                
            elseif rho ~= 0 && sigma_p == 0
                
                A1 = 2*exp(-rho^2 - (T + T_d).^2);
                A2 = sigma_p2 .* (exp(rho^2) - exp(2*rho*(T + T_d)));
                A3 = rho*(-1 + 2*sigma_p2) - 2*sigma_p2 .* (T + T_d);
                A4 = erf(rho - (T + T_d)) + erf(T + T_d);
                
                I_RESP = I_0_p*(2*n_2e + 1i*alpha_2e)*((A1 .* A2) - (A3 .* A4)*sqrt(pi))/(2*rho^2);
                
            elseif rho ~= 0 && sigma_p ~= 0

                A1 = exp((sigma_p*(sigma_p + 2*rho*(-rho + (T + T_d))))/(rho^2)) .* (sigma_p - sigma_p2);
                A2 = erf(rho - (sigma_p/rho) - (T + T_d)) + erf((sigma_p/rho) + (T + T_d));
                A3 = exp((4*sigma_p*(sigma_p + rho*(-rho + (T + T_d))))/(rho^2)) .* sigma_p2;
                A4 = erf(rho - 2*(sigma_p/rho) - (T + T_d)) + erf(2*(sigma_p/rho) + (T + T_d));
                
                I_RESP = I_0_p*(2*n_2e + 1i*alpha_2e)*((A1 .* A2) + (A3 .* A4))*sqrt(pi)/(2*rho*sigma_p);

            end
            
        case 'Collisional (NLR)' % collisional

            n_2c = param(1).value*(1e-19); % collisional isotropic nonlinear refractive index [m^2/W]
            tau_cr = param(2).value*(1e-15); % collisional rise time [s]
            tau_cf = param(3).value*(1e-15); % collisional fall time [s]

            R_n(i,:) = BD_exp_rise_decay(0,n_2c,tau_cr,tau_cf,tau_p,T1,params.diagnostics,...
                'Collisional response',{'R_c[t] (m^2/W)',''}); % calculate the response

        case 'Librational (NLR)' % librational
            
            n_2l = param(1).value*(1e-19); % librational isotropic nonlinear refractive index [m^2/W]
            omega_0 = param(2).value*(1e12); % central frequency of Gaussian distribution of quantum oscillators [s^-1]
            omega_std = param(3).value*(1e12); % standard deviation of Gaussian distribution of quantum oscillators [s^-1]
            tau_lf = param(4).value*(1e-15); % librational fall time [s]

            R_n(i,:) = BD_librational(n_2l,omega_0,omega_std,tau_lf,tau_p,T1,params.diagnostics);
    
        case 'Reorientational (NLR)' % reorientational
            
            n_2d = param(1).value*(1e-19); % diffusive (reorientational) nonlinear refractive index [m^2/W] 
            tau_rr = param(2).value*(1e-15); % diffusive rise time [s]
            tau_rf = param(3).value*(1e-15); % diffusive fall time [s]

            R_n(i,:) = BD_exp_rise_decay(0,n_2d,tau_rr,tau_rf,tau_p,T1,params.diagnostics,...
                'Reorientational response',{'R_r[t] (m^2/W)',''}); % calculate the response
            
        case 'Vibrational (NLR)' % vibrational
            
            n_2v = param(1).value*(1e-19); % vibrational Isotropic nonlinear refractive index [m^2/W]
            w_0v = param(2).value*(1e12); % vibrational Frequency [rad/s]
            tau_vf = param(3).value*(1e-12); % vibrational fall time [s]

            R_n(i,:) = BD_vibrational(n_2v,w_0v,tau_vf,tau_p,T1,params.diagnostics); % calculate the response
            
        case 'ESA - 3 level + triplet (both)' % ESA
            
            alpha_2_e = param(1).value*(1e-11);
            N0 = param(2).value*(1e+22);
            sigma_S_01 = param(3).value*(1e-21);
            sigma_S_12 = param(4).value*(1e-21);
            sigma_T_12 = param(5).value*(1e-19);
            tau_S_10 = param(6).value*(1e-12);
            tau_S_21 = param(7).value*(1e-12);
            tau_ISC = param(8).value*(1e-12);
            tau_T_21 = param(9).value*(1e-12);
            tau_TS = param(10).value*(1e-12);
            
            a_ESA = BD_3_level_triplet(alpha_2_e,I_0_p,lambda_p,ng_p,tau_p,w_0_p,...
            ng,L,N0,X0,...
            sigma_S_01,sigma_S_12,sigma_T_12,...
            tau_ISC,tau_S_10,tau_S_21,tau_TS,tau_T_21,...
            T(1,1,:,1),squeeze(T_d(1,1,1,:))*tau_p,X(1,:,1,1)*w_0,Y(:,1,1,1)*w_0,...
            params.diagnostics);
        
        case 'Exponential rise/decay (NLA)' % generic exponential rise/decay function for NLA

            alpha_2exp = param(1).value*(1e-11); % generic two-photon absorption coefficient [m/W]
            tau_expr = param(2).value*(1e-15); % rise time [s]
            tau_expf = param(3).value*(1e-15); % fall time [s]

            R_n(i,:) = BD_exp_rise_decay(alpha_2exp,0,tau_expr,tau_expf,tau_p,T1,params.diagnostics,...
                'Exponential rise/decay (NLA) response',{'','A_exp[t] (m/W)'}); % calculate the response
            
        case 'Exponential rise/decay (NLR)' % generic exponential rise/decay function for NLR

            n_2exp = param(1).value*(1e-19); % generic nonlinear index of refraction [m^2/W]
            tau_expr = param(2).value*(1e-15); % rise time [s]
            tau_expf = param(3).value*(1e-15); % fall time [s]

            R_n(i,:) = BD_exp_rise_decay(0,n_2exp,tau_expr,tau_expf,tau_p,T1,params.diagnostics,...
                'Exponential rise/decay (NLR) response',{'R_exp[t] (m^2/W)',''}); % calculate the response
            
        case 'Exponential rise/decay (both)' % generic exponential rise/decay function for NLR

            alpha_2exp = param(1).value*(1e-11); % generic two-photon absorption coefficient [m/W]
            n_2exp = param(2).value*(1e-19); % generic nonlinear index of refraction [m^2/W]
            tau_expr = param(3).value*(1e-15); % rise time [s]
            tau_expf = param(4).value*(1e-15); % fall time [s]

            R_n(i,:) = BD_exp_rise_decay(alpha_2exp,n_2exp,tau_expr,tau_expf,tau_p,T1,params.diagnostics,...
                'Exponential rise/decay (both) response',{'R_exp[t] (m^2/W)','A_exp[t] (m/W)'}); % calculate the response
            
    end
    
end

%% INDEX CHANGE DUE TO RESPONSE FUNCTIONS

is_resp = find(sum(R_n,2));
NI_RESP = zeros(size(T,3),size(T_d,4)); % T,T_d,number of independent responses

if ~isempty(is_resp)
    
    R_n = R_n(is_resp,:); % extract out the response functions

    I_p = I_0_p*exp(-(T1).^2); % normalized temporal irradiance

    % ensure that over a long time scale the gaussian turns into a delta function
    if isequal(I_p,0*T1)

        I_p(round(T1_num/2)) = 1;

    end

    % the total response function is the response convolved with the pump
    NI_resp = 0*R_n;
    for i = 1:size(R_n,1)

        NI_resp(i,:) = conv(R_n(i,:),I_p,'same')*dT1;
        % note Delta_n is coded into the real part of Resp
        % Delta_alpha is coded into the imaginary part of Resp

    end

    % determine the upper and lower ranges for the integration
    upp_base = squeeze(T(1,1,:,:) + T_d(1,1,:,:));
    low_base = squeeze((T(1,1,:,:) + T_d(1,1,:,:)) - rho);

    Resp_tot = sum(NI_resp,1); % the total index change due to all mechanisms
    Resp_tot_diag = griddedInterpolant(T1,Resp_tot); % interpolate the response function
    
    sigma_p20 = I_0_p*alpha_2d*L/2; % peak effective TPA of pump
    
    % calculate the attenuation factors
    zoom_pts = 2^8;
    att_arr = linspace(0,1,zoom_pts);
    
    if sigma_p2 == 0
       
        atten_fac = exp(-2*sigma_p*att_arr);
        
    else
        
        B1 = exp(-4*sigma_p*att_arr);
        B2 = exp(2*sigma_p*att_arr)*(sigma_p - sigma_p20) + sigma_p20;
        atten_fac = sqrt((B1 .* B2)/sigma_p);
        
    end

    for i = 1:size(NI_RESP,1)

        for j = 1:size(NI_RESP,2)

            if rho ~= 0 % integrate over the pulse overlap if rho is not zero

                T1_zoomed = linspace(low_base(i,j),upp_base(i,j),zoom_pts); % extract out the region of integration
                Resp_tot_zoomed = Resp_tot_diag(T1_zoomed) .* atten_fac; % calculate the response over the overlap
                NI_RESP(i,j) = trapz(T1_zoomed,Resp_tot_zoomed/rho); % integrate over the pulse overlap

            else % otherwise the response is just the response at one point
        
                NI_RESP(i,j) = Resp_tot_diag(low_base(i,j));
    
            end

        end

    end
        
end

%% COMBINE INSTANTANEOUS AND NON-INSTANTANEOUS RESPONSES

NI_RESP = permute(repmat(NI_RESP,[1,1,size(Y,2),size(X,2)]),[4,3,1,2]);

ALL_RESP = cat(5,I_RESP,NI_RESP); % combine all responses
TOT_RESP = sum(ALL_RESP,5);

%define non-instantaneous nonlinear parameter
eta = (2*pi/lambda)*L*real(TOT_RESP);
Gamma = L*imag(TOT_RESP);

CHI_RESP = eta + 1i*Gamma;

%% DIAGNOSTICS

if params.diagnostics %&& exist('NI_resp','var') % plot the response functions if required
    
    disp([datestr(now),': BD_Fresnel_func - generating diagnostic plots.'])
    
    if exist('n_2e','var') && exist('alpha_2e','var') && exist('NI_resp','var') % add the instantaneous response if present
    
        I_resp = (n_2e + alpha_2e*1i)*I_0_p*exp(-T1.^2); % calculate the instantaneous response
        ALL_resp_diag = cat(1,NI_resp,I_resp); % combine with the nonistantaneous response
    
    else
        
        ALL_resp_diag = NI_resp; % combine with the nonistantaneous response
    
    end
    
    Resp_tot_diag = sum(ALL_resp_diag,1); % calculate the sum
    
    NI_mech = figure('Name','Response from all mechanisms');
    
    T_plot_min = -3;
    T_plot_max = T1_max;
    
    subplot(2,1,1);
    hold on;
    plot(T1,[real(ALL_resp_diag)',real(Resp_tot_diag)']); % plot all mechanisms
    hold off
    axis('tight');
    xlim([T_plot_min,T_plot_max]);
    xlabel('t (tau_p)');
    ylabel('Delta n (unitless)');
    title('Index change from all mechanisms');
    
    subplot(2,1,2);
    hold on;
    plot(T1,[imag(ALL_resp_diag)',imag(Resp_tot_diag)']); % plot all mechanisms
    hold off
    axis('tight');
    xlim([T_plot_min,T_plot_max]);
    xlabel('t (tau_p)');
    ylabel('Delta alpha (m^-1)');
    title('Absorption change from all mechanisms');
    
    NI_mech_pos = get(NI_mech,'Position');
    set(NI_mech,'Position',[NI_mech_pos(1),NI_mech_pos(2) - NI_mech_pos(4)*0.75,NI_mech_pos(3),NI_mech_pos(4)*1.75]);
    
    drawnow;

elseif params.diagnostics && ~exist('NI_resp','var')
    
    disp([datestr(now),': BD_Fresnel_func - diagnostic plots only generated for noninstantaneous responses - no plots generated.'])
    
end

%% PROPAGATION TO DETECTOR

% calculate the field at the back of the sample

spacePart = exp(-X.^2 - Y.^2); % spatial part of field

a_out = spacePart .* exp(-sigma -((-rho + T).^2 / t^2)/2 + (1i*(a_p_0.^2) .* CHI_RESP)) .* a_ESA; % field at back of detector

a_det = fftshift(fftshift(fft2(a_out .* exp(1i*K_0/(2*D)*((X).^2+(Y).^2))),2),1); % fresnel propagate to detector

% integrate over left and right sides
E_left = squeeze(trapz(trapz(trapz(abs(a_det(:,1:ceil(end/2),:,:)).^2,1),2),3));
E_right = squeeze(trapz(trapz(trapz(abs(a_det(:,ceil(end/2):end,:,:)).^2,1),2),3));
E_tot = E_left + E_right;

%% PACKAGE OUTPUT

output.T = E_tot / E_tot(1);
output.tau_d = squeeze(T_d(1,1,1,:))*tau_p*1e15;
output.DeltaEoE = (E_left - E_right) ./ E_tot;
output.Delta_n = output.DeltaEoE / ((2*pi/lambda)*L*sqrt(2/exp(1))*(w_0 / w_0_p)*(2 / sqrt(pi)));

BD_Fresnel_func_time = toc;
if params.diagnostics
    
    disp([datestr(now),': BD_Fresnel_func - calculation time: ',num2str(BD_Fresnel_func_time),' s.'])

end

end
