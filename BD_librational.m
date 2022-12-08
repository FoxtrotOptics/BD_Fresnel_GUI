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

%%
% INPUT

function[R_l_interp] = BD_librational(n_2l,omega_0,omega_std,tau_lf,tau_p,T1,diagnostics)

if  n_2l ~= 0 && omega_0 > 0 && omega_std > 0 && tau_lf > 0; % only if parameters are valid

    omega_0 = omega_0 * tau_p; % central Frequency of Gaussian distribution of quantum oscillators [s^-1]
    omega_std = omega_std * tau_p; % standard Deviation of Gaussian distribution of quantum oscillators [s^-1]
    tau_lf = tau_lf / tau_p; % librational fall time [s]
    
    omega_pts = 2^10;
    omega = linspace(0,3*omega_0,omega_pts); % frequency distribution for integration of librational response [s^-1]
    delta_omega = omega(2) - omega(1); % frequency step [s^-1]
    omega = omega + delta_omega; % shift so first element is not zero
      
    T_pts = 2^12;
    T2 = linspace(min(T1),5*tau_lf,T_pts);
    [Tl,omegal] = meshgrid(T2,omega); % 2D matricies of time and frequency
    
    g = exp(-(omegal - omega_0).^2/2/omega_std^2) - exp(-(omegal + omega_0).^2/2/omega_std^2);
    int_l = heaviside(Tl) ./ omegal .* sin(heaviside(Tl).*omegal.*Tl) .*exp(-heaviside(Tl).*Tl/tau_lf);
    r_l = trapz(omega,g .* int_l,1);
    
    C_l = 1/trapz(T2,r_l);
    
    R_l = n_2l*C_l*r_l; % normalize and scale
    R_l_interp = interp1(T2,R_l,T1); % interpolate to requested points in time

else

    R_l_interp = 0*T1; % set the response to zero

end

if diagnostics % plot response function for diagnostics
    
    figure('Name','Librational response function R_l[t]')
    plot(T1,R_l_interp);
    axis('tight');
    xlabel('t (tau_p)');
    ylabel('R_l[t] (m^2/W)');
    title('Librational response function R_l[t]');
    
end