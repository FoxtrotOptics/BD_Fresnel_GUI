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

function[R_v] = BD_vibrational(n_2v,w_0v,tau_vf,tau_p,T1,diagnostics)

if  n_2v ~= 0 && w_0v > 0 && tau_vf > 0; % only if parameters are valid

    w_0v = w_0v * tau_p; % central Frequency of Gaussian distribution of quantum oscillators [s^-1]
    tau_vf = tau_vf / tau_p; % vibrational fall time [s]
    r_v = heaviside(T1) .* sin(heaviside(T1) .* w_0v .* T1) .* exp(-heaviside(T1) .* T1 / tau_vf);

    if max(T1)/tau_vf > 10
        
        C_v = 1/trapz(T1,r_v);
    
    else
        
        C_v = 1/(w_0v / (w_0v^2 + (1/tau_vf)^2)); % normalization factor
    
    end

    R_v = n_2v*C_v*r_v; % normalize and scale the response

else

    R_v = 0*T1; % set the response to zero

end

if diagnostics % plot response function for diagnostics
    
    figure('Name','Vibrational response function R_v[t]')
    plot(T1,R_v);
    axis('tight');
    xlabel('t (tau_p)');
    ylabel('R_v[t] (m^2/W)');
    title('Vibrational response function R_v[t]');
    
end