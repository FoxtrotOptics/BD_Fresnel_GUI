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

% irradiance of probe as function of time and space
function output = I_e_func(z,t,I_0_e,v_e,v_p,tau_e)
    
    output = I_0_e*exp(-(t - z/v_e + z/v_p)^2/tau_e^2);
    
end