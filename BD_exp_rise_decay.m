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

function[R_exp] = BD_exp_rise_decay(alpha_2_exp,n_2_exp,tau_exp_r,tau_exp_f,tau_p,T1,diagnostics,mech_name,Y_axis)

if (abs(n_2_exp) + abs(alpha_2_exp)) ~= 0 && tau_exp_r > 0 && tau_exp_f > 0 % only if parameters are valid

    t_exp_r =  tau_exp_r /tau_p; % normalized collisional rise time
    t_exp_f = tau_exp_f / tau_p; % normalized collisional fall time

    r_exp = (heaviside(T1) .* ((1 - exp(-heaviside(T1).*T1/t_exp_r)) .* exp(-heaviside(T1).*T1/t_exp_f))); % response function
        
    C_exp = 1/(t_exp_f^2/(t_exp_f + t_exp_r)); % normalization factor
   
    R_exp = (n_2_exp + alpha_2_exp*1i)*C_exp*r_exp; % scaled response function
    
else

    R_exp = 0*T1; % set the response to zero

end

if diagnostics % plot response function for diagnostics
    
    if n_2_exp ~= 0 % plot the refraction
    
        figure('Name',mech_name)
        plot(T1,real(R_exp));
        axis('tight');
        xlabel('t (tau_p)');
        ylabel(cell2mat(Y_axis(1)));
        title(mech_name);
        
    end
    
    if alpha_2_exp ~= 0 % plot the absorption
    
        figure('Name',mech_name)
        plot(T1,imag(R_exp));
        axis('tight');
        xlabel('t (tau_p)');
        ylabel(cell2mat(Y_axis(2)));
        title(mech_name);
        
    end
    
end