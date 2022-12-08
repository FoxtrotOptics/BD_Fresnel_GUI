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

function dN = excite_prop(t,y,z,...
    I_0_e,v_e,v_p,tau_e,E_photon,...
    sigma_S_01,sigma_S_12,sigma_T_12,...
    tau_ISC,tau_S_10,tau_S_21,tau_TS,tau_T_21,...
    alpha_2)

    I_e = I_e_func(z,t,I_0_e,v_e,v_p,tau_e);

    c = [sigma_S_01/E_photon*I_e*y(1),... % absorption from S_0 to S_1
        sigma_S_12/E_photon*I_e*y(2),... % absorption from S_1 to S_2
        sigma_T_12/E_photon*I_e*y(4),... % absorption from T_1 to T_2
        y(2)/tau_ISC,... % intersystem crossing from S_1 to T_1
        y(2)/tau_S_10,... % decay from S_1 to S_0
        y(3)/tau_S_21,... % decay from S_2 to S_1
        y(4)/tau_TS,... % decay from T_1 to S_0
        y(5)/tau_T_21,... % decay from T_2 to T_1
        alpha_2/(2*E_photon)*I_e^2]; % two photon absorption
    
    dN =    [-c(1) + c(5) + c(7) + c(9)         % N_S_0
            c(1) - c(2) - c(4) - c(5) + c(6);   % N_S_1
            c(2) - c(6) + c(9);                 % N_S_2
            -c(3) + c(4) - c(7) + c(8)          % N_T_1
            c(3) - c(8)];                   % N_T_2
 
end