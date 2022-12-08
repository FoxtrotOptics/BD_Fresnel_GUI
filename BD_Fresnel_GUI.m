%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Beam deflection fitting program                
%
% Version: 2.2
% Date: 03/15/17
%
% GUI calculator for the beam deflection signal calculator. This program
% allows the user to import files, specify measurement parameters, fit
% beam deflection data and export it. It takes in a data file that has
% been converted into DeltaEoE and T using the lock-in conversion program.
% It then outputs the fit file that contains the original data and the
% fitted data. It also outputs a log file that specifies all of the
% parameters used to do the fit and some ancillary data.
%
% The core algorithm the program is based on the work of Dr. Matthew
% Reichert while at the University of Central Florida in 2014. The 
% algorithm is described in the following papers:
%
% 1: Reichert, M., et al. (2014). "Temporal, spectral, and polarization
% dependence of the nonlinear optical response of carbon disulfide
% Optica 1(6): 436.
%
% 2: Reichert, M., et al. (2014). "Temporal, spectral, and polarization 
% dependence of the nonlinear optical response of carbon disulfide:
% supplementary material." Optica 1(6): 436.
%
% The beam deflection technique is based on the work of Dr. Manuel
% Ferdinandus while at the University of Central Florida. The technique 
% is described in the following paper:
%
% 3: Ferdinandus, M. R., et al. (2013). "Beam deflection measurement of 

% time and polarization resolved ultrafast nonlinear refraction." 
% Opt Lett 38(18): 3518-3521.
%
% If you are using the products of this program in your publication please
% cite the papers listed above.
%
%     Copyright (C) 2017  Manuel Ferdinandus (mrf@knights.ucf.edu)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = BD_Fresnel_GUI(varargin)
% BD_Fresnel_GUI MATLAB code for BD_Fresnel_GUI.fig
%      BD_Fresnel_GUI, by itself, creates a new BD_Fresnel_GUI or raises the existing
%      singleton*.
%
%      H = BD_Fresnel_GUI returns the handle to a new BD_Fresnel_GUI or the handle to
%      the existing singleton*.
%
%      ('CALLBACK',hObject,eventData,handles,...) calls the local BD_Fresnel_GUI
%      function named CALLBACK in BD_Fresnel_GUI.M with the given input arguments.
%
%      BD_Fresnel_GUI('Property','Value',...) creates a new BD_Fresnel_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BD_Fresnel_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BD_Fresnel_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BD_Fresnel_GUI

% Last Modified by GUIDE v2.5 09-Mar-2017 21:25:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BD_Fresnel_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BD_Fresnel_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BD_Fresnel_GUI is made visible.
function BD_Fresnel_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BD_Fresnel_GUI (see VARARGIN)

% UIWAIT makes BD_Fresnel_GUI wait for user response (see UIRESUME)
% uiwait(handles.zfit_manualfit);

clc % clear command window
disp([datestr(now),': BD_Fresnel - Begin BD_Fresnel_GUI.'])

% Choose default command line output for BD_Fresnel_GUI
handles.output = hObject;

% retrieve old settings
current_directory = pwd; % the current directory
handles.params.save_filepath = [current_directory,'\last_settings.mat'];

% get the last path used if available
if exist(handles.params.save_filepath,'file') == 2
    
    last = open(handles.params.save_filepath);
    handles.params = last.settings;
    clearvars last
    
    % check the file paths
    if ~exist([handles.params.pathname,handles.params.filename],'file') == 2
    
        % file and path name arrays for import
        handles.params.filename = '\*.txt';
        handles.params.pathname = current_directory;
    
    end
    
    disp([datestr(now),': BD_Fresnel - BD_Fresnel - Previous settings loaded from ',handles.params.save_filepath,'.']);

else % if the last path does not exist then set the parameters manually
    
    disp([datestr(now),': BD_Fresnel - BD_Fresnel - Previous settings not found, loading defaults (quartz).'])
    
    % Set the default parameters
    handles.params.E_p = 1.36; % pump beam
    handles.params.tau_p_FWHM = 55;
    handles.params.w_0p = 246;
    handles.params.lambda_p = 800;

    handles.params.tau_FWHM = 60; % probe beam
    handles.params.w_0 = 43;
    handles.params.lambda = 650;

    handles.params.R_p = 0;
    handles.params.alpha_p = 0;
    handles.params.n_p = 1.4533; % at pump wavelength
    handles.params.dndl_p = -0.017284;
    handles.params.alpha_2d = 0;

    handles.params.R = 0;
    handles.params.alpha = 0;
    handles.params.n = 1.4565; % at probe wavelength
    handles.params.dndl = -0.027158;
 
    % set up the mechanisms
    handles.params.mechanism_names = {'None (both)',...
        'Electronic (both)',...
        'Collisional (NLR)',...
        'Librational (NLR)',...
        'Reorientational (NLR)',...
        'Vibrational (NLR)',...
        'ESA - 3 level + triplet (both)',...
        'Exponential rise/decay (NLA)',...
        'Exponential rise/decay (NLR)',...
        'Exponential rise/decay (both)'};
    handles.params.n_mech = 7;
    
    % set up the structure for the mechanisms
    n_params = 10;
    param.name = '';
    param.units = '';
    param.value = 0;
    param = repmat(param,n_params);
    
    mech.mechanism = 1;
    mech.name = handles.params.mechanism_names(mech.mechanism);
    mech.param = param;
    handles.params.mech = repmat(mech,handles.params.n_mech);    

    handles.params.l = 1; % system parameters
    handles.params.d = 15;
    handles.params.shift = 0;
    handles.params.x_0 = 1/2;

    handles.params.diagnostics = 0;
    
    handles.params.data_imported = 0; % whether valid data has been imported
    handles.params.data_fit = 0; % whether valid data has been fit

    handles.params.tau_p = handles.params.tau_p_FWHM /(2*sqrt(log(2))); % convert between FWHM to HW1/e [s]
    handles.params.tau = handles.params.tau_FWHM /(2*sqrt(log(2)));

    handles.params.X_max = 5; % maximum of normalized spatial vector in X/Y-direction [w]
    handles.params.X_num = 50; % number of points in spatial direction

    handles.params.T_min = -5; % minimum of normalized temporal vector [tau_p]
    handles.params.T_max = 5; % maximum of normalized temporal vector [tau_p]
    handles.params.T_num = 100; % number of points in temporal vector

    handles.params.axesFontSize = 10; % plot options
    handles.params.dataMarkerSize = 5;
    
    % file and path name arrays for import
    handles.params.filename = '\*.txt';
    handles.params.pathname = current_directory;
           
end

% load default values into text boxes
set(handles.et_E_p,'String',num2str(handles.params.E_p));
set(handles.et_tau_p_FWHM,'String',num2str(handles.params.tau_p_FWHM));
set(handles.et_w_0p,'String',num2str(handles.params.w_0p));
set(handles.et_lambda_p,'String',num2str(handles.params.lambda_p));

set(handles.et_tau_FWHM,'String',num2str(handles.params.tau_FWHM));
set(handles.et_lambda,'String',num2str(handles.params.lambda));
set(handles.et_w_0,'String',num2str(handles.params.w_0));

set(handles.et_R_p,'String',num2str(handles.params.R_p));
set(handles.et_alpha_p,'String',num2str(handles.params.alpha_p));
set(handles.et_n_p,'String',num2str(handles.params.n_p));
set(handles.et_dndl_p,'String',num2str(handles.params.dndl_p));
set(handles.et_alpha_2d,'String',num2str(handles.params.alpha_2d));

set(handles.et_R,'String',num2str(handles.params.R));
set(handles.et_alpha,'String',num2str(handles.params.alpha));
set(handles.et_n,'String',num2str(handles.params.n));
set(handles.et_dndl,'String',num2str(handles.params.dndl));

set(handles.et_l,'String',num2str(handles.params.l));
set(handles.et_d,'String',num2str(handles.params.d));
set(handles.et_shift,'String',num2str(handles.params.shift));
set(handles.et_x_0,'String',num2str(handles.params.x_0));

handles = calc_beam_params(hObject,handles);

% other parameters
set(handles.st_rho_val,'String',num2str(handles.params.rho,'%6.3f'));

for i = 1:handles.params.n_mech % for all the input boxes

    % set up the input boxes
    set(eval(['handles.pu_mech_',num2str(i)]),'String',handles.params.mechanism_names);
    set(eval(['handles.pu_mech_',num2str(i)]),'Value',handles.params.mech(i).mechanism);

end

% Disable the Fit and Export buttons
set(handles.pb_Export,'Enable','off');
set(handles.pb_Fit,'Enable','off');

% set diagnostics button
set(handles.cb_diagnostics,'Value',handles.params.diagnostics);

% Set up the transmission plot
xlabel(handles.T_plot,'tau_d (ps)','FontSize',handles.params.axesFontSize); % transmission plot
ylabel(handles.T_plot,'T (%)','FontSize',handles.params.axesFontSize);

% set up the deflection plot
xlabel(handles.DeltaEoE_plot,'tau_d (ps)','FontSize',handles.params.axesFontSize); % DeltaEoE plot
ylabel(handles.DeltaEoE_plot,'DeltaEoE (%)','FontSize',handles.params.axesFontSize);

guidata(hObject,handles); % Update the handles


% --- Outputs from this function are returned to the command line.
function varargout = BD_Fresnel_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function et_E_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_E_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_E_p as text
%        str2double(get(hObject,'String')) returns contents of et_E_p as a double

temp = str2double(get(handles.et_E_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.E_p = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - BD_Fresnel - Value of E_p updated to ',num2str(handles.params.E_p),' uJ.'])
    
else
    
    disp([datestr(now),': BD_Fresnel - BD_Fresnel - E_p must be a real number greater than 0.'])
    handles.params.E_p = 0; % Update the variable
        
end

set(handles.et_E_p,'String',num2str(handles.params.E_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_E_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_E_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_tau_p_FWHM_Callback(hObject, eventdata, handles)
% hObject    handle to et_tau_p_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_tau_p_FWHM as text
%        str2double(get(hObject,'String')) returns contents of et_tau_p_FWHM as a double

temp = str2double(get(handles.et_tau_p_FWHM,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.tau_p_FWHM = abs(temp); % Update the variable
    handles.params.tau_p = handles.params.tau_p_FWHM/(2*sqrt(log(2))); % convert to HW/1e^2M
    disp([datestr(now),': BD_Fresnel - Value of tau_p_FWHM updated to ',num2str(handles.params.tau_p_FWHM),' fs.'])
    
else
    
    disp([datestr(now),': BD_Fresnel - tau_p_FWHM must be a real number greater than 0.'])
    handles.params.tau_p_FWHM = 0; % Update the variable
    
end

set(handles.et_tau_p_FWHM,'String',num2str(handles.params.tau_p_FWHM)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_tau_p_FWHM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_tau_p_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_w_0p_Callback(hObject, eventdata, handles)
% hObject    handle to et_w_0p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_w_0p as text
%        str2double(get(hObject,'String')) returns contents of et_w_0p as a double

temp = str2double(get(handles.et_w_0p,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.w_0p = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of w_0p updated to ',num2str(handles.params.w_0p),' um.'])

else
    
    disp([datestr(now),': BD_Fresnel - w_0p must be a real number greater than 0.'])
    handles.params.w_0p = 0; % Update the variable
    
end

set(handles.et_w_0p,'String',num2str(handles.params.w_0p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_w_0p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_w_0p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_lambda_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_lambda_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_lambda_p as text
%        str2double(get(hObject,'String')) returns contents of et_lambda_p as a double

temp = str2double(get(handles.et_lambda_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.lambda_p = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of lambda_p updated to ',num2str(handles.params.lambda_p),' nm.'])
    
else
    
    disp([datestr(now),': BD_Fresnel - lambda_p must be a real number greater than 0.'])
    handles.params.lambda_p = 0; % Update the variable
    
end

set(handles.et_lambda_p,'String',num2str(handles.params.lambda_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_lambda_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_lambda_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_tau_FWHM_Callback(hObject, eventdata, handles)
% hObject    handle to et_tau_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_tau_FWHM as text
%        str2double(get(hObject,'String')) returns contents of et_tau_FWHM as a double

temp = str2double(get(handles.et_tau_FWHM,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.tau_FWHM = abs(temp); % Update the variable
    handles.params.tau = handles.params.tau_FWHM/(2*sqrt(log(2))); % convert to HW1/e^2M
    disp([datestr(now),': BD_Fresnel - Value of tau_FWHM updated to ',num2str(handles.params.tau_FWHM),' fs.'])
    
else
    
    disp([datestr(now),': BD_Fresnel - tau_FWHM must be a real number greater than 0.'])
    handles.params.tau_FWHM = 0; % Update the variable
    
end

set(handles.et_tau_FWHM,'String',num2str(handles.params.tau_FWHM)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid
    
else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_tau_FWHM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_tau_FWHM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_w_0_Callback(hObject, eventdata, handles)
% hObject    handle to et_w_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_w_0 as text
%        str2double(get(hObject,'String')) returns contents of et_w_0 as a double

temp = str2double(get(handles.et_w_0,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.w_0 = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of w_0 updated to ',num2str(handles.params.w_0),' um.'])

else
    
    disp([datestr(now),': BD_Fresnel - w_0 must be a real number greater than 0.'])
    handles.params.w_0 = 0; % Update the variable
    
end

set(handles.et_w_0,'String',num2str(handles.params.w_0)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_w_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_w_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to et_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_lambda as text
%        str2double(get(hObject,'String')) returns contents of et_lambda as a double

temp = str2double(get(handles.et_lambda,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.lambda = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of lambda updated to ',num2str(handles.params.lambda),' nm.'])

else
    
    disp([datestr(now),': BD_Fresnel - lambda must be a real number greater than 0.'])
    handles.params.lambda = 0; % Update the variable
    
end

set(handles.et_lambda,'String',num2str(handles.params.lambda)) % update the display
handles = calc_beam_params(hObject,handles);

% check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_R_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_R_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_R_p as text
%        str2double(get(hObject,'String')) returns contents of et_R_p as a double

temp = str2double(get(handles.et_R_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp >= 0 && temp < 1 % if the value is a number between 0 and 1

    handles.params.R_p = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of R_p updated to ',num2str(handles.params.R_p),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - R_p must be a real number greater than 0 and less than 1.'])
    handles.params.R_p = 0; % Update the variable
    
end

set(handles.et_R_p,'String',num2str(handles.params.R_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_R_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_R_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_R_Callback(hObject, eventdata, handles)
% hObject    handle to et_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_R as text
%        str2double(get(hObject,'String')) returns contents of et_R as a double

temp = str2double(get(handles.et_R,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp >= 0 && temp < 1 % if the value is a number between 0 and 1

    handles.params.R = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of R updated to ',num2str(handles.params.R),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - R must be a real number greater than 0 and less than 1.'])
    handles.params.R = 0; % Update the variable
    
end

set(handles.et_R,'String',num2str(handles.params.R)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to et_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_alpha as text
%        str2double(get(hObject,'String')) returns contents of et_alpha as a double

temp = str2double(get(handles.et_alpha,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp >= 0 % if the value is a number between 0 and 1

    handles.params.alpha = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of alpha updated to ',num2str(handles.params.alpha),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - alpha must be a real number greater than 0.'])
    handles.params.alpha = 0; % Update the variable
    
end

set(handles.et_alpha,'String',num2str(handles.params.alpha)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_alpha_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_alpha_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_alpha_p as text
%        str2double(get(hObject,'String')) returns contents of et_alpha_p as a double

temp = str2double(get(handles.et_alpha_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp >= 0 % if the value is a number between 0 and 1

    handles.params.alpha_p = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of alpha_p updated to ',num2str(handles.params.alpha_p),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - alpha_p must be a real number greater than 0.'])
    handles.params.alpha_p = 0; % Update the variable
    
end

set(handles.et_alpha_p,'String',num2str(handles.params.alpha_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_alpha_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_alpha_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_n_Callback(hObject, eventdata, handles)
% hObject    handle to et_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_n as text
%        str2double(get(hObject,'String')) returns contents of et_n as a double

temp = str2double(get(handles.et_n,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp > 0 % if the value is a number and between 0 and 1

    handles.params.n = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of n updated to ',num2str(handles.params.n),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - n must be a real number greater than or equal to 1.'])
    handles.params.n = 0; % Update the variable
    
end

set(handles.et_n,'String',num2str(handles.params.n)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_dndl_Callback(hObject, eventdata, handles)
% hObject    handle to et_dndl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_dndl as text
%        str2double(get(hObject,'String')) returns contents of et_dndl as a double

temp = str2double(get(handles.et_dndl,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.dndl = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of dndl updated to ',num2str(handles.params.dndl),' um^-1.'])

else
    
    disp([datestr(now),': BD_Fresnel - dndl must be a real number.'])
    handles.params.dndl = 0; % Update the variable
    
end

set(handles.et_dndl,'String',num2str(handles.params.dndl)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_dndl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_dndl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_n_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_n_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_n_p as text
%        str2double(get(hObject,'String')) returns contents of et_n_p as a double

temp = str2double(get(handles.et_n_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) && temp > 0 % if the value is a number greater than 0

    handles.params.n_p = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of n_p updated to ',num2str(handles.params.n_p),'.'])

else
    
    disp([datestr(now),': BD_Fresnel - n_p must be a real number greater than or equal to 1.'])
    handles.params.n_p = 0; % Update the variable
    
end

set(handles.et_n_p,'String',num2str(handles.params.n_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_n_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_n_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_dndl_p_Callback(hObject, eventdata, handles)
% hObject    handle to et_dndl_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_dndl_p as text
%        str2double(get(hObject,'String')) returns contents of et_dndl_p as a double

temp = str2double(get(handles.et_dndl_p,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.dndl_p = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of dndl_p updated to ',num2str(handles.params.dndl_p),' um^-1.'])   

else
    
    disp([datestr(now),': BD_Fresnel - dndl_p must be a real number.'])
    handles.params.dndl_p = 0; % Update the variable
    
end

set(handles.et_dndl_p,'String',num2str(handles.params.dndl_p)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_dndl_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_dndl_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_alpha_2d_Callback(hObject, eventdata, handles)
% hObject    handle to et_alpha_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_alpha_2d as text
%        str2double(get(hObject,'String')) returns contents of et_alpha_2d as a double

temp = str2double(get(handles.et_alpha_2d,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.alpha_2d = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of alpha_2d updated to ',num2str(handles.params.alpha_2d),' e-11 m/W.'])   

else
    
    disp([datestr(now),': BD_Fresnel - alpha_2d must be a real number.'])
    handles.params.alpha_2d = 0; % Update the variable
    
end

set(handles.et_alpha_2d,'String',num2str(handles.params.alpha_2d)) % update the display
handles = calc_beam_params(hObject,handles);

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_alpha_2d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_alpha_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pu_mech_1.
function pu_mech_1_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_1

handles.params.mech(1).mechanism = get(handles.pu_mech_1,'Value');
handles.params.mech(1) = set_params_GUI(handles.params.mech(1),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_1_params.
function pb_mech_1_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_1_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(1) = set_params_GUI(handles.params.mech(1),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_2.
function pu_mech_2_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_2

handles.params.mech(2).mechanism = get(handles.pu_mech_2,'Value');
handles.params.mech(2) = set_params_GUI(handles.params.mech(2),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_2_params.
function pb_mech_2_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_2_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(2) = set_params_GUI(handles.params.mech(2),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_3.
function pu_mech_3_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_3

handles.params.mech(3).mechanism = get(handles.pu_mech_3,'Value');
handles.params.mech(3) = set_params_GUI(handles.params.mech(3),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_3_params.
function pb_mech_3_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_3_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(3) = set_params_GUI(handles.params.mech(3),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_3.
function pu_mech_4_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        pu_mech_4

handles.params.mech(4).mechanism = get(handles.pu_mech_4,'Value');
handles.params.mech(4) = set_params_GUI(handles.params.mech(4),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_3_params.
function pb_mech_4_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_4_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(4) = set_params_GUI(handles.params.mech(4),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_5.
function pu_mech_5_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_5

handles.params.mech(5).mechanism = get(handles.pu_mech_5,'Value');
handles.params.mech(5) = set_params_GUI(handles.params.mech(5),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_5_params.
function pb_mech_5_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_5_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(5) = set_params_GUI(handles.params.mech(5),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_6.
function pu_mech_6_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_6

handles.params.mech(6).mechanism = get(handles.pu_mech_6,'Value');
handles.params.mech(6) = set_params_GUI(handles.params.mech(6),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_6_params.
function pb_mech_6_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_6_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(6) = set_params_GUI(handles.params.mech(6),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes on selection change in pu_mech_7.
function pu_mech_7_Callback(hObject, eventdata, handles)
% hObject    handle to pu_mech_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pu_mech_7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pu_mech_7

handles.params.mech(7).mechanism = get(handles.pu_mech_7,'Value');
handles.params.mech(7) = set_params_GUI(handles.params.mech(7),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function pu_mech_7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pu_mech_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_mech_7_params.
function pb_mech_7_params_Callback(hObject, eventdata, handles)
% hObject    handle to pb_mech_7_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.params.mech(7) = set_params_GUI(handles.params.mech(7),handles.params.mechanism_names);

guidata(hObject,handles); % Update the handles


function et_l_Callback(hObject, eventdata, handles)
% hObject    handle to et_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_l as text
%        str2double(get(hObject,'String')) returns contents of et_l as a double

temp = str2double(get(handles.et_l,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.l = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of l updated to ',num2str(handles.params.l),' mm.'])
    
else
    
    disp([datestr(now),': BD_Fresnel - l must be a real number greater than 0.'])
    handles.params.l = 0; % Update the variable
    
end

set(handles.et_l,'String',num2str(handles.params.l)) % update the display
guidata(hObject,handles); % Update the handles

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end


% --- Executes during object creation, after setting all properties.
function et_l_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_d_Callback(hObject, eventdata, handles)
% hObject    handle to et_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temp = str2double(get(handles.et_d,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.d = abs(temp); % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of d updated to ',num2str(handles.params.d),' cm.'])

else
    
    disp([datestr(now),': BD_Fresnel - d must be a real number greater than 0.'])
    handles.params.d = 0; % Update the variable
    
end

set(handles.et_d,'String',num2str(handles.params.d)) % update the display
guidata(hObject,handles); % Update the handles

%check to see if the inputs are valid
if inputs_valid(handles.params)
    
    set(handles.pb_Fit,'Enable','on'); % enable the button if valid

else
    
    set(handles.pb_Fit,'Enable','off'); % disable the button if invalid
    
end


% --- Executes during object creation, after setting all properties.
function et_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_shift_Callback(hObject, eventdata, handles)
% hObject    handle to et_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_shift as text
%        str2double(get(hObject,'String')) returns contents of et_shift as a double

temp = str2double(get(handles.et_shift,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.shift = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of shift updated to ',num2str(handles.params.shift),' fs.'])

else
    
    disp([datestr(now),': BD_Fresnel - Shift must be a real number.'])
    handles.params.shift = 0; % Update the variable
    
end

% update display
set(handles.et_shift,'String',num2str(handles.params.shift)) % update the display
handles.data.tau_d_shifted = handles.data.tau_d + handles.params.shift; % the new shifted x axis
refresh_plots(hObject,handles);
    
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function et_x_0_Callback(hObject, eventdata, handles)
% hObject    handle to et_x_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_x_0 as text
%        str2double(get(hObject,'String')) returns contents of et_x_0 as a double

temp = str2double(get(handles.et_x_0,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.params.x_0 = temp; % Update the variable
    disp([datestr(now),': BD_Fresnel - Value of x_0 updated to ',num2str(handles.params.shift),' w_0_p.'])

else
    
    disp([datestr(now),': BD_Fresnel - x_0 must be a real number.'])
    handles.params.x_0 = 0; % Update the variable
    
end

% update display
set(handles.et_x_0,'String',num2str(handles.params.x_0)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_x_0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_x_0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_Import.
function pb_Import_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp([datestr(now),': BD_Fresnel - Begin file import.'])

[handles.params.filename,handles.params.pathname] = uigetfile(strcat(handles.params.pathname,handles.params.filename),'Import transmission/deflection data'); % get the file path

% import and update display
if ~strcmp(num2str(handles.params.filename),'0') % only if the filename is valid
    
    handles.params.import_filepath = [handles.params.pathname,handles.params.filename]; % construct the file path

    % Import and parse the data file
    DATA = importdata(handles.params.import_filepath);

    % check to see if the data is valid
    arraySize = size(DATA.data);
    if arraySize(2) >= 3 && sum(sum(isnan(DATA.data))) == 0

        tau_d_col = 1; % the column of delay
        DeltaEoE_col = 2; % the column of deflection
        T_col = 3; % the column of transmission
        
        % Extract the data and populate the structure
        handles.data.tau_d = DATA.data(:,tau_d_col); % extract the data from the array and put it in the right structure
        handles.data.tau_d_shifted = handles.data.tau_d + handles.params.shift;
        handles.data.T = DATA.data(:,T_col);
        handles.data.DeltaEoE = DATA.data(:,DeltaEoE_col);
        
        % calculate Delta_n
        lambda = handles.params.lambda * 1e-9;
        l = handles.params.l * 1e-3;
        w_0 = handles.params.w_0 * 1e-6;
        w_0_p = handles.params.w_0p * 1e-6;
        k_0 = 2 * pi / lambda;
        handles.data.Delta_n = handles.data.DeltaEoE / (k_0*l*sqrt(2/exp(1))*(w_0 / w_0_p)*(2 / sqrt(pi)));

        clearvars DATA % This is no longer needed, so we can delete it to save space

        % set flags to indicate valid data has been imported but not fit
        handles.params.data_imported = 1;
        handles.params.data_fit = 0;
        
        set(handles.pb_Export,'Enable','off'); % disable the export button until a new fit has been performed
        
        % refresh the plots
        refresh_plots(hObject,handles);
                
        % enable the fit button if the inputs are valid
        if inputs_valid(handles.params)
            
           set(handles.pb_Fit,'Enable','on');

        else
            
           set(handles.pb_Fit,'Enable','off');
            
        end
        
        % update the file name on the display
        temp = strsplit(handles.params.pathname,'\');
        dir_cell = cellstr(temp(end - 1));
        dir_str = dir_cell{1};
        file_display_name = ['\',dir_str,'\',handles.params.filename];
        set(handles.st_filename,'String',file_display_name);
        disp([datestr(now),': BD_Fresnel - Import of ',handles.params.pathname,handles.params.filename,' complete.'])
        
        save_settings(hObject,handles); % save settings
                     
    else
        
        % show error if file is invalid
        msgbox(['File: ',handles.params.filename,' is invalid. Please select a valid data file.'],'File Import Error');
        
    end

else
    
    disp([datestr(now),': BD_Fresnel - File import cancelled.'])
    
end
guidata(hObject, handles); % Update the handles.
    

% --- Executes on button press in pb_Fit.
function pb_Fit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic;
disp([datestr(now),': BD_Fresnel - Calulating fit.'])

% Disable the buttons
set(handles.pb_Fit,'String','Processing...');
set(handles.pb_Fit,'Enable','off');
set(handles.pb_Import,'Enable','off');
set(handles.pb_Export,'Enable','off');
guidata(hObject, handles);
drawnow

% calculate fitting parameters
handles.params.T_d_min = min(handles.data.tau_d)/handles.params.tau_p;
handles.params.T_d_max = max(handles.data.tau_d)/handles.params.tau_p;
handles.params.T_d_num = 150; % number of points in T and T_d

guidata(hObject, handles);

% run the fitting program
handles.fit = BD_Fresnel_func_ni(handles.params);
handles.params.data_fit = 1;

% refresh the plots
refresh_plots(hObject,handles);

% Enable the buttons
set(handles.pb_Fit,'String','Fit');
set(handles.pb_Fit,'Enable','on');
set(handles.pb_Import,'Enable','on');
set(handles.pb_Export,'Enable','on');

elapsedTime = toc;
disp([datestr(now),': BD_Fresnel - Fitting complete. Elapsed time: ',num2str(elapsedTime),' s.'])

save_settings(hObject,handles); % save settings

guidata(hObject, handles);


% --- Executes on button press in pb_Export.
function pb_Export_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp([datestr(now),': BD_Fresnel - Begin fit export.'])

% get the export filename
splitInd = regexp(handles.params.filename,'.txt');
exportFileName = strcat(handles.params.filename(1:splitInd-1),'_fit.txt');
export_filepath = [handles.params.pathname,exportFileName];
[exportFileName,exportPath] = uiputfile(export_filepath,'Select filename for fit export');

total_spacing = 20;

if ~strcmp(num2str(exportFileName),'0') % just in case the user cancels the operation.

    export_filepath = [exportPath,exportFileName]; % the combined export filepath

    % interpolate the fit to the data grid
    T_interp = interp1(handles.fit.tau_d,handles.fit.T,handles.data.tau_d_shifted);
    DeltaEoE_interp = interp1(handles.fit.tau_d,handles.fit.DeltaEoE,handles.data.tau_d_shifted);
    Delta_n_interp = interp1(handles.fit.tau_d,handles.fit.Delta_n,handles.data.tau_d_shifted);

    % write the export data to the file
    export_data = cat(2,handles.data.tau_d_shifted,handles.data.DeltaEoE,handles.data.T,handles.data.Delta_n,...
        DeltaEoE_interp,T_interp,Delta_n_interp);
    
    export_fid = fopen(export_filepath,'w+');
    fprintf(export_fid,'Delay (fs), DeltaE / E (Data), T (Data), Delta_n (data), DeltaE / E (Fit), T (Fit), Delta_n (fit)');
    dlmwrite(export_filepath,export_data,'newline','pc','delimiter',',','precision',8,'-append','roffset',1);
    disp([datestr(now),': BD_Fresnel - Export of ',exportFileName,' complete.'])

    % Write Log
    logfile_exportname = strrep(exportFileName,'.txt','_log.txt');
    logfile_exportpath = [exportPath,logfile_exportname];
    logfile_fid = fopen(logfile_exportpath,'w+');
    fprintf(logfile_fid,'%s\r\n\r\n',[datestr(now),': Begin beam deflection program: ',mfilename('fullpath'),'.']);
    fprintf(logfile_fid,'%s\r\n\r\n',['Import data file ',handles.params.import_filepath,'.']);

    % Write Data Offsets to Log
    fprintf(logfile_fid,'%s\r\n','Data Offsets:');
    fprintf(logfile_fid,'%s\r\n',['    shift:               ',num2str(handles.params.shift),' fs']);
    fprintf(logfile_fid,'%s\r\n','');

    % Write Simulation Parameters to Log
    fprintf(logfile_fid,'%s\r\n','Simulation Parameters:');
    fprintf(logfile_fid,'%s\r\n',['    X_num:               ',num2str(handles.params.X_num),]);
    fprintf(logfile_fid,'%s\r\n',['    X_max:               ',num2str(handles.params.X_max),' w_0']); 
    fprintf(logfile_fid,'%s\r\n',['    T_num:               ',num2str(handles.params.T_num)]);
    fprintf(logfile_fid,'%s\r\n',['    T_max:               ',num2str(handles.params.T_max),' tau_p']);
    fprintf(logfile_fid,'%s\r\n','');

    % Write Fit Parameters to Log
    fprintf(logfile_fid,'%s\r\n','Excitation parameters:');
    fprintf(logfile_fid,'%s\r\n',['    E_e:                 ',num2str(handles.params.E_p),' uJ']);  
    fprintf(logfile_fid,'%s\r\n',['    lambda_e:            ',num2str(handles.params.lambda_p),' nm']);
    fprintf(logfile_fid,'%s\r\n',['    tau_e (FWHM):        ',num2str(handles.params.tau_p_FWHM'),' fs']);    
    fprintf(logfile_fid,'%s\r\n',['    tau_e (HW1/eM):      ',num2str(handles.params.tau_p,'%6.3f'),' fs']);
    fprintf(logfile_fid,'%s\r\n',['    w_0e:                ',num2str(handles.params.w_0p),' um']);
    fprintf(logfile_fid,'%s\r\n',['    I_0e:                ',num2str(handles.params.I_0_p*1e-13,'%6.3f'),' GW/cm^2']);
    fprintf(logfile_fid,'%s\r\n','');

    fprintf(logfile_fid,'%s\r\n','Probe parameters:');
    fprintf(logfile_fid,'%s\r\n',['    lambda_p:            ',num2str(handles.params.lambda),' nm']);
    fprintf(logfile_fid,'%s\r\n',['    tau_p (FWHM):        ',num2str(handles.params.tau_FWHM),' fs']);
    fprintf(logfile_fid,'%s\r\n',['    tau_p (HW1/eM):      ',num2str(handles.params.tau,'%6.3f'),' fs']);
    fprintf(logfile_fid,'%s\r\n',['    w_0p:                ',num2str(handles.params.w_0),' um']);
    fprintf(logfile_fid,'%s\r\n','');

    fprintf(logfile_fid,'%s\r\n','System parameters:');
    fprintf(logfile_fid,'%s\r\n',['    d:                   ',num2str(handles.params.d),' cm']);
    fprintf(logfile_fid,'%s\r\n',['    x_0:                 ',num2str(handles.params.x_0),' w_0_p']);
    fprintf(logfile_fid,'%s\r\n','');

    fprintf(logfile_fid,'%s\r\n','Material parameters:');
    fprintf(logfile_fid,'%s\r\n',['    L:                   ',num2str(handles.params.l),' mm']);
    fprintf(logfile_fid,'%s\r\n',['    rho:                 ',num2str(handles.params.rho,'%6.3f')]);
    fprintf(logfile_fid,'%s\r\n','');

    fprintf(logfile_fid,'%s\r\n','  At excitation wavelength:');
    fprintf(logfile_fid,'%s\r\n',['    R_e:                 ',num2str(handles.params.R_p*100,'%6.3f'),' %']);
    fprintf(logfile_fid,'%s\r\n',['    alpha_e:             ',num2str(handles.params.alpha_p,'%6.3f'),' cm^-1']);
    fprintf(logfile_fid,'%s\r\n',['    n_e:                 ',num2str(handles.params.n_p)]);
    fprintf(logfile_fid,'%s\r\n',['    dndl_e:              ',num2str(handles.params.dndl_p,'%6.6f'),' um^-1']);
    fprintf(logfile_fid,'%s\r\n',['    ng_e:                ',num2str(handles.params.ng_p,'%6.3f')]);
    fprintf(logfile_fid,'%s\r\n',['    Tlin_e:              ',num2str(handles.params.Tlin_p*100,'%6.3f'),' %']);
    fprintf(logfile_fid,'%s\r\n','');
    
    fprintf(logfile_fid,'%s\r\n','  At probe wavelength:');
    fprintf(logfile_fid,'%s\r\n',['    R_p:                 ',num2str(handles.params.R*100,'%6.3f'),' %']);
    fprintf(logfile_fid,'%s\r\n',['    alpha_p:             ',num2str(handles.params.alpha,'%6.3f'),' cm^-1']);
    fprintf(logfile_fid,'%s\r\n',['    n_p:                 ',num2str(handles.params.n)]);
    fprintf(logfile_fid,'%s\r\n',['    dndl_p:              ',num2str(handles.params.dndl,'%6.6f'),' um^-1']);
    fprintf(logfile_fid,'%s\r\n',['    ng_p:                ',num2str(handles.params.ng,'%6.3f')]);
    fprintf(logfile_fid,'%s\r\n',['    Tlin_p:              ',num2str(handles.params.Tlin*100,'%6.3f'),' %']);
    fprintf(logfile_fid,'%s\r\n','');   

    for i = 1:length(handles.params.mech) % write parameters for each mechanism
        
        mech = handles.params.mech(i);
        if ~strcmp(mech.name,'None (both)') % only if the mechanism is used
        
            fprintf(logfile_fid,'%s\r\n',['Mechanism ',num2str(i),':']);
            fprintf(logfile_fid,'%s\r\n',['    Type:                ',cell2mat(mech.name)]);

            param = mech.param;
            for j = 1:length(param) % for each parameter within each mechanism

                parameter_name = param(j).name;
                if ~isempty(parameter_name) % for the parameters used
                
                    spacing = blanks(total_spacing - length(parameter_name)); % calculate the number of spaces required
                    fprintf(logfile_fid,'%s\r\n',['    ',parameter_name,':',spacing,num2str(param(j).value),' ',param(j).units]);

                end
                                
            end

            fprintf(logfile_fid,'%s\r\n','');
        
        end
       
    end
    
    fprintf(logfile_fid,'%s\r\n\r\n',['Exported fitting data to ',export_filepath,'.']);
    
    fprintf(logfile_fid,'%s\r\n\r\n','Please cite the following references when using products from this program.');
    fprintf(logfile_fid,'%s\r\n\r\n','    1: Reichert, M., et al. (2014). "Temporal, spectral, and polarization dependence of the nonlinear optical response of carbon disulfide." Optica 1(6): 436.');
    fprintf(logfile_fid,'%s\r\n\r\n','    2: Reichert, M., et al. (2014). "Temporal, spectral, and polarization dependence of the nonlinear optical response of carbon disulfide: supplementary material." Optica 1(6): 436.');
    fprintf(logfile_fid,'%s\r\n\r\n','    3: Ferdinandus, M. R., et al. (2013). "Beam deflection measurement of time and polarization resolved ultrafast nonlinear refraction." Opt Lett 38(18): 3518-3521.');
    
    fclose(logfile_fid); % Close log file.
    
    disp([datestr(now),': BD_Fresnel - Export of ',logfile_exportname,' complete.'])

else
    
    disp([datestr(now),': BD_Fresnel - Fit export cancelled.'])

end
    
guidata(hObject,handles); % Update the handles


% --- Executes on button press in pb_Exit.
function pb_Exit_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

save_settings(hObject,handles); % save settings
disp([datestr(now),': BD_Fresnel - Saving settings to ',handles.params.save_filepath,'.']);

disp([datestr(now),': BD_Fresnel - End BD_Fresnel_GUI.'])
delete(handles.zfit_manualfit)


% --- Executes on button press in cb_diagnostics.
function cb_diagnostics_Callback(hObject, eventdata, handles)
% hObject    handle to cb_diagnostics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_diagnostics

handles.params.diagnostics = get(handles.cb_diagnostics,'Value');

guidata(hObject,handles); % Update the handles


% --- Check to see if the inputs are valid
function output = inputs_valid(params)

    if params.E_p > 0 &&....
        params.tau_p_FWHM > 0 &&....
        params.w_0p > 0 &&....
        params.lambda_p > 0 &&....
        params.tau_FWHM > 0 &&....
        params.w_0 > 0 &&....
        params.lambda > 0 &&....
        params.n_p >= 1 &&....
        params.n >= 1 &&....
        params.R >= 0 && params.R < 1 &&....
        params.R_p >= 0 && params.R_p < 1 &&....
        params.alpha >= 0 &&....
        params.alpha_p >= 0 &&....
        params.d > 0 &&....
        params.l > 0 &&....
        params.data_imported == 1
    
        output = 1;
            
    else
        
        output = 0;
        
    end


% --- Refresh the plots.
function handles = refresh_plots(hObject,handles)

% set the data marker size
n_tau_d = length(handles.data.tau_d);
if n_tau_d > 100

    handles.params.dataMarkerSize = ceil(6*100/n_tau_d);

else
    
    handles.params.dataMarkerSize = 6;
    
end
    
% determine the timescale
Delta_tau_d =  abs(max(handles.data.tau_d) - min(handles.data.tau_d)); % timescale

if Delta_tau_d > 1000^2 % if on the nanosecond timescale

    handles.data.tau_d_shifted_display = handles.data.tau_d_shifted/(1000^2);
    handles.params.tau_d_scale = 'ns';

elseif Delta_tau_d > 1000 % if on the picosecond timescale

    handles.data.tau_d_shifted_display = handles.data.tau_d_shifted/1000;
    handles.params.tau_d_scale = 'ps';

else % if on the femtosecond timescale

    handles.data.tau_d_shifted_display = handles.data.tau_d_shifted;
    handles.params.tau_d_scale = 'fs';

end

% plot a bit different depending on whether there is an active fit
if handles.params.data_fit % if fit has been performed
    
    if Delta_tau_d > 1000^2 % if on the nanosecond timescale

        handles.fit.tau_d_display = handles.fit.tau_d/(1000^2);

    elseif Delta_tau_d > 1000 % if on the picosecond timescale

        handles.fit.tau_d_display = handles.fit.tau_d/1000;

    else % if on the femtosecond timescale

        handles.fit.tau_d_display = handles.fit.tau_d;

    end
    
    % combine the points for the data and the fit
    Tpoints = [handles.data.T;handles.fit.T];
    DeltaEoEpoints = [handles.data.DeltaEoE;handles.fit.DeltaEoE];
    
    % draw the transmission plot
    plot(handles.T_plot,handles.fit.tau_d_display,100*handles.fit.T,'-r','MarkerSize',handles.params.dataMarkerSize); % plot the fit
    hold(handles.T_plot,'on');

    % draw the refraction plot
    plot(handles.DeltaEoE_plot,handles.fit.tau_d_display,100*handles.fit.DeltaEoE,'-r','MarkerSize',handles.params.dataMarkerSize); % plot the fit
    hold(handles.DeltaEoE_plot,'on');
    
else % if no fit has been performed
   
    Tpoints = handles.data.T;
    DeltaEoEpoints = handles.data.DeltaEoE;
    
end

% determine mins and maxes of x and y axes
xMin = min(handles.data.tau_d_shifted_display);
xMax = max(handles.data.tau_d_shifted_display);

TyMin = min(Tpoints);
TyMax = max(Tpoints);
DeltaEoEyMin = min(DeltaEoEpoints);
DeltaEoEyMax = max(DeltaEoEpoints);

% plot the transmission data
plot(handles.T_plot,handles.data.tau_d_shifted_display,100*handles.data.T,'ob','MarkerSize',handles.params.dataMarkerSize);
grid(handles.T_plot,'on');
hold(handles.T_plot,'off');

% plot the deflection data
plot(handles.DeltaEoE_plot,handles.data.tau_d_shifted_display,100*handles.data.DeltaEoE,'ob','MarkerSize',handles.params.dataMarkerSize); % plot the fit
grid(handles.DeltaEoE_plot,'on');
hold(handles.DeltaEoE_plot,'off');

% set limits and labels for the transmission plot
xlim(handles.T_plot,[xMin,xMax]);
ylim(handles.T_plot,100*[TyMin,TyMax]);
xlabel(handles.T_plot,['tau_d (',handles.params.tau_d_scale,')'],'FontSize',handles.params.axesFontSize);
ylabel(handles.T_plot,'T (%)','FontSize',handles.params.axesFontSize);

% set limits and labels for the deflection plot
xlim(handles.DeltaEoE_plot,[xMin,xMax]);
ylim(handles.DeltaEoE_plot,100*[DeltaEoEyMin,DeltaEoEyMax]);
xlabel(handles.DeltaEoE_plot,['tau_d (',handles.params.tau_d_scale,')'],'FontSize',handles.params.axesFontSize);
ylabel(handles.DeltaEoE_plot,'DeltaEoE (%)','FontSize',handles.params.axesFontSize);
set(handles.DeltaEoE_plot,'YColor','Black'); % color for the right axis

guidata(hObject,handles); % Update the handles


% --- Refresh the calculated beam parameters
function handles = calc_beam_params(hObject,handles)

c = 299797548; % speed of light (m/s)

% pump parameters
E_p = handles.params.E_p*(1e-6); % pump input energy [J]
lambda_p = handles.params.lambda_p*(1e-9); % pump wavelength [m]
tau_p = handles.params.tau_p*(1e-15); % pump temporal pulse width (FWHM) [s]
w_0_p = handles.params.w_0p*(1e-6); % pump spot size (HW1/e^2 of I) [m]

% probe parameters
lambda = handles.params.lambda*(1e-9); % probe wavelength [m]
w_0 = handles.params.w_0*(1e-6); % probe spot size (HW1/e^2 or I) [m]

L = handles.params.l*(1e-3); % thickness [m]

% material parameters at probe wavelength
alpha = handles.params.alpha*100; % absorption coefficient of probe
n = handles.params.n; % index at probe - can be found on www.refractiveindex.org
dn_dl = handles.params.dndl*(1e+6); % dispersion at Probe [m^-1]

% material parameters at pump wavelength
R_p = handles.params.R_p; % reflectivity of pump
alpha_p = handles.params.alpha_p*100; % absorption coefficient of pump
n_p = handles.params.n_p; % index at pump wavelength
dn_dlp = handles.params.dndl_p*(1e+6); % dispersion at Pump [m^-1]
alpha_2d = handles.params.alpha_2d*(1e-11); % effective degenerate two photon absorption of pump

% calculated beam parameters
I_0_p = 2 * E_p/(pi^(3/2)*w_0_p^2*tau_p)*(1 - R_p); % peak irradiance of pump [W/m^2]
handles.params.I_0_p = I_0_p;

handles.params.sigma_p = alpha_p*L/2; % linear absorption parameter of pump
sigma_p = handles.params.sigma_p;

sigma_p2 = I_0_p*alpha_2d*L/2; % effective d-2PA coefficient of pump

% group velocity calculation
handles.params.ng_p = n_p - lambda_p*dn_dlp; % group index at pump
ng_p = handles.params.ng_p;

handles.params.ng = n - lambda*dn_dl; % group index at probe
ng = handles.params.ng;

% Rayleigh ranges
handles.params.z_0 = pi*w_0^2/lambda; % rayleigh Range [m];
handles.params.z_0_p = pi*w_0_p^2/lambda_p; % rayleigh Range [m]

% transmissions
handles.params.Tlin = exp(-alpha*L); % transmission of probe

if sigma_p2 == 0 % transmission of pump
    
    handles.params.Tlin_p = (1 - R_p)*exp(-2*sigma_p);
    
elseif sigma_p == 0 && sigma_p2 ~= 0
    
    handles.params.Tlin_p = (1 - R_p)*(-PolyLog(3/2,-2*sigma_p2)/(2*sigma_p2));
    
else
    
    A1 = -1 + exp(2*sigma_p);
    
    handles.params.Tlin_p = (1 - R_p)*(-sigma_p*PolyLog(3/2,-(A1*exp(-2*sigma_p)*sigma_p2)/sigma_p)/...
        (A1*sigma_p2));

end
    
handles.params.rho = L/(tau_p*c)*(ng - ng_p); % GVM Parameter

% calculated probe parameters
set(handles.st_ng_val,'String',num2str(handles.params.ng,'%6.3f'));
set(handles.st_Tlin_val,'String',num2str(handles.params.Tlin,'%6.3f'));
set(handles.st_z_0_val,'String',[num2str(handles.params.z_0*1e3,'%6.3f'),' mm']);

% calculated pump parameters
set(handles.st_ng_p_val,'String',num2str(handles.params.ng_p,'%6.3f'));
set(handles.st_Tlin_p_val,'String',num2str(handles.params.Tlin_p,'%6.3f'));
set(handles.st_z_0_p_val,'String',[num2str(handles.params.z_0_p*1e3,'%6.3f'),' mm']);
set(handles.st_I_0_p_val,'String',[num2str(handles.params.I_0_p*1e-13,'%6.3f'),' GW/cm^2']);

% other parameters
set(handles.st_rho_val,'String',num2str(handles.params.rho,'%6.3f'));

guidata(hObject,handles); % Update the handles


% --- Save settings.
function handles = save_settings(hObject,handles)

settings = handles.params;
save(handles.params.save_filepath,'settings'); % save the file path

guidata(hObject,handles); % Update the handles
