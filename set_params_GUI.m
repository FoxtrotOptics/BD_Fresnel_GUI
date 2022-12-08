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

function varargout = set_params_GUI(varargin)
% SET_PARAMS_GUI MATLAB code for set_params_GUI.fig
%      SET_PARAMS_GUI, by itself, creates a new SET_PARAMS_GUI or raises the existing
%      singleton*.
%
%      H = SET_PARAMS_GUI returns the handle to a new SET_PARAMS_GUI or the handle to
%      the existing singleton*.
%
%      SET_PARAMS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SET_PARAMS_GUI.M with the given input arguments.
%
%      SET_PARAMS_GUI('Property','Value',...) creates a new SET_PARAMS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before set_params_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to set_params_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help set_params_GUI

% Last Modified by GUIDE v2.5 07-Mar-2016 15:04:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @set_params_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @set_params_GUI_OutputFcn, ...
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


% --- Executes just before set_params_GUI is made visible.
function set_params_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to set_params_GUI (see VARARGIN)

% Choose default command line output for set_params_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.input = varargin{1}; % parse out the inputs
handles.output = handles.input;
mechanism = handles.input.mechanism;

disp([datestr(now),': set_params_GUI - Input parameters for mechanism: ',cell2mat(handles.input.name),'.'])

n_params = length(handles.input.param);

% restore previous settings if required
mechanism_names = varargin(2);
mechanism_names = mechanism_names{1};

if strcmp(handles.input.name,cell2mat(mechanism_names(mechanism)))

    handles.output.param = handles.input.param;

else

    for i = 1:n_params % set the values to zero for new mech

        handles.output.param(i).value = 0;

    end

end

set(handles.st_mechanism_parameters,'String',[cell2mat(mechanism_names(mechanism)),' parameters']); % set title

% initialize values
param_names = {'','','','','','','','','',''};
param_units = {'','','','','','','','','',''};

% set up the interface
switch mechanism
    
    case 1
        
        % delete the parameters from the structure except the mechanism
        handles.output.name = mechanism_names(1);
        param_names = {'','','','','','','','','',''};
        param_units = {'','','','','','','','','',''};
        
    case 2
               
        % set up the parameter fields
        handles.output.name = mechanism_names(2);
        param_names = {'alpha_2e','n_2e','','','','','','','',''};
        param_units = {'e-11 m/W','e-19 m^2/W','','','','','','','',''};
       
    case 3
    
        % set up the parameter fields
        handles.output.name = mechanism_names(3);
        param_names = {'n_2c','tau_cr','tau_cf','','','','','','',''};
        param_units = {'e-19 m^2/W','fs','fs','','','','','','',''};
        
    case 4

        % set up the parameter fields
        handles.output.name = mechanism_names(4);
        param_names = {'n_2l','omega_0','sigma','tau_lf','','','','','',''};
        param_units = {'e-19 m^2/W','ps^-1','ps^-1','fs','','','','','',''};
      
    case 5
        
        % set up the parameter fields
        handles.output.name = mechanism_names(5);
        param_names = {'n_2r','tau_rr','tau_rf','','','','','','',''};
        param_units = {'e-19 m^2/W','fs','fs','','','','','','',''};
        
    case 6
        
        % set up the parameter fields
        handles.output.name = mechanism_names(6);
        param_names = {'n_2v','w_0v','tau_vf','','','','','','',''};
        param_units = {'e-19 m^2/W','ps^-1','fs','','','','','','',''};
        
    case 7
        
        % set up the parameter fields
        handles.output.name = mechanism_names(7);
        param_names = {'alpha_2','N0','sigma_S_01','sigma_S_12','sigma_T_12','tau_S_10','tau_S_21','tau_ISC','tau_T_21','tau_TS'};
        param_units = {'e-11 m/W','e+22 m^-3','e-21 m^2','e-21 m^2','e-21 m^2','ps','ps','ps','ps','ps'};
        
    case 8
        
        % set up the parameter fields
        handles.output.name = mechanism_names(8);
        param_names = {'alpha_2_exp','tau_expr','tau_expf','','','','','','',''};
        param_units = {'e-11 m/W','fs','fs','','','','','','',''};

    case 9
        
        % set up the parameter fields
        handles.output.name = mechanism_names(9);
        param_names = {'n_2_exp','tau_expr','tau_expf','','','','','','',''};
        param_units = {'e-19 m^2/W','fs','fs','','','','','','',''};
        
    case 10
        
        % set up the parameter fields
        handles.output.name = mechanism_names(10);
        param_names = {'alpha_2_exp','n_2_exp','tau_expr','tau_expf','','','','','',''};
        param_units = {'e-11 m/W','e-19 m^2/W','fs','fs','','','','','',''};
        
end

% find the inputs that are unused
inactive_cells = cellfun(@isempty,param_names);
active_cells = ~inactive_cells;

inactive_start = find(inactive_cells,1,'first');
active_end = find(active_cells,1,'last');

if ~isempty(active_end)

    for i = 1:active_end % for all the input boxes

        % load the parameter names and units
        handles.output.param(i).name = cell2mat(param_names(i));
        handles.output.param(i).units = cell2mat(param_units(i));

        % set up the input boxes
        set(eval(['handles.st_param_',num2str(i),'_name']),'Visible','On');
        set(eval(['handles.st_param_',num2str(i),'_name']),'String',[handles.output.param(i).name,' (',handles.output.param(i).units,')']);
        set(eval(['handles.et_param_',num2str(i),'_value']),'Visible','On');
        set(eval(['handles.et_param_',num2str(i),'_value']),'String',num2str(handles.output.param(i).value));

    end

end

if ~isempty(inactive_start) % if there are empty input cells to be blanked

    for j = inactive_start:n_params % blank th unused inputs

        % blank the unused input boxes
        set(eval(['handles.st_param_',num2str(j),'_name']),'Visible','Off');
        set(eval(['handles.et_param_',num2str(j),'_value']),'Visible','Off');

    end

end

guidata(hObject,handles); % Update the handles

% UIWAIT makes set_params_GUI wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = set_params_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1)



function et_param_1_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_1_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_1_value as a double

temp = str2double(get(handles.et_param_1_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(1).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(1).name,' updated to ',num2str(handles.output.param(1).value),' ',handles.output.param(1).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(1).name,' must be a real number.'])
    handles.output.param(1).name = 0; % Update the variable
    
end

set(handles.et_param_1_value,'String',num2str(handles.output.param(1).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_1_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_1_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_2_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_2_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_2_value as a double

temp = str2double(get(handles.et_param_2_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(2).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(2).name,' updated to ',num2str(handles.output.param(2).value),' ',handles.output.param(2).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(2).name,' must be a real number.'])
    handles.output.param(2).name = 0; % Update the variable
    
end

set(handles.et_param_2_value,'String',num2str(handles.output.param(2).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_2_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_2_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_3_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_3_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_3_value as a double

temp = str2double(get(handles.et_param_3_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(3).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(3).name,' updated to ',num2str(handles.output.param(3).value),' ',handles.output.param(3).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(3).name,' must be a real number.'])
    handles.output.param(3).name = 0; % Update the variable
    
end

set(handles.et_param_3_value,'String',num2str(handles.output.param(3).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_3_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_3_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_4_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_4_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_4_value as a double

temp = str2double(get(handles.et_param_4_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(4).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(4).name,' updated to ',num2str(handles.output.param(4).value),' ',handles.output.param(4).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(4).name,' must be a real number.'])
    handles.output.param(4).name = 0; % Update the variable
    
end

set(handles.et_param_4_value,'String',num2str(handles.output.param(4).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_4_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_4_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_5_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_5_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_5_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_5_value as a double

temp = str2double(get(handles.et_param_5_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(5).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(5).name,' updated to ',num2str(handles.output.param(5).value),' ',handles.output.param(5).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(5).name,' must be a real number.'])
    handles.output.param(5).name = 0; % Update the variable
    
end

set(handles.et_param_5_value,'String',num2str(handles.output.param(5).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_5_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_5_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_6_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_6_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_6_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_6_value as a double

temp = str2double(get(handles.et_param_6_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(6).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(6).name,' updated to ',num2str(handles.output.param(6).value),' ',handles.output.param(6).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(6).name,' must be a real number.'])
    handles.output.param(6).name = 0; % Update the variable
    
end

set(handles.et_param_6_value,'String',num2str(handles.output.param(6).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_6_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_6_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_7_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_7_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_7_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_7_value as a double

temp = str2double(get(handles.et_param_7_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(7).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(7).name,' updated to ',num2str(handles.output.param(7).value),' ',handles.output.param(7).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(7).name,' must be a real number.'])
    handles.output.param(7).name = 0; % Update the variable
    
end

set(handles.et_param_7_value,'String',num2str(handles.output.param(7).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_7_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_7_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_8_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_8_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_8_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_8_value as a double

temp = str2double(get(handles.et_param_8_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(8).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(8).name,' updated to ',num2str(handles.output.param(8).value),' ',handles.output.param(8).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(8).name,' must be a real number.'])
    handles.output.param(8).name = 0; % Update the variable
    
end

set(handles.et_param_8_value,'String',num2str(handles.output.param(8).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_8_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_8_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_9_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_9_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_9_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_9_value as a double

temp = str2double(get(handles.et_param_9_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(9).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(9).name,' updated to ',num2str(handles.output.param(9).value),' ',handles.output.param(9).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(9).name,' must be a real number.'])
    handles.output.param(9).name = 0; % Update the variable
    
end

set(handles.et_param_9_value,'String',num2str(handles.output.param(9).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_9_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_9_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_param_10_value_Callback(hObject, eventdata, handles)
% hObject    handle to et_param_10_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_param_10_value as text
%        str2double(get(hObject,'String')) returns contents of et_param_10_value as a double

temp = str2double(get(handles.et_param_10_value,'String')); % get the string in the box and convert to double
if ~isnan(temp) % if the value is a number

    handles.output.param(10).value = temp; % Update the variable
    disp([datestr(now),': set_params_GUI - Value of ',handles.output.param(10).name,' updated to ',num2str(handles.output.param(10).value),' ',handles.output.param(10).units,'.'])

else
    
    disp([datestr(now),': set_params_GUI - ',handles.output.param(10).name,' must be a real number.'])
    handles.output.param(10).name = 0; % Update the variable
    
end

set(handles.et_param_10_value,'String',num2str(handles.output.param(10).value)) % update the display
guidata(hObject,handles); % Update the handles


% --- Executes during object creation, after setting all properties.
function et_param_10_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_param_10_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp([datestr(now),': set_params_GUI - Parameters updated.'])

uiresume(handles.figure1); % resume the UI and return value
