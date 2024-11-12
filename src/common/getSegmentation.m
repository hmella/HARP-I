function varargout = getSegmentation(varargin)
% getSegmentation MATLAB code for getSegmentation.fig
%      getSegmentation, by itself, creates a new getSegmentation or raises the existing
%      singleton*.
%
%      H = getSegmentation returns the handle to a new getSegmentation or the handle to
%      the existing singleton*.
%
%      getSegmentation('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in getSegmentation.M with the given input arguments.
%
%      getSegmentation('Property','Value',...) creates a new getSegmentation or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getSegmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getSegmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools meinu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getSegmentation

% Last Modified by GUIDE v2.5 26-Jul-2020 21:47:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getSegmentation_OpeningFcn, ...
                   'gui_OutputFcn',  @getSegmentation_OutputFcn, ...
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


% --- Executes just before getSegmentation is made visible.
function getSegmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getSegmentation (see VARARGIN)

  %% Create the data to plot.
  % Mask, phase and contours
  handles.I = varargin{1}.Image;
  handles.Pha = varargin{1}.Phase;

  % Image size and number of frames
  handles.Isz = size(handles.I,[1 2]);
  handles.Nfr = size(handles.I,3);

  % Segmentation
  handles.mask = false([handles.Isz handles.Nfr]);
  
  % Graph-cuts positions
  handles.positions = {};
  handles.curves = {};
  handles.contours = cell([1 handles.Nfr]);
  handles.contours_number = 2;
  handles.contours_positions = cell([1 handles.Nfr]);
  handles.contours_isclosed = cell([1 handles.Nfr]);
  handles.contours_iscorner = cell([1 handles.Nfr]);
  handles.contours_iscurved = cell([1 handles.Nfr]);
  try 
      handles.contours = {varargin{1}.Contours};
      for i=1:handles.Nfr
          for j=1:handles.contours_number
              handles.contours_positions{i}{j} = handles.contours{i}.Position{j};
              handles.contours_isclosed{i}{j} = handles.contours{i}.IsClosed{j};
              handles.contours_iscorner{i}{j} = handles.contours{i}.IsCorner{j};
              handles.contours_iscurved{i}{j} = handles.contours{i}.IsCurved{j};              
          end
      end
  catch
      fprintf('\n No contours were found on inputs arguments')
  end
  
  % Initial frame
  handles.frame = 1;
  
  % Initial axis
  try
      handles.axis = varargin{1}.Axis;
  catch
      handles.axis = [1 handles.Isz(1) 1 handles.Isz(2)];
  end
  handles.edit2.String = num2str(handles.axis(1));
  handles.edit3.String = num2str(handles.axis(2));
  handles.edit4.String = num2str(handles.axis(3));
  handles.edit5.String = num2str(handles.axis(4));  
  
  % Initial color axis
  handles.caxis = [0 0.5];

  % Image positions
  [X,Y] = meshgrid(0:(handles.Isz(2)-1),0:(handles.Isz(1)-1));
  handles.X = X;
  handles.Y = Y;

  % Update frame indicator
  set(handles.text4,'String',sprintf('Processing frame %.0d',handles.frame));  
  
  % Plot wrapped phase
  axes(handles.axes1)
  h=imagesc(handles.axes1,handles.X(:),handles.Y(:),handles.I(:,:,handles.frame));
  colormap gray; caxis(handles.axes1,handles.caxis);
  set(handles.axes1,'YDir','Normal');
  set(handles.axes1,'visible','off');
  axis(handles.axes1,handles.axis)
%   set(gcf,'toolbar','figure');  
  
  % Store image handle for masks creation
  handles.image_handle = h;

  % Plot Phases
  imagesc(handles.axes2,handles.X(:),handles.Y(:),handles.Pha(:,:,1,handles.frame));
  set(handles.axes2,'YDir','Normal');
  set(handles.axes2,'visible','off');
  axis(handles.axes2,handles.axis);

  imagesc(handles.axes3,handles.X(:),handles.Y(:),handles.Pha(:,:,2,handles.frame));
  set(handles.axes3,'YDir','Normal');
  set(handles.axes3,'visible','off');
  axis(handles.axes3,handles.axis);  
  
  % Plot the contours
  if ~isempty(handles.contours{1})
      % Edit contour lines
      for H = [handles.axes1,handles.axes2,handles.axes3]
          h = imcline(handles.contours{handles.frame},H);
          h.Enable = 'on';
          h.Visible = 'on';
          h.IndependentDrag{2} = 'on';
          [h.Appearance(1:2).Color] = deal('r','g');
          [h.Appearance(1:2).MarkerFaceColor] = deal('r','g');
          [h.Appearance(1:2).MarkerSize] = deal(10, 10);
          [h.Appearance(1:2).Marker] = deal('o', 'o');
          [h.Appearance(1:2).LineWidth] = deal(2, 2);
      end
      linkaxes([handles.axes1,handles.axes2,handles.axes3]);
      iptPointerManager(handles.figure1,'enable')
  end
  
  % Choose default command line output for contours
  handles.output = hObject;

  % Update handles structure
  guidata(hObject, handles);

  % UIWAIT makes getSegmentation wait for user response (see UIRESUME)
%   uiwait(handles.figure1);
  uiwait();


% --- Outputs from this function are returned to the command line.
function varargout = getSegmentation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Outputs for save dialog
mask = handles.mask;
contours = handles.contours;

% Store segmentation
[filename, pathname] = uiputfile('*.mat','Export segmentation As');
try
  save(fullfile(pathname, filename),'mask','contours')
catch
  fprintf('\n The segmentation was not saved!\n')
end

% Hint: delete(hObject) closes the ficlosegure
% delete(hObject);
handles.output = handles;
guidata(hObject, handles);  % Store the outputs in the GUI
uiresume()                  % resume UI which will trigger the OutputFcn


% --- Executes on button press in Unwrap phase
function unwrap_phase_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Unwrap phases
XPha = unwrap2(handles.Pha(:,:,1,handles.frame),'Mask',handles.mask(:,:,handles.frame),'Connectivity',4,'Seed','auto');
YPha = unwrap2(handles.Pha(:,:,2,handles.frame),'Mask',handles.mask(:,:,handles.frame),'Connectivity',4,'Seed','auto');

% Plot images
imagesc(handles.axes2,handles.X(:),handles.Y(:),XPha,'AlphaData',handles.mask(:,:,handles.frame));
plot(handles.contours{handles.frame},'LineWidth',2,'Color','r','Parent',handles.axes2);
set(handles.axes2,'YDir','Normal')
set(handles.axes2,'visible','off');
axis(handles.axes2,handles.axis)

imagesc(handles.axes3,handles.X(:),handles.Y(:),YPha,'AlphaData',handles.mask(:,:,handles.frame));
plot(handles.contours{handles.frame},'LineWidth',2,'Color','r','Parent',handles.axes3);
set(handles.axes3,'YDir','Normal')
set(handles.axes3,'visible','off');
axis(handles.axes3,handles.axis)


% --- Executes on button press in pushbutton3.
function preview_mask_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Generate mask from current graph-cuts
args = struct(...
    'MaskSize',     handles.Isz,...
    'Contours',     handles.contours{handles.frame},...
    'Resolution',   0.5);
handles.mask(:,:,handles.frame) = getMask(args); 

% Update plots
imagesc(handles.axes2,handles.X(:),handles.Y(:),handles.Pha(:,:,1,handles.frame),'AlphaData',handles.mask(:,:,handles.frame));
set(handles.axes2,'YDir','Normal')
set(handles.axes2,'visible','off');
axis(handles.axes2,handles.axis)

imagesc(handles.axes3,handles.X(:),handles.Y(:),handles.Pha(:,:,2,handles.frame),'AlphaData',handles.mask(:,:,handles.frame));
set(handles.axes3,'YDir','Normal')
set(handles.axes3,'visible','off');
axis(handles.axes3,handles.axis)

% Update figure and handles
guidata(hObject,handles)


% --- Executes on button press in pushbutton10.
function generate_masks_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Generate mask from current graph-cuts
for frame=1:handles.Nfr
    if ~isempty(handles.contours{frame})
        args = struct(...
            'MaskSize',     handles.Isz,...
            'Contours',     handles.contours{frame},...
            'Resolution',   0.5);
        handles.mask(:,:,frame) = getMask(args);
    end
end
fprintf('\nMasks generated succesfully!\n')

% Update figure and handles
guidata(hObject,handles)


% --- Executes on button press in pushbutton4.
function next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Before updating the frame, the getSegmentation positions are stored and
% getSegmentation are resetted

% Store contours positions in the previous frame
if ~isempty(handles.contours{handles.frame})
    for i=1:handles.contours_number
        handles.contours_positions{handles.frame}{i} = handles.contours{handles.frame}.Position{i};
        handles.contours_isclosed{handles.frame}{i} = handles.contours{handles.frame}.IsClosed{i};
        handles.contours_iscorner{handles.frame}{i} = handles.contours{handles.frame}.IsCorner{i};
        handles.contours_iscurved{handles.frame}{i} = handles.contours{handles.frame}.IsCurved{i};
    end
end

% Update frame
if handles.frame == handles.Nfr
    handles.frame = 1;
else
    handles.frame = handles.frame + 1;
end

% Update frame indicator
set(handles.text4,'String',sprintf('Processing frame %.0d',handles.frame));

% Update plots
imagesc(handles.axes1,handles.X(:),handles.Y(:),handles.I(:,:,handles.frame));
colormap gray; caxis(handles.axes1,handles.caxis);
set(handles.axes1,'YDir','Normal');
set(handles.axes1,'visible','off');
axis(handles.axes1,handles.axis)

imagesc(handles.axes2,handles.X(:),handles.Y(:),handles.Pha(:,:,1,handles.frame));
set(handles.axes2,'YDir','Normal')
set(handles.axes2,'visible','off');
axis(handles.axes2,handles.axis)

imagesc(handles.axes3,handles.X(:),handles.Y(:),handles.Pha(:,:,2,handles.frame));
set(handles.axes3,'YDir','Normal')
set(handles.axes3,'visible','off');
axis(handles.axes3,handles.axis)  

% Check for previously stored contours. If there are not previous 
% contours the user can draw new ones
if ~isempty(handles.contours_positions{handles.frame})
    % Get contours
    handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame});
    handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame};
    handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame};
    handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame};
else
    % Retrieve and edit contours from previous frame
    if handles.frame ~= 1
        if and(~isempty(handles.contours_positions{handles.frame-1}), handles.frame ~= 1)
            % Get contours
            handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame-1});
            handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame-1};
            handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame-1};
            handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame-1};
        end
    end
    if handles.frame ~= handles.Nfr
        if ~isempty(handles.contours_positions{handles.frame+1})
            % Get contours
            handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame+1});
            handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame+1};
            handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame+1};
            handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame+1};
        end  
    end
end

% Show and edit contours
try
    % Edit contour lines
    for H = [handles.axes1,handles.axes2,handles.axes3]
        h = imcline(handles.contours{handles.frame},H);
        h.Enable = 'on';
        h.Visible = 'on';
        h.IndependentDrag{2} = 'on';
          [h.Appearance(1:2).Color] = deal('r','g');
          [h.Appearance(1:2).MarkerFaceColor] = deal('r','g');
          [h.Appearance(1:2).MarkerSize] = deal(10, 10);
          [h.Appearance(1:2).Marker] = deal('o', 'o');
          [h.Appearance(1:2).LineWidth] = deal(2, 2);
    end
    linkaxes([handles.axes1,handles.axes2,handles.axes3]);
    iptPointerManager(handles.figure1,'enable')
catch
end

% Update handles object
guidata(hObject, handles);


% --- Executes on button press in pushbutton5.
function previous_frame_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Store contours positions in the previouse frame
if ~isempty(handles.contours{handles.frame})
    fprintf('\n    Updating contours positions for frame %.0d',handles.frame)
    for i=1:handles.contours_number
        handles.contours_positions{handles.frame}{i} = handles.contours{handles.frame}.Position{i};
        handles.contours_isclosed{handles.frame}{i} = handles.contours{handles.frame}.IsClosed{i};
        handles.contours_iscorner{handles.frame}{i} = handles.contours{handles.frame}.IsCorner{i};
        handles.contours_iscurved{handles.frame}{i} = handles.contours{handles.frame}.IsCurved{i};
    end
end

% Update frame
if handles.frame == 1
    handles.frame = handles.Nfr;
else
    handles.frame = handles.frame - 1;
end

% Update frame indicator
set(handles.text4,'String',sprintf('Processing frame %.0d',handles.frame));

% update plots
imagesc(handles.axes1,handles.X(:),handles.Y(:),handles.I(:,:,handles.frame));
colormap gray; caxis(handles.axes1,handles.caxis);
set(handles.axes1,'YDir','Normal');
set(handles.axes1,'visible','off');
axis(handles.axes1,handles.axis)

imagesc(handles.axes2,handles.X(:),handles.Y(:),handles.Pha(:,:,1,handles.frame));
set(handles.axes2,'YDir','Normal')
set(handles.axes2,'visible','off');
axis(handles.axes2,handles.axis)

imagesc(handles.axes3,handles.X(:),handles.Y(:),handles.Pha(:,:,2,handles.frame));
set(handles.axes3,'YDir','Normal')
set(handles.axes3,'visible','off');
axis(handles.axes3,handles.axis)  

% Check for previously stored contours. If there are not previous 
% contours the user can draw new ones
if ~isempty(handles.contours_positions{handles.frame})
    % Get contours
    handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame});
    handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame};
    handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame};
    handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame};    
else
    % Retrieve and edit contours from previous frame
    if handles.frame ~= 1
        if ~isempty(handles.contours_positions{handles.frame-1})
            % Get contours
            handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame-1});
            handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame-1};
            handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame-1};
            handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame-1};            
        end
    end
    if handles.frame ~= handles.Nfr
        if ~isempty(handles.contours_positions{handles.frame+1})
            % Get contours
            handles.contours{handles.frame} = cline(handles.contours_positions{handles.frame+1});
            handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{handles.frame+1};
            handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{handles.frame+1};
            handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{handles.frame+1};
        end
    end  
end

% Show and edit contours
try
    % Edit contour lines
    for H = [handles.axes1,handles.axes2,handles.axes3]
        h = imcline(handles.contours{handles.frame},H);
        h.Enable = 'on';
        h.Visible = 'on';
        h.IndependentDrag{2} = 'on';
        [h.Appearance(1:2).Color] = deal('r','g');
        [h.Appearance(1:2).MarkerFaceColor] = deal('r','g');
        [h.Appearance(1:2).MarkerSize] = deal(10, 10);
        [h.Appearance(1:2).Marker] = deal('o', 'o');
        [h.Appearance(1:2).LineWidth] = deal(2, 2);
    end
    linkaxes([handles.axes1,handles.axes2,handles.axes3]);
    iptPointerManager(handles.figure1,'enable')
catch
end

% Update handles object
guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function add_contours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get contours curves
for c=1:handles.contours_number
    curves = getcline(handles.axes1);
    handles.positions{c} = curves.Position{1,1};
end

% create contour line (cline) object
handles.contours{handles.frame} = cline(handles.positions);

% Edit contour lines
for H = [handles.axes1,handles.axes2,handles.axes3]
    h = imcline(handles.contours{handles.frame},H);
    h.Enable = 'on';
    h.Visible = 'on';
    h.IndependentDrag{2} = 'on';
    [h.Appearance(1:2).Color] = deal('r','g');
    [h.Appearance(1:2).MarkerFaceColor] = deal('r','g');
    [h.Appearance(1:2).MarkerSize] = deal(10, 10);
    [h.Appearance(1:2).Marker] = deal('o', 'o');
    [h.Appearance(1:2).LineWidth] = deal(2, 2);
end
linkaxes([handles.axes1,handles.axes2,handles.axes3]);
iptPointerManager(handles.figure1,'enable')

% Update handles object
guidata(hObject, handles);
uiwait()


function x_min = edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double  
x_min = str2double(get(handles.edit2,'String'));

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_max = edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
x_max = str2double(get(handles.edit3,'String'));
  
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y_min = edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

y_min = str2double(get(handles.edit4,'String'));

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y_max = edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
y_max = str2double(get(handles.edit5,'String'));

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function update_axis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update axis
x_min = edit2_Callback(hObject, eventdata, handles);
x_max = edit3_Callback(hObject, eventdata, handles); 
y_min = edit4_Callback(hObject, eventdata, handles);
y_max = edit5_Callback(hObject, eventdata, handles);
handles.axis = [x_min x_max y_min y_max];

% Set axis on axes1 and axes 2
axis([handles.axes1,handles.axes2,handles.axes3],handles.axis)

% Update handles object
guidata(hObject, handles);



function c_min = edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
c_min = str2double(get(handles.edit6,'String'));


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c_max = edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
c_max = str2double(get(handles.edit7,'String'));


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function update_caxis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update axis
c_min = edit6_Callback(hObject, eventdata, handles);
c_max = edit7_Callback(hObject, eventdata, handles); 
handles.caxis = [c_min c_max];

% Set axis on axes1 and axes 2
caxis(handles.axes1,handles.caxis);

% Update handles object
guidata(hObject, handles);


% --- Executes on button press in pushbutton11.
function duplicate_contours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get frame from which contours will be duplicated
c_frame = edit8_Callback(hObject, eventdata, handles);

% Replace contours
if ~isempty(handles.contours_positions{c_frame})
    % Get contours
    handles.contours{handles.frame} = cline(handles.contours_positions{c_frame});
    handles.contours{handles.frame}.IsClosed = handles.contours_isclosed{c_frame};
    handles.contours{handles.frame}.IsCorner = handles.contours_iscorner{c_frame};
    handles.contours{handles.frame}.IsCurved = handles.contours_iscurved{c_frame};
end

% Show and edit contours
try
    % Edit contour lines
    for H = [handles.axes1,handles.axes2,handles.axes3]
        h = imcline(handles.contours{handles.frame},H);
        h.Enable = 'on';
        h.Visible = 'on';
        h.IndependentDrag{2} = 'on';
        [h.Appearance(1:2).Color] = deal('r','g');
        [h.Appearance(1:2).MarkerFaceColor] = deal('r','g');
        [h.Appearance(1:2).MarkerSize] = deal(10, 10);
    end
    linkaxes([handles.axes1,handles.axes2,handles.axes3]);
    iptPointerManager(handles.figure1,'enable')
catch
end

% Update handles object
guidata(hObject, handles);


function c = edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
c = str2double(get(handles.edit8,'String'));


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
