function varargout = plotIQMP(varargin)
% plotIQMP - allows to compare simulated experiments data to measurements
%
% USAGE:
% ======
% [] = plotIQMP(plotdata)
%
% plotdata: This datastructure is the output argument of the function 
% IQMcomparemeasurements. 

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plotIQMP_OpeningFcn, ...
                   'gui_OutputFcn',  @plotIQMP_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERFACE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before plotIQMP is made visible.
function plotIQMP_OpeningFcn(hObject, eventdata, handles, varargin)

if nargin ~= 4,
    error('Incorrect number of input arguments.');
end
handles.plotdata = varargin{1};
% set modelselection and choose first model
set(handles.modelselection,'String',{handles.plotdata.model.name});
set(handles.modelselection,'Value',1);
% set experimentselection for first model and first experiment
set(handles.experimentselection,'String',{handles.plotdata.model(1).experiment.name});
set(handles.experimentselection,'Value',1);
% select plottype 
handles.dataPlotType = 'plot';     
% set errorbarflag to 1
handles.errorbars = 1;
% Initialize export figure handle and grid flag
handles.exportFigureHandle = [];
handles.grid = 0;
% Set up axes
plotdata = handles.plotdata;
m = get(handles.modelselection,'Value');
e = get(handles.experimentselection,'Value');
componentNames = cell(0, 0);
for experimentIndex = 1 : numel(e)
    edata = plotdata.model(m).experiment(e(experimentIndex));
    componentNames = [componentNames, edata.componentnames];
end
handles.componentNames = unique(componentNames);
handles.maxComponents = numel(componentNames);
cols = ceil(sqrt(handles.maxComponents));
rows = ceil(handles.maxComponents/cols);
position = get(handles.plotarea, 'position');
xStart = position(1);
yStart = position(2);
width = position(3);
height = position(4);
row = 1;
col = 1;
handles.axesArray = cell(handles.maxComponents, 1);
for i = 1 : handles.maxComponents
   row = ceil(i/cols);
   col = i - cols * (row-1);
   handles.axesArray{i} = subplot('Position', [ (width/cols)*(col-1)+xStart, (height/rows)*(rows-row)+yStart, (width/cols)*0.85, (height/rows)*0.82]);
end
% Doing a first plot
doPlot(handles);
% Choose default command line output for plotIQMP
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
return

% --- Outputs from this function are returned to the command line.
function varargout = plotIQMP_OutputFcn(hObject, eventdata, handles) 
%varargout{1} = handles.output;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function doPlot(handles)
if isequal(get(handles.ThreeDMode, 'Value'), true)
   doPlot3D(handles); 
else
   doPlot2D(handles);
end
return

function doPlot3D(handles)

return

function doPlot2D(handles)

colorvector = {'b','g','r','c','m','y','k'};
markervector = {'o','x','+','*','s','d','v','^','<','>','p'};
warning off;
% get the data to plot
plotdata = handles.plotdata;
m = get(handles.modelselection,'Value');
e = get(handles.experimentselection,'Value');

for experimentIndex = 1 : numel(e)
    edata = plotdata.model(m).experiment(e(experimentIndex));
    % general information
    titletext = regexprep(edata.name,'_',' ');
    xlabeltext = 'Time';
    % simulated data information
    sim_timevector = edata.timevector;
    sim_componentnames = edata.componentnames;
    sim_componentvalues = edata.componentvalues;
    if isempty(sim_componentvalues),
        error('Please check if the names of the measured data appear in the model.');
    end
    % plot simulated data
    for k=1:size(sim_componentvalues,2),
        subplot(handles.axesArray{k})
        try
            feval(handles.dataPlotType,sim_timevector,sim_componentvalues(:,k),'linewidth',2,'color',colorvector{ mod(experimentIndex,numel(colorvector)) + 1} ); hold on;
        catch MException
            ''
        end
    end
    % measured data information 
    for meas=1:length(edata.measurement),
        meas_timevector = edata.measurement(meas).timevector;
        meas_componentnames = edata.measurement(meas).componentnames;
        meas_componentvalues = edata.measurement(meas).componentvalues;   
        meas_maxvalues = edata.measurement(meas).maxvalues;   
        meas_minvalues = edata.measurement(meas).minvalues;   
    %     marker = markervector{mod(meas-1,length(markervector))+1};
    %     feval(handles.dataPlotType,meas_timevector,meas_componentvalues,['--' marker]); hold on;
        for k=1:size(sim_componentvalues,2),
            subplot(handles.axesArray{k});
            feval(handles.dataPlotType,meas_timevector,meas_componentvalues(:,k),['*:'],'linewidth',2,'color',colorvector{mod(experimentIndex,numel(colorvector)) + 1},  'marker',markervector{mod(experimentIndex,numel(markervector)) + 1}); hold on;
        end
        if handles.errorbars == 1 && strcmp(handles.dataPlotType,'plot'),
            % plot error bounds
            for k=1:handles.maxComponents,
                subplot(handles.axesArray{k})
                color = colorvector{mod(experimentIndex,numel(colorvector)) + 1};
    %             for k1 = 1:size(meas_timevector,1),
                for k1 = 1:length(meas_timevector),
                    feval(handles.dataPlotType,[meas_timevector(k1),meas_timevector(k1)],[meas_minvalues(k1,k),meas_maxvalues(k1,k)],['.:',color]);
                end
            end
        end
    end
end
titleTexts = cell(numel(e), 1);
for experimentIndex = 1 : numel(e)
    
    edata = plotdata.model(m).experiment(e(experimentIndex));
    % general information
    titleTexts{experimentIndex*2-1} = [regexprep(edata.name,'_',' '), ' fit'];
    titleTexts{experimentIndex*2} = [regexprep(edata.name,'_',' '), ' measurement'];

end

for k = 1 : numel(handles.axesArray)
    subplot(handles.axesArray{k});
    hold off;

    hlhlx = legend(titleTexts);
    set(hlhlx,'Interpreter','none');

    hlhlx = title(sim_componentnames{k});
    
    set(hlhlx,'Interpreter','none');
    hlhlx = xlabel(xlabeltext);
    set(hlhlx,'Interpreter','none');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT FIGURE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export the figure
function export_Callback(hObject, eventdata, handles)
warning off;
if isempty(handles.exportFigureHandle),
    figH = figure;
    handles.exportFigureHandle = figH;
    % Update handles structure
    guidata(hObject, handles);
else
    figH = handles.exportFigureHandle;
    figure(figH);
end
nrow = str2num(get(handles.nrow,'String'));
ncol = str2num(get(handles.ncol,'String'));
nnumber = str2num(get(handles.nnumber,'String'));
subplot(nrow,ncol,nnumber);
doPlot(handles);
if handles.grid == 1,
    grid;
end
% set axes
XLim = get(handles.plotarea,'Xlim');
YLim = get(handles.plotarea,'Ylim');
axis([XLim, YLim]);
return

% Request new figure for export
function newexportfigure_Callback(hObject, eventdata, handles)
handles.exportFigureHandle = [];
% Update handles structure
guidata(hObject, handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIVIAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in zoombutton.
function errorbarbutton_Callback(hObject, eventdata, handles)
% toogle the errorbars in the figure
if handles.errorbars == 0,
    handles.errorbars = 1;
else
    handles.errorbars = 0;
end
% plot
doPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% toogle the zoom in the figure
zoom
return

% --- Executes on button press in gridbutton.
function gridbutton_Callback(hObject, eventdata, handles)
% toogle the grid in the figure
grid
if handles.grid == 1,
    handles.grid = 0;
else
    handles.grid = 1;
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in modelselection.
function modelselection_Callback(hObject, eventdata, handles)
try
    modelindex = get(handles.modelselection,'Value');
    set(handles.experimentselection,'String',{handles.plotdata.model(modelindex).experiment.name});
    doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
return

% --- Executes on selection change in experimentselection.
function experimentselection_Callback(hObject, eventdata, handles)
try
    doPlot(handles);
catch MException
    errordlg('This selection is not possible.','Error','on');               
end
return

% --- From R2014B the radiobutton groups are handled differently ...
function plotAxesSelection_SelectionChangeFcn(hObject, eventdata, handles)
handles.dataPlotType = eventdata.NewValue.String;
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'plot';
% disable errorbarbutton
set(handles.errorbarbutton,'Visible','on');
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in loglog.
function semilogx_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogx';
% disable errorbarbutton
set(handles.errorbarbutton,'Visible','off');
handles.errorbars = 0;
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in semilogx.
function semilogy_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogy';
% disable errorbarbutton
set(handles.errorbarbutton,'Visible','off');
handles.errorbars = 0;
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in loglog.
function loglog_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'loglog';
% disable errorbarbutton
set(handles.errorbarbutton,'Visible','off');
handles.errorbars = 0;
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

function multiExperimentSelect_Callback(hObject, eventdata, handles)
if isequal(get(hObject, 'Value'), true)
    set(handles.experimentselection, 'max', 2);
else
    set(handles.experimentselection, 'max', 1);
end
guidata(hObject, handles);
return