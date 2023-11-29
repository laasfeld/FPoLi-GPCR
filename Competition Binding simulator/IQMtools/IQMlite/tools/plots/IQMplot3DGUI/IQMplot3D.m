function varargout = IQMplot3D(varargin)
% IQMplot3D - plots given data.
%
% USAGE:
% ======
% [] = IQMplot3D(time,data)
% [] = IQMplot3D(time,data,names)
% [] = IQMplot3D(time,data,names,name)
% [] = IQMplot3D(time,data,names,legendtext,name)
% [] = IQMplot3D(time,data,names,legendtext,marker,name)
% [] = IQMplot3D(time,data,names,errorindices,minvalues,maxvalues,legendtext,marker,name)
%
% [] = IQMplot3D(datastruct1)
% [] = IQMplot3D(datastruct1,datastruct2)
% [] = IQMplot3D(datastruct1,datastruct2, ..., datastructN)
%
% The datastructures are created most easily using the function
% createdatastructIQMplotIQM.
%
% time: column vector with time information
% data: matrix with data where each row corresponds to one time point and
%   each column to a different variable
% names: cell-array with the names of the data variables
% legendtext: cell-array of same length as names with text to be used for
%   the legend.
% marker: marker and line style for plot
% errorindices: indices of the data for which errorbounds are available
% minvalues: error bounds for data ... to be shown by error bars
% maxvalues: error bounds for data ... to be shown by error bars
% name: name describing the datastruct
%
% datastruct: datastructure with all the plotting data (allows for
%   displaying several datastructs at a time in the same GUI).
%
% DEFAULT VALUES:
% ===============
% names: the plotted variables obtain the name 'x1', 'x2', ...
% legendtext: same as names
% marker: '-'
% min/maxvalues: no errorbars shown

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IQMplot3D_OpeningFcn, ...
                   'gui_OutputFcn',  @IQMplot3D_OutputFcn, ...
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

% --- Executes just before IQMplot3D is made visible.
function IQMplot3D_OpeningFcn(hObject, eventdata, handles, varargin)
% check if datastructure or normal data as input
handles.dataSets = varargin{1};
if nargin > 4
    handles.measureSets = varargin{2};
end
handles = switchDataSet(handles,1);     % switch to first datastruct
% Initialize datastructs pulldown menu
datastructnames = {};
for k = 1:length(handles.dataSets),
    datastructnames{k} = handles.dataSets{k}.name;
end
set(handles.datastructs,'String',datastructnames);
% select plottype to start with
handles.dataPlotType = 'plot';          
% Initialize export figure handle
handles.exportFigureHandle = [];
handles.grid = 0;
% Doing a first plot
doPlot(handles);
% Choose default command line output for IQMplot3D
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
return

% --- Executes just before IQMplot3D is made visible.
function Exit_Callback(hObject, eventdata, handles, varargin)
clear global doRemoveZeroComponentFlag
closereq
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWITCH GIVEN DATASTRUCTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles] = switchDataSet(handles,indexDataSet)
dataSet = handles.dataSets(indexDataSet);
measureSet = handles.measureSets(indexDataSet);
% Set all the plot data also in the handles structure to be accessed by 
% all callback functions
% get number of yaxisdata in old plot
if isfield(handles,'dataNames'),
    znumberold = length(handles.dataNames);
else
    znumberold = 0;
end
handles.measureTime = cell(size(measureSet));
handles.measureData = cell(size(measureSet));

handles.yAxisTimepoint = 0;
handles.measureTime{1} = measureSet{1}.timevector';
handles.measureData{1} = measureSet{1}.componentvalues;
handles.minvalues = measureSet{1}.minvalues;
handles.maxvalues = measureSet{1}.maxvalues;
handles.time = dataSet{1}.time;
handles.data = dataSet{1}.data;
handles.dataNames = dataSet{1}.dataNames;
handles.legentext = dataSet{1}.legendtext;
handles.marker = '*'%dataSet{1}.marker;
handles.errorindices = dataSet{1}.errorindices;
handles.name = dataSet{1}.name;

for dataSetIndex = 2 : numel(dataSet)
    handles.time = [handles.time, dataSet{dataSetIndex}.time];
    handles.data = cat(3, handles.data, dataSet{dataSetIndex}.data);
    handles.measureTime{dataSetIndex} = measureSet{dataSetIndex}.timevector';
    handles.measureData{dataSetIndex} = measureSet{dataSetIndex}.componentvalues;
    %handles.dataNames = [handles.dataNames, dataSet{dataSetIndex}.dataNames];
    %handles.legentext = [handles.legentext, dataSet{dataSetIndex}.legendtext];
    %handles.marker = [handles.marker, dataSet{dataSetIndex}.marker];
    handles.minvalues = [handles.minvalues; measureSet{dataSetIndex}.minvalues];
    handles.maxvalues = [handles.maxvalues; measureSet{dataSetIndex}.maxvalues];
    handles.errorindices = [handles.errorindices; dataSet{dataSetIndex}.errorindices];
    %handles.name = [handles.name, dataSet{dataSetIndex}.name];
end
% update selection menu
set(handles.xaxisselection,'String',{'TIME',handles.dataNames{:}});
set(handles.zaxisselection,'String',handles.dataNames);
stateNames = cell(0,0);
for dataName = 1 : numel(handles.dataNames);
    if strfind(handles.dataNames{dataName}, '(state)')
        stateNames{end + 1} = handles.dataNames{dataName};
    end
end
%set(handles.yaxisselection, 'String', stateNames);
set(handles.yaxisselection,'String',handles.dataNames);
set(handles.yaxisselection, 'Value', 1);
set(handles.xaxisselection,'Value',1);
% change selection only if unequal numbers of data in the sets
if znumberold ~= length(handles.dataNames),
    set(handles.zaxisselection,'Value',1);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATASTRUCTS SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function datastructs_Callback(hObject, eventdata, handles)
dataSetIndex = get(handles.datastructs,'Value');
handles = switchDataSet(handles,dataSetIndex);
%doPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

% --- Outputs from this function are returned to the command line.
function varargout = IQMplot3D_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
return

% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% toogle the zoom in the figure
zoom
return

% --- Executes on selection change in xaxisselection.
function xaxisselection_Callback(hObject, eventdata, handles)
try
    %doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in yaxisselection.
function yaxisselection_Callback(hObject, eventdata, handles)
try
    %doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in zaxisselection.
function zaxisselection_Callback(hObject, eventdata, handles)
try
    %doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on edit of yAxisTimepointSelection.
function yAxisTimepointSelection_Callback(hObject, eventdata, handles)
    handles.yAxisTimepoint = str2num(get(hObject, 'String'));
    guidata(hObject, handles);
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

% --- From R2014B the radiobutton groups are handled differently ...
function plotAxesSelection_SelectionChangeFcn(hObject, eventdata, handles)
handles.dataPlotType = eventdata.NewValue.String;
% Update handles structure
guidata(hObject, handles);
%doPlot(handles);
return


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'plot';
% Update handles structure
guidata(hObject, handles);
%doPlot(handles);
return

% --- Executes on button press in loglog.
function semilogx_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogx';
% Update handles structure
guidata(hObject, handles);
%doPlot(handles);
return

% --- Executes on button press in semilogx.
function semilogy_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogy';
% Update handles structure
guidata(hObject, handles);
%doPlot(handles);
return

% --- Executes on button press in loglog.
function loglog_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'loglog';
% Update handles structure
guidata(hObject, handles);
%doPlot(handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export_Callback(hObject, eventdata, handles)
warning off;
if isempty(handles.exportFigureHandle),
    figH = figure;
    handles.exportFigureHandle = figH;
    copyobj(handles.plotarea, figH);
    % Update handles structure
    guidata(hObject, handles);
else
    figH = handles.exportFigureHandle;
    figure(figH);
    copyobj(handles.plotarea, figH);
end
%nrow = str2num(get(handles.nrow,'String'));
%ncol = str2num(get(handles.ncol,'String'));
%nnumber = str2num(get(handles.nnumber,'String'));
%subplot(nrow,ncol,nnumber);
%doPlot(handles);
%handles.exportFigureHandle.Children(3).Children(1).FaceColor = 'interp';
%handles.exportFigureHandle.Children(3).Children(1).FaceAlpha = 0.8;
%if handles.grid == 1,
%    grid;
%end
% set axes
%XLim = get(handles.plotarea,'Xlim');
%YLim = get(handles.plotarea,'Ylim');
%axis([XLim, YLim]);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUEST NEW EXPORT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newexportfigure_Callback(hObject, eventdata, handles)
handles.exportFigureHandle = [];
% Update handles structure
guidata(hObject, handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotNow_Callback(hObject, eventdata, handles)
doPlot(handles);

return

function pushbutton9_Callback(hObject, eventdata, handles)
[name, path] = uiputfile('.mat');
save([path, name], 'handles');
    
return

function doPlot(handles)
warning off;
colorvector = {'b','g','r','c','m','y','k'};
time = handles.time;
data = handles.data;
measureTime = handles.measureTime;
measureData = handles.measureData;
dataNames = handles.dataNames;
errorindices = handles.errorindices;
maxvalues = handles.maxvalues;
minvalues = handles.minvalues;
xaxis = handles.xaxisselection;
zaxis = handles.zaxisselection;



yAxisParams = get(handles.yaxisselection, 'String');
yAxisParam = yAxisParams{get(handles.yaxisselection, 'Value')};
paramIndexY = [];

zParams = get(zaxis, 'String');

for componentIndex = 1 : numel(get(zaxis, 'String'))
    if strcmp(yAxisParam, zParams{componentIndex})
        paramIndexY = componentIndex;
        break;
    end
end



% get variable that is chosen for the x-axis
indexX = get(xaxis,'Value');
% get variables that are chosen for the z-axis
indexZ = get(zaxis,'Value');
yvariables = zeros(size(time, 1), size(time, 2));
counter = 1;
if isequal(get(handles.holdOnOption, 'Value'), true)
    hold on;
else
    hold off;
end
xvariables = time;
%for expIndex = 1 : size(time, 2)
%    yvariables(:, expIndex) = ones(numel(1 : size(time, 1)),1)*data(handles.yAxisTimepoint+1, paramIndexY, expIndex);
%end
yvariables = squeeze(data(:,paramIndexY, :));
zvariables = squeeze(data(:,indexZ, :));




for expIndex = 1 : size(time, 2)
    plot3(time(1 : size(time, 1) , expIndex), ones(numel(1 : size(time, 1)),1)*data(handles.yAxisTimepoint+1, paramIndexY, expIndex), data(1 : size(time, 1), indexZ, expIndex), '-b','MarkerEdgeColor', colorvector{1} );
    hold on;
end

zAxisParams = get(zaxis, 'String');
zAxisParam = strsplit(zAxisParams{indexZ}, ' ');
zAxisParam = zAxisParam(1);
measureIndexZ = [];
for componentIndex = 1 : numel(handles.measureSets{1}.componentnames)
    if strcmp(zAxisParam, handles.measureSets{1}.componentnames{componentIndex})
        measureIndexZ = componentIndex;
        break;
    end
end

times = [];
measurements = [];
yAxisData = [];

if ~isempty(measureIndexZ)
    for expIndex = 1 : numel(measureTime)
        times = [times; measureTime{expIndex}];
        yAxisData = [yAxisData; ones(numel(measureTime{expIndex}), 1).*data(handles.yAxisTimepoint+1, paramIndexY, expIndex)];
        measurements = [measurements; measureData{expIndex}(:, measureIndexZ)];
        %plot3(measureTime{expIndex}, ones(numel(measureTime{expIndex}))*data(1, paramIndexY, expIndex), measureData{expIndex}(:, measureIndexZ), '*g','MarkerEdgeColor', colorvector{2} );
        %plot3(measureTime(1 : size(measureTime, 1) , expIndex), ones(numel(1 : size(measureTime, 1)),1)*data(1, paramIndexY, expIndex), measureData(1 : size(measureTime, 1), measureIndexZ, expIndex), '*g','MarkerEdgeColor', colorvector{2} );
        hold on;
    end
end
plot3(times, yAxisData, measurements, '*g','MarkerEdgeColor', colorvector{2});

% select linewidth
if ~isempty(errorindices),
    % wider line in case of data with error bounds
    addOption = sprintf(',''linewidth'',2');
else
    addOption = '';
end

% plot
%plot3(xvariables, yvariables, zvariables, '*','MarkerEdgeColor', colorvector{1});
try
    surf(xvariables, yvariables, zvariables, 'edgecolor','none')
    %children = get(handles.plotarea.Children);
    %children(1).FaceColor;
    handles.plotarea.Children(1).FaceColor = 'interp';
    handles.plotarea.Children(1).FaceAlpha = 0.8;
catch
    plot3(xvariables, yvariables, zvariables, '*','MarkerEdgeColor', colorvector{1});
end
%eval(sprintf('feval(handles.dataPlotType,xvariable,zvariables,handles.marker%s);',addOption))


% plot error bounds 


hlhlx = legend(handles.legentext(indexZ)); 
%set(hlhlx,'Interpreter','none');
% write axis labels
if indexX == 1,
    xlabel('Time');
else
    hlhlx = xlabel(dataNames(indexX-1));
 %   set(hlhlx,'Interpreter','none');
end
ylabel(yAxisParam);
zlabel(zAxisParam);
% write title (name)
hlhlx = title(handles.name);
%set(hlhlx,'Interpreter','none');
return

function plotTypeSelection_Callback(hObject,eventdata,handles)
handles.dataPlotType = get(hObject, 'String');
%doPlot(handles);

return

