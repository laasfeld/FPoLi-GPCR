function varargout = AdvancedGroupingGUI(varargin)
% ADVANCEDGROUPINGGUI MATLAB code for AdvancedGroupingGUI.fig
%      ADVANCEDGROUPINGGUI, by itself, creates a new ADVANCEDGROUPINGGUI or raises the existing
%      singleton*.
%
%      H = ADVANCEDGROUPINGGUI returns the handle to a new ADVANCEDGROUPINGGUI or the handle to
%      the existing singleton*.
%
%      ADVANCEDGROUPINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDGROUPINGGUI.M with the given input arguments.
%
%      ADVANCEDGROUPINGGUI('Property','Value',...) creates a new ADVANCEDGROUPINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdvancedGroupingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdvancedGroupingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdvancedGroupingGUI

% Last Modified by GUIDE v2.5 13-May-2020 14:43:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdvancedGroupingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AdvancedGroupingGUI_OutputFcn, ...
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


% --- Executes just before AdvancedGroupingGUI is made visible.
function AdvancedGroupingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AdvancedGroupingGUI (see VARARGIN)

% Choose default command line output for AdvancedGroupingGUI
handles.modelStates = varargin{1};
handles.experimentStates = varargin{2};
handles.parameterNames = varargin{3};
handles.reactionNames = varargin{4};
handles.reactionFormulae = varargin{5};
handles.ODE = varargin{6};
handles.experimentNames = varargin{7};
handles.paramInfo = varargin{8};
handles.parameterNominalValues = varargin{9};
handles.statesInSpecificExperiments = varargin{10};
modelStateArray = cell(numel(handles.experimentStates), 1);
modelStateArray(:) = handles.modelStates(1);
set(handles.stateTypeTable, 'data', [handles.experimentStates', modelStateArray]);
set(handles.stateTypeTable,'ColumnFormat', {'numeric', handles.modelStates'});
handles.output = hObject;
set(handles.stateListbox, 'String', handles.modelStates);
handles.currentState = 1;
handles.relevantParamDependance = cell(numel(handles.modelStates), 1);
handles.selectedCells = [];
for i = 1 : numel(handles.modelStates)
    handles.relevantParamDependance{i} = true(1, numel(getRelevantParameterIndices(handles, i)));
end
stateListbox_Callback(handles.stateListbox, eventdata, handles);
% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);
% UIWAIT makes AdvancedGroupingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function relevantParamIndices = getRelevantParameterIndices(handles, stateIndex)

ODEParams = symvar(handles.ODE{stateIndex});
[unneeded, presentEquationIndices, equationIndicesInODE] = intersect(handles.reactionNames, ODEParams);
paramsInEquations = cell(1, 0);
for reactionIndex = presentEquationIndices'
    % only parameters which are multiplied by a particular state conc are
    % relevant
    equationParts = strsplit(handles.reactionFormulae{reactionIndex}, '-'); % this logic might not be completely general for all cases
    for equationPart = 1 : numel(equationParts)
        vars = symvar(equationParts{equationPart});
        intersection = intersect(handles.modelStates{stateIndex}, vars);
        if strcmp(intersection, handles.modelStates{stateIndex});
             paramsInEquations = [paramsInEquations; vars];
        end
    end
   
end
%remove states
[unneeded, presentStateIndices, unneeded2] = intersect(paramsInEquations, handles.modelStates);
paramsInEquations(presentStateIndices) = [];
paramsInEquations = unique(paramsInEquations);
paramsInODE = ODEParams;
paramsInODE(equationIndicesInODE) = [];
relevantParams = unique([paramsInODE; paramsInEquations]);
[unneeded, unneeded2, relevantParamIndices] = intersect(relevantParams, handles.parameterNames);



% --- Outputs from this function are returned to the command line.
function varargout = AdvancedGroupingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
data = get(handles.parameterGroupingTable, 'data');
data(:, 1) = [];
stateData = get(handles.stateTypeTable, 'data');
varargout{1} = data;
varargout{2} = stateData;
varargout{3} = handles.minMaxTables;
delete(handles.figure1);

% --- Executes on button press in generateParameterGroupingTable.
function generateParameterGroupingTable_Callback(hObject, eventdata, handles)
% hObject    handle to generateParameterGroupingTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.parameterGroupingTable.set('ColumnName', ['Experiment'; handles.parameterNames]);
numData = cell(numel(handles.experimentNames), numel(handles.parameterNames));
numData(:) = {'1'};
relevantStates = cell(numel(handles.parameterNames), 1);
% find out which states nature the parameter depends upon
for parameter = 1 : numel(handles.parameterNames)
   relevantStates{parameter} = cell(1, 0);
   for state = 1 : numel(handles.relevantParamDependance)
       relevantParamIndices = getRelevantParameterIndices(handles, state);
       paramTruthValues = handles.relevantParamDependance{state};
       correctParamIndices = relevantParamIndices(paramTruthValues);
       for i = 1 : numel(correctParamIndices)
          if isequal(parameter, correctParamIndices(i))
              relevantStates{parameter}{end + 1} = handles.modelStates{state};
              break; 
          end
       end
   end
end

stateTypeTableData = get(handles.stateTypeTable, 'data');
relevantExperimentStates = cell(1, numel(handles.parameterNames));
allFoundCombinations = cell(numel(handles.parameterNames), 1);
for parameter = 1 : numel(handles.parameterNames)
    relevantExperimentStates{parameter} = cell(numel(relevantStates{parameter}), 1);
    for relevantStateIndex = 1 : numel(relevantStates{parameter})
        relevantExperimentStates{parameter}{relevantStateIndex} = cell(1, 0);
        for experimentStateIndex = 1 : size(stateTypeTableData, 1)
            if strcmp(stateTypeTableData{experimentStateIndex, 2}, relevantStates{parameter}{relevantStateIndex})
                relevantExperimentStates{parameter}{relevantStateIndex}{end + 1} = stateTypeTableData{experimentStateIndex, 1};
            end
        end
        relevantExperimentStates{parameter}{relevantStateIndex} = unique(relevantExperimentStates{parameter}{relevantStateIndex});
        if isempty(relevantExperimentStates{parameter}{relevantStateIndex})
           relevantExperimentStates{parameter}{relevantStateIndex} = cell(1, 1);
        end
    end
    if ~isempty(relevantExperimentStates{parameter})
        allFoundCombinations{parameter} = allcomb(relevantExperimentStates{parameter}{:});
    end
end

%find how many different groups should exist for each parameter
for experimentIndex = 1 : numel(handles.statesInSpecificExperiments)
    correspondingStatesInModel = cell(numel(handles.statesInSpecificExperiments{experimentIndex}), 1);      
    for stateInExperiment = 1 : numel(handles.statesInSpecificExperiments{experimentIndex})
        for experimentStateIndex = 1 : size(stateTypeTableData, 1)
            if strcmp(stateTypeTableData{experimentStateIndex, 1}, handles.statesInSpecificExperiments{experimentIndex}{stateInExperiment})
                correspondingStatesInModel{stateInExperiment} = stateTypeTableData{experimentStateIndex, 2};
                break;
            end 
        end
    end
    for parameter = 1 : numel(handles.parameterNames)
        relevantStateArray = cell(numel(relevantStates{parameter}), 1);
        for i = 1 : numel(relevantStates{parameter})
           for k = 1 : numel(correspondingStatesInModel)
               if strcmp(relevantStates{parameter}{i}, correspondingStatesInModel{k})
                   relevantStateArray{i} = handles.statesInSpecificExperiments{experimentIndex}{k};
               end
           end
        end
        indices = 1 : numel(allFoundCombinations{parameter});
        for i = 1 : numel(relevantStates{parameter})
            if isempty(relevantStateArray{i})

            else
                indices = intersect(indices, find(cellfun(@isempty, strfind(allFoundCombinations{parameter}(:, i), relevantStateArray{i})) == 0));
            end
        end
        
        if(isempty(indices))
            numData{experimentIndex, parameter} = 'standard';
        else
            groupIndex = indices(1); % if there are more indices then this experiment does not depend on the parameter value at all and this can be added to any group
            foundCombinations = allFoundCombinations{parameter};
            for i = 1 : numel(foundCombinations(indices(1), :))
                if isempty(foundCombinations{indices(1), i})
                    foundCombinations{indices(1), i} = relevantStates{parameter}{i};
                end
            end
            groupName = foundCombinations{indices(1), 1};           
            for i = 2 : numel(foundCombinations(indices(1), :))
                groupName = [groupName, ';',foundCombinations{indices(1), i}];
            end
            numData{experimentIndex, parameter} = groupName;%num2str(groupIndex);
        end
    end
end

data = [handles.experimentNames', numData];
handles.parameterGroupingTable.set('data', data);
handles.parameterGroupingTable.set('columnEditable', [false, true(1, size(data, 2))]);
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);
if isequal(get(handles.figure1, 'waitstatus'),'waiting')
    uiresume(handles.figure1)
else
    delete(handles.figure1);
end


% --- Executes on selection change in stateListbox.
function stateListbox_Callback(hObject, eventdata, handles)
% hObject    handle to stateListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stateListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stateListbox
contents = cellstr(get(hObject,'String'));
selectedState = get(hObject,'Value');
relevantParameterIndices = getRelevantParameterIndices(handles, selectedState);

data = get(handles.stateTypeTable, 'data');
columnFormat = cell(1, numel(handles.parameterNames));
columnFormat(:) = {'logical'};
columnEditable = true(1, numel(handles.parameterNames));
handles.parameterStateDependanceTable.set('ColumnName', handles.parameterNames(relevantParameterIndices));
handles.parameterStateDependanceTable.set('ColumnEditable', true(1, numel(relevantParameterIndices)));
format = cell(1, numel(relevantParameterIndices));
format(:) = {'logical'};
%handles.parameterStateDependanceTable.set('ColumnFormat', format);
%handles.parameterStateDependanceTable.set('data', true(numel(handles.modelStates), numel(handles.parameterNames)))
%handles.parameterStateDependanceTable.set('data', true(numel(handles.modelStates), numel(handles.parameterNames)))
handles.currentState = selectedState;
data = handles.relevantParamDependance{handles.currentState};

handles.parameterStateDependanceTable.set('data', data)%, 'RowName', handles.modelStates, 'ColumnName', handles.parameterNames, 'ColumnFormat', columnFormat, 'ColumnEditable', columnEditable);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stateListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stateListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in parameterStateDependanceTable.
function parameterStateDependanceTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to parameterStateDependanceTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.relevantParamDependance{handles.currentState} = get(hObject, 'data');
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in stateTypeTable.
function stateTypeTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to stateTypeTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setGroupOfSelectedWells.
function setGroupOfSelectedWells_Callback(hObject, eventdata, handles)
% hObject    handle to setGroupOfSelectedWells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
'a'
data = get(handles.parameterGroupingTable, 'data');
for cell = 1 : size(handles.selectedCells, 1);
    data{handles.selectedCells(cell, 1), handles.selectedCells(cell, 2)} = get(handles.setSelectedGroupBox, 'String');
end
set(handles.parameterGroupingTable, 'data', data);
guidata(hObject, handles);


function setSelectedGroupBox_Callback(hObject, eventdata, handles)
% hObject    handle to setSelectedGroupBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setSelectedGroupBox as text
%        str2double(get(hObject,'String')) returns contents of setSelectedGroupBox as a double


% --- Executes during object creation, after setting all properties.
function setSelectedGroupBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setSelectedGroupBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected cell(s) is changed in parameterGroupingTable.
function parameterGroupingTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to parameterGroupingTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.selectedCells = eventdata.Indices;
guidata(hObject, handles);


% --- Executes on selection change in chooseParamForMinMax.
function chooseParamForMinMax_Callback(hObject, eventdata, handles)
% hObject    handle to chooseParamForMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chooseParamForMinMax contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chooseParamForMinMax
set(handles.parameterMinMaxTable, 'data', handles.minMaxTables{get(hObject,'Value')});
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function chooseParamForMinMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooseParamForMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in chooseGroupForNominality.
function chooseGroupForNominality_Callback(hObject, eventdata, handles)
% hObject    handle to chooseGroupForNominality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chooseGroupForNominality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chooseGroupForNominality


% --- Executes during object creation, after setting all properties.
function chooseGroupForNominality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooseGroupForNominality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lockAndConfirm.
function lockAndConfirm_Callback(hObject, eventdata, handles)
% hObject    handle to lockAndConfirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.parameterNominalValues
groupingData = get(handles.parameterGroupingTable, 'data');
groupingData(:, 1) = [];
handles.minMaxTables = cell(1, size(groupingData, 2));
for parameterIndex = 1 : size(groupingData, 2)
   groupIndices = groupingData(:, parameterIndex);
   groupNamesAsVector = cell(1, numel(groupIndices));
   for groupIndex = 1 : numel(groupIndices)
      groupNamesAsVector{groupIndex} = groupIndices{groupIndex};
   end
   uniqueNames = unique(groupNamesAsVector);
   handles.minMaxTables{parameterIndex} = [uniqueNames', num2cell([ones(numel(uniqueNames), 1)*handles.paramInfo.paramlow(parameterIndex), ones(numel(uniqueNames), 1)*handles.paramInfo.paramhigh(parameterIndex), ones(numel(uniqueNames), 1)*handles.parameterNominalValues(parameterIndex)])];
end
set(handles.parameterMinMaxTable, 'data', handles.minMaxTables{1});
set(handles.chooseParamForMinMax, 'String', handles.parameterNames);

%chooseGroupForNominality
guidata(hObject, handles);



% --- Executes when entered data in editable cell(s) in parameterMinMaxTable.
function parameterMinMaxTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to parameterMinMaxTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

handles.minMaxTables{get(handles.chooseParamForMinMax, 'Value')} = get(hObject, 'data');
handles = updateFullParameterTable(handles);
guidata(hObject, handles);

function returnFocus(hObject)
set(hObject, 'Enable', 'off');
drawnow;
set(hObject, 'Enable', 'on');


% --- Executes on button press in initializeFactor.
function initializeFactor_Callback(hObject, eventdata, handles)
% hObject    handle to initializeFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.forSelectedParams, 'Value')
    data = get(handles.parameterMinMaxTable, 'data');
    nominalValues = data(:, 4);
    minValues = num2cell(cell2num(nominalValues) * str2num(get(handles.factorLowBound, 'String')));
    maxValues = num2cell(cell2num(nominalValues) * str2num(get(handles.factorHighBound, 'String')));
    data(:, 2) = minValues;
    data(:, 3) = maxValues;
    set(handles.parameterMinMaxTable, 'data', data);
    handles.minMaxTables{get(handles.chooseParamForMinMax, 'Value')} = get(handles.parameterMinMaxTable, 'data');
else
    for parameter = 1 : numel(handles.minMaxTables)
        data = handles.minMaxTables{parameter};
        nominalValues = data(:, 4);
        minValues = num2cell(cell2num(nominalValues) * str2num(get(handles.factorLowBound, 'String')));
        maxValues = num2cell(cell2num(nominalValues) * str2num(get(handles.factorHighBound, 'String')));
        data(:, 2) = minValues;
        data(:, 3) = maxValues;
        set(handles.parameterMinMaxTable, 'data', data);
        handles.minMaxTables{parameter} = get(handles.parameterMinMaxTable, 'data');
    end
end
handles = updateFullParameterTable(handles);
guidata(hObject, handles);

function handles = updateFullParameterTable(handles)
finalTable = handles.minMaxTables{1};

for parameter = 2 : numel(handles.minMaxTables)
    finalTable = [finalTable; handles.minMaxTables{parameter}];
end
set(handles.fullParameterTable, 'data', finalTable);


function factorLowBound_Callback(hObject, eventdata, handles)
% hObject    handle to factorLowBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of factorLowBound as text
%        str2double(get(hObject,'String')) returns contents of factorLowBound as a double


% --- Executes during object creation, after setting all properties.
function factorLowBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to factorLowBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function factorHighBound_Callback(hObject, eventdata, handles)
% hObject    handle to factorHighBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of factorHighBound as text
%        str2double(get(hObject,'String')) returns contents of factorHighBound as a double


% --- Executes during object creation, after setting all properties.
function factorHighBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to factorHighBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in initializeAbsolute.
function initializeAbsolute_Callback(hObject, eventdata, handles)
% hObject    handle to initializeAbsolute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.forSelectedParams, 'Value')
    data = get(handles.parameterMinMaxTable, 'data');
    nominalValues = data(:, 4);
    minValues = num2cell(str2num(get(handles.absoluteLowBound, 'String')));
    maxValues = num2cell(str2num(get(handles.absoluteHighBound, 'String')));
    data(:, 2) = minValues;
    data(:, 3) = maxValues;
    set(handles.parameterMinMaxTable, 'data', data);
    handles.minMaxTables{get(handles.chooseParamForMinMax, 'Value')} = get(handles.parameterMinMaxTable, 'data');
else
    for parameter = 1 : numel(handles.minMaxTables)
        data = handles.minMaxTables{parameter};
        nominalValues = data(:, 4);
        minValues = num2cell(str2num(get(handles.absoluteLowBound, 'String')));
        maxValues = num2cell(str2num(get(handles.absoluteHighBound, 'String')));
        data(:, 2) = minValues;
        data(:, 3) = maxValues;
        set(handles.parameterMinMaxTable, 'data', data);
        handles.minMaxTables{parameter} = get(handles.parameterMinMaxTable, 'data');
    end
end
handles = updateFullParameterTable(handles);
guidata(hObject, handles);


function absoluteLowBound_Callback(hObject, eventdata, handles)
% hObject    handle to absoluteLowBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of absoluteLowBound as text
%        str2double(get(hObject,'String')) returns contents of absoluteLowBound as a double


% --- Executes during object creation, after setting all properties.
function absoluteLowBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absoluteLowBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function absoluteHighBound_Callback(hObject, eventdata, handles)
% hObject    handle to absoluteHighBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of absoluteHighBound as text
%        str2double(get(hObject,'String')) returns contents of absoluteHighBound as a double


% --- Executes during object creation, after setting all properties.
function absoluteHighBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absoluteHighBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
