function varargout = IQMmanualtuning(varargin)
% IQMmanualtuning: Allows to compare simulated experiments data to
% measurements and perform manual parameter tuning. Additionally it is
% possible to run simulations of selected experiments directly from the 
% IQMmanualtuning window, allowing for a deeper analysis of what happens 
% in the system.
%
% USAGE:
% ======
% [project] = IQMmanualtuning(project)
% [project] = IQMmanualtuning(project,modelindex)
% [project] = IQMmanualtuning(project,estimation)
% [project] = IQMmanualtuning(project,modelindex,additionalcomponents)
% [project] =
% IQMmanualtuning(project,estimation,additionalcomponents)
%
% project: IQMprojectSB to tune manually. 
% modelindex: The index of the model in the project to use for tuning.
% estimation: Same structure as used for parameter estimation
%   (IQMparameterestimation), defining the (global) parameters and their
%   bounds, and modelindex.
% additionalcomponents: Cell-array containing the names of components
%   (states or variables) that have not been measured but which should also
%   be plotted in IQMmanualtuning.
%
% DEFAULT VALUES:
% ===============
% modelindex: 1:n, where n is the number of the models in the project. This
%             means that tuning is carried out for all models in the
%             project if modelindex is not given
% estimation: all global parameters and all models, bounds: nominalvalue*[0.01 100]
% additionalcomponents: {} (none)
%
% Output Arguments:
% =================
% The manually tuned project is given back.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS PRO>>>


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IQMmanualtuning_OpeningFcn, ...
                   'gui_OutputFcn',  @IQMmanualtuning_OutputFcn, ...
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
% --- Executes just before IQMmanualtuning is made visible.
function IQMmanualtuning_OpeningFcn(hObject, eventdata, handles, varargin)

% Handle variable input arguments
estimation = [];
if nargin >= 4,
    project = varargin{1};
    if ~isIQMprojectSB(project), 
        error('Input argument is not an IQMprojectSB.'); 
    end
    projectstruct = IQMstruct(project);
    modelindex = [1:length(projectstruct.models)];
end
if nargin >= 5,
    if ~isstruct(varargin{2}),
        modelindex = varargin{2};
    else
        estimation = varargin{2};
        modelindex = 1; % default value
        try modelindex = estimation.modelindex; catch, end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HANDLE PARAM DATA MATRIX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iscell(estimation.parameters),
            if ~isempty(estimation.parameters),
                paramnames = estimation.parameters(:,1);
                paramlowbounds = cell2mat(estimation.parameters(:,2));
                paramhighbounds = cell2mat(estimation.parameters(:,3));
            else
                paramnames = {};
                paramlowbounds = [];
                paramhighbounds = [];
            end
            estimation.parameters = [];
            estimation.parameters.names = paramnames;
            estimation.parameters.lowbounds = paramlowbounds;
            estimation.parameters.highbounds = paramhighbounds;
        end
        if iscell(estimation.parameterslocal),
            if ~isempty(estimation.parameterslocal),
                paramnameslocal = estimation.parameterslocal(:,1);
                paramlowboundslocal = cell2mat(estimation.parameterslocal(:,2));
                paramhighboundslocal = cell2mat(estimation.parameterslocal(:,3));
            else
                paramnameslocal = {};
                paramlowboundslocal = [];
                paramhighboundslocal = [];
            end
            estimation.parameterslocal = [];
            estimation.parameterslocal.names = paramnameslocal;
            estimation.parameterslocal.lowbounds = paramlowboundslocal;
            estimation.parameterslocal.highbounds = paramhighboundslocal;
        end
        if iscell(estimation.initialconditions),
            if ~isempty(estimation.initialconditions),
                icnames = estimation.initialconditions(:,1);
                iclowbounds = cell2mat(estimation.initialconditions(:,2));
                ichighbounds = cell2mat(estimation.initialconditions(:,3));
            else
                icnames = {};
                iclowbounds = [];
                ichighbounds = [];
            end
            estimation.initialconditions = [];
            estimation.initialconditions.names = icnames;
            estimation.initialconditions.lowbounds = iclowbounds;
            estimation.initialconditions.highbounds = ichighbounds;
        end
    end
end
if nargin >= 6,
    additionalPlotNames = varargin{3};
else 
    additionalPlotNames = {};
end
% if nargin<3 || nargin>6,
%     error('Incorrect number of input arguments.');
% end

% get integrator options
try
    OPTIONS = estimation.integrator.options; 
catch
    OPTIONS = [];
end

handles.modelindex = modelindex;
handles.projectstruct = struct(project);
handles.estimation = estimation;
handles.integratoroptions = OPTIONS;
handles.boundsaxis = NaN;
set(handles.bounds,'String','Nominal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine ALL parameters changed by events (in all models and in all experiments)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that are changed by events are not allowed to be changed by the user during tuning
pnames = {};
for k=1:length(handles.projectstruct.models),
    model = handles.projectstruct.models{k};
    ms = struct(model);
    pnames = {pnames{:} ms.parameters.name};
end
pnames = unique(pnames);
apnames = {};
% MODELS FIRST:
for k1=1:length(handles.projectstruct.models),
    model = handles.projectstruct.models{k1};
    ms = struct(model);
    % collect all assignment variables in all events that are parameters
    for k=1:length(ms.events),
        for k2=1:length(ms.events(k).assignment),
            vname = ms.events(k).assignment(k2).variable;
            if ~isempty(strmatchIQM(vname,pnames,'exact')),
                apnames{end+1} = vname;
            end
        end
    end
end
% EXPERIMENTS SECOND:
for k1=1:length(handles.projectstruct.experiments),
    experiment = handles.projectstruct.experiments(k1).experiment;
    es = struct(experiment);
    % collect all assignment variables in all events that are parameters
    for k=1:length(es.stateevents),
        for k2=1:length(es.stateevents(k).assignment),
            vname = es.stateevents(k).assignment(k2).variable;
            if ~isempty(strmatchIQM(vname,pnames,'exact')),
                apnames{end+1} = vname;
            end
        end
    end
end
% apnames are the names of the parameters that are changed by events.
% these are not allowed to be tuned
apnames = unique(apnames);
% Save the info in the 
handles.event_param_names = apnames;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle experimentindices
try 
    experimentindices = estimation.experiments.indices; 
catch
    experimentindices = [1:length(handles.projectstruct.experiments)];
end
handles.nrsliders = 7;  % not really used everywhere yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD EXTRA PLOT NAMES TO THE MEASUREMENT STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additionally save the unchanged project structure for output
handles.saveORIGprojectstruct = handles.projectstruct;
for e=1:length(handles.projectstruct.experiments),
    measurements = handles.projectstruct.experiments(e).measurements;
    for m=1:length(measurements),
        msm = struct(measurements{m});
        for k=1:length(additionalPlotNames),
            % check if additional plot name is already present
            if isempty(strmatchIQM(additionalPlotNames{k},{msm.data(1:end).name},'exact')),
                msm.data(end+1).name = additionalPlotNames{k};
                msm.data(end).notes = 'User-selected additional component for plotting';
                msm.data(end).values = NaN*ones(length(msm.time),1);
                msm.data(end).maxvalues = NaN*ones(length(msm.time),1);
                msm.data(end).minvalues = NaN*ones(length(msm.time),1);
            end
        end
        measurements{m} = IQMmeasurement(msm);
    end
    handles.projectstruct.experiments(e).measurements = measurements;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK WHICH EXPERIMENTS DO NOT HAVE MEASUREMENT DATA ASSIGNED TO
% THESE ARE SKIPPED FROM THE CONSIDERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = '';
useexperimentindices = [];
for e=1:length(experimentindices),
    if isempty(handles.projectstruct.experiments(experimentindices(e)).measurements),
        text = sprintf('%sExperiment %d has no measurements assigned to ... not considered here.\n',text,experimentindices(e));
    else
        useexperimentindices = [useexperimentindices experimentindices(e)];
    end
end
if ~isempty(text),
    disp(text);
end
if isempty(useexperimentindices),
    error('No measurements present in the project. No comparison possible.');
end
handles.experimentindices = useexperimentindices;
[handles.plotdata, handles.work] = initializePlotdata(handles.projectstruct,modelindex,handles);
% set modelselection and choose first model
set(handles.modelselection,'String',{handles.plotdata.model.name});
set(handles.modelselection,'Value',1);
% set experimentselection for first model and first experiment
set(handles.experimentselection,'String',{handles.plotdata.model(1).experiment.name});
set(handles.experimentselection,'Value',1);
% set component selection

% modified by Tõnis Laasfeld 21.11.2018
allComponentNames = handles.plotdata.model(1).allmeascomponents;
residualNames = cell(size(allComponentNames));
for component = 1 : numel(allComponentNames)
   residualNames{component} = [allComponentNames{component}, '_residuals']; 
end
allNewComponentNames = [allComponentNames, residualNames];
set(handles.componentselection,'String',allNewComponentNames);
%set(handles.componentselection,'String',handles.plotdata.model(1).allmeascomponents);
set(handles.componentselection,'Value',[1:length(handles.plotdata.model(1).allmeascomponents)]);
% determine parameter data
handles = getparamdata(handles);
% set all parameter selection boxes with parameter names (model dependent)
initializeParameterLists(handles);
% select plottype
handles.dataPlotType = 'plot';     
% set errorbarflag to 1
handles.errorbars = 1;
% set legendflag to 1
handles.legendflag = 1;
% Initialize export figure handle and grid flag
handles.exportFigureHandle = [];
handles.grid = 0;
% Generate arrays of parameter sliders
handles.parameterMaxValuesArray = cell(handles.nrsliders, 1);
handles.parameterMinValuesArray = cell(handles.nrsliders, 1);
handles.parameterValuesArray = cell(handles.nrsliders, 1);
handles.parameterSlidersArray = cell(handles.nrsliders, 1);
handles.parameterDropdownArray = cell(handles.nrsliders, 1);
for parameterIndex = 1 : handles.nrsliders
   handles.parameterMaxValuesArray{parameterIndex} = eval(['handles.manualmax', num2str(parameterIndex)]);
   handles.parameterMinValuesArray{parameterIndex} = eval(['handles.manualmin', num2str(parameterIndex)]);
   handles.parameterValuesArray {parameterIndex} = eval(['handles.value', num2str(parameterIndex)]);
   handles.parameterSlidersArray{parameterIndex} = eval(['handles.manualslider', num2str(parameterIndex)]);
   handles.parameterDropdownArray{parameterIndex} = eval(['handles.manualparam', num2str(parameterIndex)]);
end
% Doing a first plot
doPlot(handles);
handles=doSimAndPlot(handles);
% Choose default command line output for IQMmanualtuning
handles.output = hObject;
% set up individual weight settings for each measurement
%handles.estimation.experiments.measurementweight = cell(numel(handles.estimation.experiments.weight), 1);
% Update handles structure

guidata(hObject, handles);

uiwait;
return

% --- Outputs from this function are returned to the command line.
function varargout = IQMmanualtuning_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
varargout{1} = handles.output;
% delete temporary MEX models
global compiledExpModelsIQMparamestGUI % if not empty then models are precompiled and should not be deleted
if isempty(compiledExpModelsIQMparamestGUI),
    clear mex
    for m=1:length(handles.plotdata.model)
        for e=1:length(handles.plotdata.model(m).experiment),
            delete(handles.plotdata.model(m).experiment(e).mexfullpath);
        end
    end
end
%varargout{2} = handles.measurementweight;
% close the GUI
delete(hObject);
return

% --- Closed by user
function IQMmanualtuning_CloseRequestFcn(hObject, eventdata, handles)
exit_Callback(hObject, eventdata, handles)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXIT IQMmanualtuning CALLBACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exit_Callback(hObject, eventdata, handles)
    % construct the output (its just an updated IQMprojectSB in which the 
    % manually tuned parametervalues are entered)
    % Modified by Tõnis Laasfeld
    % Here also the experiment weight will be modified with new features of
    % manual tuning
    % FOR THIS USE THE SAVED ORIG PROJECT!
    for k=1:length(handles.modelindex),
        m = handles.modelindex(k);
        % get model to update
        model = handles.saveORIGprojectstruct.models{m};
        % set manually tuned parameters
        % take away the things that are no parameters
        noParamIndices = strmatchIQM('No Parameter',handles.parammodel(k).names,'exact');
        paramIndices = setdiff([1:length(handles.parammodel(k).names)],noParamIndices);
        paramnames = handles.parammodel(k).names(paramIndices);
        paramvalues = handles.parammodel(k).values(paramIndices);
        model = IQMparameters(model,paramnames,paramvalues);
        % add model to project
        handles.saveORIGprojectstruct.models{m} = model;
        %estimationIndex = get(handles.listofestimations
        handles.saveORIGprojectstruct.estimations{1} = handles.estimation;
        %
        for i = 1 : numel(handles.experimentindices)
            expIndex = handles.experimentindices(i);
            measstruct = handles.plotdata.model.experiment(i).measurement;
            meas = struct(handles.saveORIGprojectstruct.experiments(expIndex).measurements{1});
            meas.time = measstruct.timevector';

            for d = 1 : numel(meas.data)
                meas.data(d).values = measstruct.componentvalues(:, d);
                meas.data(d).minvalues = measstruct.minvalues(:, d);
                meas.data(d).maxvalues = measstruct.maxvalues(:, d);
            end

            measurement = IQMmeasurement(meas);
            handles.saveORIGprojectstruct.experiments(expIndex).measurements{1} = measurement;
        end
    end
    handles.output = IQMprojectSB(handles.saveORIGprojectstruct);
    % Update handles structure
    guidata(hObject, handles);
    uiresume;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles=doSimAndPlot(handles)
% get the parameter values to use for simulation
m = get(handles.modelselection,'Value');
%pv = handles.parammodel(m).values;
% run through all experiments and perform them for the parameter settings
residualStack = [];

prob = zeros(numel(handles.plotdata.model(m).experiment), 1);
for e = 1:length(handles.plotdata.model(m).experiment),
    mexmodel = handles.plotdata.model(m).experiment(e).mexmodel;
    timevector = handles.plotdata.model(m).experiment(e).timevector;
    ic = handles.plotdata.model(m).experiment(e).initialconditions;
    % take away the things that are no parameters
    noParamIndices = strmatchIQM('No Parameter',handles.parammodel(m).names,'exact');
    paramIndices = setdiff([1:length(handles.parammodel(m).names)],noParamIndices);
    paramnames = handles.parammodel(m).names(paramIndices);
    paramvalues = handles.parammodel(m).values(paramIndices);
    % construct parameter vector for simulation
    pv = makeparamvecIQM(mexmodel,paramnames,paramvalues);
    try
        tmin = timevector(1); tmax = timevector(end);
        tvec = [tmin:(tmax-tmin)/1000:tmax];
        simdata = feval(mexmodel,tvec,ic,pv,handles.integratoroptions);        
    catch
        simdata = [];
        break;
    end
    stateindices = handles.plotdata.model(m).experiment(e).stateindices;
    statevalues = simdata.statevalues(:,stateindices);
    variableindices = handles.plotdata.model(m).experiment(e).variableindices;
    variablevalues = simdata.variablevalues(:,variableindices);
    handles.plotdata.model(m).experiment(e).componentvalues = [statevalues variablevalues];
    
    % calculate residuals
    handles.work(e).timevector = handles.plotdata.model(m).experiment(e).measurement.timevector;
    handles.work(e).statevalues = statevalues;
    handles.work(e).variablevalues = variablevalues;
    try
    residuals = calculateResiduals(handles.integratoroptions, handles.work(e), ic, pv);
    residValues = fliplr(residuals.measurement.residuals);
    [h, p] = runstest(residValues);
    prob(e) = -log(signtest(residValues) * p);
    numcomponents = numel(handles.plotdata.model.allmeascomponents);
    if numel(handles.plotdata.model(m).experiment(e).measurement.componentnames) < numel(handles.plotdata.model.allmeascomponents)*2
        addNames = cell(1, numcomponents);
        for comp = 1 : numcomponents
            addNames{comp} = [handles.plotdata.model(m).experiment(e).measurement.componentnames{comp}, '_residuals'];
        end
        handles.plotdata.model(m).experiment(e).measurement.componentnames = [handles.plotdata.model(m).experiment(e).measurement.componentnames, addNames];
        handles.plotdata.model(m).experiment(e).measurement.componentvalues = [handles.plotdata.model(m).experiment(e).measurement.componentvalues, residValues];
        handles.plotdata.model(m).experiment(e).measurement.maxvalues = [handles.plotdata.model(m).experiment(e).measurement.maxvalues, nan(size(handles.plotdata.model(m).experiment(e).measurement.maxvalues))]; 
        handles.plotdata.model(m).experiment(e).measurement.minvalues = [handles.plotdata.model(m).experiment(e).measurement.minvalues, nan(size(handles.plotdata.model(m).experiment(e).measurement.minvalues))];

        handles.plotdata.model(m).experiment(e).componentnames = [handles.plotdata.model(m).experiment(e).componentnames, addNames];
        handles.plotdata.model(m).experiment(e).componentvalues = [handles.plotdata.model(m).experiment(e).componentvalues, zeros(size(handles.plotdata.model(m).experiment(e).componentvalues))];
    else
        handles.plotdata.model(m).experiment(e).measurement.componentvalues(:,numcomponents +1 : numcomponents*2) = residValues;
        handles.plotdata.model(m).experiment(e).measurement.maxvalues(:,numcomponents + 1 : numcomponents*2) = nan(size(handles.plotdata.model(m).experiment(e).measurement.maxvalues, 1),numcomponents); 
        handles.plotdata.model(m).experiment(e).measurement.minvalues(:,numcomponents + 1 : numcomponents*2) = nan(size(handles.plotdata.model(m).experiment(e).measurement.minvalues, 1),numcomponents);

        handles.plotdata.model(m).experiment(e).componentvalues(:,numcomponents + 1 : numcomponents*2) = zeros(size(handles.plotdata.model(m).experiment(e).componentvalues, 1),numcomponents);
    end
    catch MException
        ''
    end
end

set(handles.runsTestProbability, 'String', num2str(sum(prob)));


if ~isempty(simdata),
    doPlot(handles);
else
    errordlg('Parameter setting leads to error.');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPlot(handles)
try
    if isequal(get(handles.experimentWise, 'Value'), true)
        doPlotExperimentWise(handles); 
    else
        doPlotComponentWise(handles);
    end
catch MException
   errordlg('This selection is not possible.','Error','on');   
end
return

function doPlotComponentWise(handles)
colorvector = {'b','g','r','c','m','y','k'};
% markervector = {'o','x','+','*','s','d','v','^','<','>','p','h'};
warning off;
% get the data to plot
plotdata = handles.plotdata;
m = get(handles.modelselection,'Value');
eall = get(handles.experimentselection,'Value');
callnames = get(handles.componentselection,'String');
callselected = get(handles.componentselection,'Value');
NTOTAL = length(callselected);
NROW = ceil(sqrt(NTOTAL));
NCOL = ceil(NTOTAL/NROW);

for sp = 1:numel(eall),
    e = eall(sp);
    
    edata = plotdata.model(m).experiment(e);
    % general information
    titletext = regexprep(edata.name,'_',' ');
    xlabeltext = 'Time';
    % simulated data information
    sim_timevector = edata.timevector;
    sim_componentnames = edata.componentnames;
    sim_componentvalues = edata.componentvalues;
    for k=1:length(edata.componentnames),
        if isempty(strmatchIQM(edata.componentnames{k},callnames(callselected),'exact')),
            sim_componentvalues(:,k) = NaN;
        end
    end
    % plot simulated data
    if ~isnan(handles.boundsaxis),
        [sim_timevector,sim_componentvalues] = stripSmallerTimes(sim_timevector,sim_componentvalues,handles.boundsaxis(1));    
    end

    for kplotsim=1:numel(callselected)
        subplot(NROW,NCOL,kplotsim,'Parent',handles.plotpanel);
        feval(handles.dataPlotType,sim_timevector,sim_componentvalues(:,callselected(kplotsim)),'linewidth',2,'color',colorvector{mod(sp-1,7)+1}); hold on;
    end
    
    
    % measured data information
    for meas=1:length(edata.measurement),
        meas_timevector = edata.measurement(meas).timevector;
        meas_componentnames = edata.measurement(meas).componentnames;
        meas_componentvalues = edata.measurement(meas).componentvalues;
        meas_maxvalues = edata.measurement(meas).maxvalues;
        meas_minvalues = edata.measurement(meas).minvalues;
        for k=1:length(edata.measurement(meas).componentnames),
            if isempty(strmatchIQM(edata.measurement(meas).componentnames{k},callnames(callselected),'exact')),
                meas_componentvalues(:,k) = NaN;
                meas_maxvalues(:,k) = NaN;
                meas_minvalues(:,k) = NaN;
            end
        end
        % continue
%         marker = markervector{mod(meas-1,length(markervector))+1};
%         feval(handles.dataPlotType,meas_timevector,meas_componentvalues,['--' marker]); hold on;
        color4meas = colorvector{mod(meas-1,length(colorvector))+1};
        if ~isnan(handles.boundsaxis),
            [meas_timevector,meas_componentvalues] = stripSmallerTimes(meas_timevector,meas_componentvalues,handles.boundsaxis(1));   
        end

        for kplotsim=1:numel(callselected)
            subplot(NROW,NCOL,kplotsim,'Parent',handles.plotpanel);
            feval(handles.dataPlotType,meas_timevector,meas_componentvalues(:,callselected(kplotsim)),'*','linewidth',2,'color',colorvector{mod(sp-1,7)+1}); hold on;
        end

        if handles.errorbars == 1 && strcmp(handles.dataPlotType,'plot'),
            % plot error bounds
            for k=1:numel(callselected),
                subplot(NROW,NCOL,k,'Parent',handles.plotpanel);
                hold on;
                color = colorvector{mod(sp-1,7)+1};
%                 for k1 = 1:size(meas_timevector,1),
                for k1 = 1:length(meas_timevector),
                    if ~isnan(handles.boundsaxis),
                        if meas_timevector(k1) >= handles.boundsaxis(1),
                                        feval(handles.dataPlotType,[meas_timevector(k1),meas_timevector(k1)],[meas_minvalues(k1,k),meas_maxvalues(k1,k)],['.:',color]);
                        end
                    else
                                        feval(handles.dataPlotType,[meas_timevector(k1),meas_timevector(k1)],[meas_minvalues(k1,k),meas_maxvalues(k1,k)],['.:',color]);
                    end
                end
            end
        end
    end


if ~isnan(handles.boundsaxis),
    axis(handles.boundsaxis)
end

end

for k=1:numel(callselected)
    subplot(NROW,NCOL,k,'Parent',handles.plotpanel);
    expNames = cell(numel(eall)*2, 1);
    for sp = 1:numel(eall),
        e = eall(sp);    
        edata = plotdata.model(m).experiment(e);
        expNames{sp*2} = edata.name;
        expNames{sp*2-1} = [edata.name, ' fit'];
    end
    if handles.legendflag == 1,
        hlhlx = legend(expNames);
        set(hlhlx,'Interpreter','none');
    end
    hlhlx = title(sim_componentnames{callselected(k)});
    set(hlhlx,'Interpreter','none');
    hlhlx = xlabel(xlabeltext);
    set(hlhlx,'Interpreter','none');    
    hold off;
end

return

function doPlotExperimentWise(handles)
colorvector = {'b','g','r','c','m','y','k'};
% markervector = {'o','x','+','*','s','d','v','^','<','>','p','h'};
warning off;
% get the data to plot
plotdata = handles.plotdata;
m = get(handles.modelselection,'Value');
eall = get(handles.experimentselection,'Value');
callnames = get(handles.componentselection,'String');
callselected = get(handles.componentselection,'Value');
NTOTAL = length(eall);
NROW = ceil(sqrt(NTOTAL));
NCOL = ceil(NTOTAL/NROW);
for sp = 1:NTOTAL,
    e = eall(sp);
    subplot(NROW,NCOL,sp,'Parent',handles.plotpanel);
    edata = plotdata.model(m).experiment(e);
    % general information
    titletext = regexprep(edata.name,'_',' ');
    xlabeltext = 'Time';
    % simulated data information
    sim_timevector = edata.timevector;
    sim_componentnames = edata.componentnames;
    sim_componentvalues = edata.componentvalues;
    for k=1:length(edata.componentnames),
        if isempty(strmatchIQM(edata.componentnames{k},callnames(callselected),'exact')),
            sim_componentvalues(:,k) = NaN;
        end
    end
    % plot simulated data
if ~isnan(handles.boundsaxis),
    [sim_timevector,sim_componentvalues] = stripSmallerTimes(sim_timevector,sim_componentvalues,handles.boundsaxis(1));    
end

    for kplotsim=1:size(sim_componentvalues,2)
        feval(handles.dataPlotType,sim_timevector,sim_componentvalues(:,kplotsim),'linewidth',2,'color',colorvector{mod(kplotsim-1,7)+1}); hold on;
    end
    
    
    % measured data information
    for meas=1:length(edata.measurement),
        meas_timevector = edata.measurement(meas).timevector;
        meas_componentnames = edata.measurement(meas).componentnames;
        meas_componentvalues = edata.measurement(meas).componentvalues;
        meas_maxvalues = edata.measurement(meas).maxvalues;
        meas_minvalues = edata.measurement(meas).minvalues;
        for k=1:length(edata.measurement(meas).componentnames),
            if isempty(strmatchIQM(edata.measurement(meas).componentnames{k},callnames(callselected),'exact')),
                meas_componentvalues(:,k) = NaN;
                meas_maxvalues(:,k) = NaN;
                meas_minvalues(:,k) = NaN;
            end
        end
        % continue
%         marker = markervector{mod(meas-1,length(markervector))+1};
%         feval(handles.dataPlotType,meas_timevector,meas_componentvalues,['--' marker]); hold on;
        color4meas = colorvector{mod(meas-1,length(colorvector))+1};
if ~isnan(handles.boundsaxis),
    [meas_timevector,meas_componentvalues] = stripSmallerTimes(meas_timevector,meas_componentvalues,handles.boundsaxis(1));   
end

        for kplotsim=1:size(meas_componentvalues,2)
            feval(handles.dataPlotType,meas_timevector,meas_componentvalues(:,kplotsim),'*','linewidth',2,'color',colorvector{mod(kplotsim-1,7)+1}); hold on;
        end

        if handles.errorbars == 1 && strcmp(handles.dataPlotType,'plot'),
            % plot error bounds
            for k=1:length(meas_componentnames),
                color = colorvector{mod(k-1,7)+1};
%                 for k1 = 1:size(meas_timevector,1),
                for k1 = 1:length(meas_timevector),
if ~isnan(handles.boundsaxis),
    if meas_timevector(k1) >= handles.boundsaxis(1),
                    feval(handles.dataPlotType,[meas_timevector(k1),meas_timevector(k1)],[meas_minvalues(k1,k),meas_maxvalues(k1,k)],['.:',color]);
    end
else
                    feval(handles.dataPlotType,[meas_timevector(k1),meas_timevector(k1)],[meas_minvalues(k1,k),meas_maxvalues(k1,k)],['.:',color]);
end
                end
            end
        end
    end
    hold off;
    if handles.legendflag == 1,
        hlhlx = legend(sim_componentnames);
        set(hlhlx,'Interpreter','none');
    end
    hlhlx = title(titletext);
    set(hlhlx,'Interpreter','none');
    hlhlx = xlabel(xlabeltext);
    set(hlhlx,'Interpreter','none');

if ~isnan(handles.boundsaxis),
    axis(handles.boundsaxis)
end

end


return

function [t,y] = stripSmallerTimes(t,y,minT)
    indexdel = find(t<minT);
    t(indexdel) = [];
    y(indexdel,:) = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine parameter data and save in handles structure
% Only parameters that appear in all expmmexmodels are considered.
% Furthermore, they need to have the same values (so that they are not
% modified due to experimental settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = getparamdata(handles)
if isempty(handles.estimation),
    for m = 1:length(handles.plotdata.model),
        % get parameter names that appear in all expmexmodels
        for e = 1:length(handles.plotdata.model(m).experiment),
            model = handles.plotdata.model(m).experiment(e).mexmodel;
            if e == 1,
                parameters = feval(model,'parameters');
            else
                parameters = intersect(parameters, feval(model,'parameters'));
            end
        end
        % keep only those parameters that have the same value (not changed by
        % experiment or changed equally)
        indexdelete = [];
        for e = 1:length(handles.plotdata.model(m).experiment),
            model = handles.plotdata.model(m).experiment(e).mexmodel;
            if e == 1,
                % first run just get the parametervalues
                parametervalues = IQMparameters(model,parameters);
            else
                % subsequent runs ... check if equal and get the indices of the
                % non equal entries
                pv = IQMparameters(model,parameters);
                indexdelete = unique([indexdelete(:)', find(pv(:)'~=parametervalues(:)')]);
            end
        end
        indexkeep = setdiff([1:length(parameters)],indexdelete);
        parameters = parameters(indexkeep);
        % FINALLY ALSO REMOVE THE PARAMETERS THAT ARE CHANGED BY EVENTS DURING
        % THE SIMULATION (DEFINED IN: handles.event_param_names)
        parameters = setdiff(parameters,handles.event_param_names);
        parametervalues = IQMparameters(model,parameters);
        % add param info to structure
        handles.parammodel(m).names = parameters;
        handles.parammodel(m).values = parametervalues;
        handles.parammodel(m).startvalues = parametervalues;
        max = 100*parametervalues;
        max(find(max==0)) = 1;
        handles.parammodel(m).max = max;
        handles.parammodel(m).startmax = max;
        handles.parammodel(m).min = 0.01*parametervalues;
        handles.parammodel(m).startmin = 0.01*parametervalues;
    end
else
    % add param info to structure
    % directly from estimation structure (just one model possible)
    model = handles.plotdata.model(1).experiment(1).mexmodel;
    handles.parammodel.names = handles.estimation.parameters.names;
    parametervalues = IQMparameters(model,handles.parammodel.names);
    handles.parammodel.values = parametervalues;
    handles.parammodel.startvalues = parametervalues;
    if ~isfield(handles.estimation.parameters,'lowbounds'),
        handles.estimation.parameters.lowbounds = [];
    end
    if ~isfield(handles.estimation.parameters,'highbounds'),
        handles.estimation.parameters.highbounds = [];
    end
    if isempty(handles.estimation.parameters.lowbounds),
        min = 0.01*parametervalues;
    else
        min = handles.estimation.parameters.lowbounds;
    end
    if isempty(handles.estimation.parameters.highbounds),
        max = 100*parametervalues;
    else
        max = handles.estimation.parameters.highbounds;
    end
    max(find(max==0)) = 1;
    handles.parammodel.max = max;
    handles.parammodel.startmax = max;
    handles.parammodel.min = min;
    handles.parammodel.startmin = min;
end
% if fewer parameters than handles.nrsliders then add dummy ones
if length(parametervalues) < handles.nrsliders,
    for k=1:(handles.nrsliders-length(parametervalues)),
        handles.parammodel.names{end+1} = 'No Parameter';
        handles.parammodel.values(end+1) = NaN;
        handles.parammodel.startvalues(end+1) = NaN;
        handles.parammodel.max(end+1) = NaN;
        handles.parammodel.startmax(end+1) = NaN;
        handles.parammodel.min(end+1) = NaN;
        handles.parammodel.startmin(end+1) = NaN;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize parameter lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = initializeParameterLists(handles,varargin)
if nargin == 1,
    i1 = 1;
    i2 = 2;
    i3 = 3;
    i4 = 4;
    i5 = 5;
    i6 = 6;
    i7 = 7;
elseif nargin == 2,
    if strcmp(varargin{1},'resetall'),
        i1 = 1;
        i2 = 2;
        i3 = 3;
        i4 = 4;
        i5 = 5;
        i6 = 6;
        i7 = 7;
    else
        i1 = get(handles.manualparam1,'Value');
        i2 = get(handles.manualparam2,'Value');
        i3 = get(handles.manualparam3,'Value');
        i4 = get(handles.manualparam4,'Value');
        i5 = get(handles.manualparam5,'Value');
        i6 = get(handles.manualparam6,'Value');
        i7 = get(handles.manualparam7,'Value');
    end
end

m = get(handles.modelselection,'Value');
% get the models parameters and values
parameters = handles.parammodel(m).names;
parametervalues = handles.parammodel(m).values;
parammax = handles.parammodel(m).max;
parammin = handles.parammodel(m).min;
% set the 7 parameter lists (adjust list of parameters)
nrtakeaway = length(strmatchIQM('No Parameter',parameters,'exact'));
if nrtakeaway > 0,
    endindex = handles.nrsliders-nrtakeaway;
    set(handles.manualparam1,'String',parameters(1:endindex));
    set(handles.manualparam2,'String',parameters(1:endindex));
    set(handles.manualparam3,'String',parameters(1:endindex));
    set(handles.manualparam4,'String',parameters(1:endindex));
    set(handles.manualparam5,'String',parameters(1:endindex));
    set(handles.manualparam6,'String',parameters(1:endindex));
    set(handles.manualparam7,'String',parameters(1:endindex));
else
    set(handles.manualparam1,'String',parameters);
    set(handles.manualparam2,'String',parameters);
    set(handles.manualparam3,'String',parameters);
    set(handles.manualparam4,'String',parameters);
    set(handles.manualparam5,'String',parameters);
    set(handles.manualparam6,'String',parameters);
    set(handles.manualparam7,'String',parameters);
end    
% set selected values
set(handles.manualparam1,'Value',i1);
set(handles.manualparam2,'Value',i2);
set(handles.manualparam3,'Value',i3);
set(handles.manualparam4,'Value',i4);
set(handles.manualparam5,'Value',i5);
set(handles.manualparam6,'Value',i6);
set(handles.manualparam7,'Value',i7);
% set current parameter values
set(handles.value1,'String',parametervalues(i1));
set(handles.value2,'String',parametervalues(i2));
set(handles.value3,'String',parametervalues(i3));
set(handles.value4,'String',parametervalues(i4));
set(handles.value5,'String',parametervalues(i5));
set(handles.value6,'String',parametervalues(i6));
set(handles.value7,'String',parametervalues(i7));
% set max and min values per default to *100 / *0.01
set(handles.manualmax1,'String',parammax(i1));
set(handles.manualmax2,'String',parammax(i2));
set(handles.manualmax3,'String',parammax(i3));
set(handles.manualmax4,'String',parammax(i4));
set(handles.manualmax5,'String',parammax(i5));
set(handles.manualmax6,'String',parammax(i6));
set(handles.manualmax7,'String',parammax(i7));
set(handles.manualmin1,'String',parammin(i1));
set(handles.manualmin2,'String',parammin(i2));
set(handles.manualmin3,'String',parammin(i3));
set(handles.manualmin4,'String',parammin(i4));
set(handles.manualmin5,'String',parammin(i5));
set(handles.manualmin6,'String',parammin(i6));
set(handles.manualmin7,'String',parammin(i7));
% set slider min max and value
set(handles.manualslider1,'Max',101);
set(handles.manualslider2,'Max',101);
set(handles.manualslider3,'Max',101);
set(handles.manualslider4,'Max',101);
set(handles.manualslider5,'Max',101);
set(handles.manualslider6,'Max',101);
set(handles.manualslider7,'Max',101);
set(handles.manualslider1,'Min',1);
set(handles.manualslider2,'Min',1);
set(handles.manualslider3,'Min',1);
set(handles.manualslider4,'Min',1);
set(handles.manualslider5,'Min',1);
set(handles.manualslider6,'Min',1);
set(handles.manualslider7,'Min',1);
% construct vectors
if parammin(i1) > 0,
    vector1 = logspace(log(parammin(i1))/log(10),log(parammax(i1))/log(10),101);
else
    vector1 = [parammin(i1):(parammax(i1)-parammin(i1))/100:parammax(i1)];
end
if parammin(i2) > 0,
    vector2 = logspace(log(parammin(i2))/log(10),log(parammax(i2))/log(10),101);
else
    vector2 = [parammin(i2):(parammax(i2)-parammin(i2))/100:parammax(i2)];
end
if parammin(i3) > 0,
    vector3 = logspace(log(parammin(i3))/log(10),log(parammax(i3))/log(10),101);
else
    vector3 = [parammin(i3):(parammax(i3)-parammin(i3))/100:parammax(i3)];
end
if parammin(i4) > 0,
    vector4 = logspace(log(parammin(i4))/log(10),log(parammax(i4))/log(10),101);
else
    vector4 = [parammin(i4):(parammax(i4)-parammin(i4))/100:parammax(i4)];
end
if parammin(i5) > 0,
    vector5 = logspace(log(parammin(i5))/log(10),log(parammax(i5))/log(10),101);
else
    vector5 = [parammin(i5):(parammax(i5)-parammin(i5))/100:parammax(i5)];
end
if parammin(i6) > 0,
    vector6 = logspace(log(parammin(i6))/log(10),log(parammax(i6))/log(10),101);
else
    vector6 = [parammin(i6):(parammax(i6)-parammin(i6))/100:parammax(i6)];
end
if parammin(i7) > 0,
    vector7 = logspace(log(parammin(i7))/log(10),log(parammax(i7))/log(10),101);
else
    vector7 = [parammin(i7):(parammax(i7)-parammin(i7))/100:parammax(i7)];
end
% set sliders
[dummy,index] = min(abs(vector1-parametervalues(i1)));
set(handles.manualslider1,'Value',index);
[dummy,index] = min(abs(vector2-parametervalues(i2)));
set(handles.manualslider2,'Value',index);
[dummy,index] = min(abs(vector3-parametervalues(i3)));
set(handles.manualslider3,'Value',index);
[dummy,index] = min(abs(vector4-parametervalues(i4)));
set(handles.manualslider4,'Value',index);
[dummy,index] = min(abs(vector5-parametervalues(i5)));
set(handles.manualslider5,'Value',index);
[dummy,index] = min(abs(vector6-parametervalues(i6)));
set(handles.manualslider6,'Value',index);
[dummy,index] = min(abs(vector7-parametervalues(i7)));
set(handles.manualslider7,'Value',index);
% set slider steps
set(handles.manualslider1,'SliderStep',[0.01 0.01]);
set(handles.manualslider2,'SliderStep',[0.01 0.01]);
set(handles.manualslider3,'SliderStep',[0.01 0.01]);
set(handles.manualslider4,'SliderStep',[0.01 0.01]);
set(handles.manualslider5,'SliderStep',[0.01 0.01]);
set(handles.manualslider6,'SliderStep',[0.01 0.01]);
set(handles.manualslider7,'SliderStep',[0.01 0.01]);
% set Enable to off for the elements
nrtakeaway = length(strmatchIQM('No Parameter',parameters,'exact'));
if nrtakeaway >= 7,
    set(handles.manualslider1,'Enable','off');
    set(handles.manualmax1,'Enable','off');
    set(handles.manualmin1,'Enable','off');
    set(handles.manualparam1,'Enable','off');
end
if nrtakeaway >= 6,
    set(handles.manualslider2,'Enable','off');
    set(handles.manualmax2,'Enable','off');
    set(handles.manualmin2,'Enable','off');
    set(handles.manualparam2,'Enable','off');
end
if nrtakeaway >= 5,
    set(handles.manualslider3,'Enable','off');
    set(handles.manualmax3,'Enable','off');
    set(handles.manualmin3,'Enable','off');
    set(handles.manualparam3,'Enable','off');
end
if nrtakeaway >= 4,
    set(handles.manualslider4,'Enable','off');
    set(handles.manualmax4,'Enable','off');
    set(handles.manualmin4,'Enable','off');
    set(handles.manualparam4,'Enable','off');
end
if nrtakeaway >= 3,
    set(handles.manualslider5,'Enable','off');
    set(handles.manualmax5,'Enable','off');
    set(handles.manualmin5,'Enable','off');
    set(handles.manualparam5,'Enable','off');
end
if nrtakeaway >= 2,
    set(handles.manualslider6,'Enable','off');
    set(handles.manualmax6,'Enable','off');
    set(handles.manualmin6,'Enable','off');
    set(handles.manualparam6,'Enable','off');
end
if nrtakeaway >= 1,
    set(handles.manualslider7,'Enable','off');
    set(handles.manualmax7,'Enable','off');
    set(handles.manualmin7,'Enable','off');
    set(handles.manualparam7,'Enable','off');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider1 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function sliderGeneralCallback(hObject, eventdata, handles, index)

% get slider settings
v = ceil(get(hObject,'Value'));                       %%%
% get min max settings
maxV = str2double(get(handles.parameterMaxValuesArray{index},'String'));                %%%
minV = str2double(get(handles.parameterMinValuesArray{index},'String'));                %%%
% construct vectors
if minV > 0,
    vector = logspace(log(minV)/log(10),log(maxV)/log(10),101);
else
    vector = [minV:(maxV-minV)/100:maxV];
end
% get paramvalues
value = vector(v);
% set parameter value
set(handles.parameterValuesArray{index},'String',value);                                 %%%
% update parameter information structure in handles with new values
m = get(handles.modelselection,'Value');
pindex = get(handles.parameterDropdownArray{index},'Value');                         %%%
handles.parammodel(m).values(pindex) = value;
% sim and plot
handles = doSimAndPlot(handles);
% Update handles structure
guidata(hObject, handles);

return

function slider1_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider2 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider2_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider3 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider3_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 3);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider4 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider4_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 4);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider5 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider5_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 5);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider6 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider6_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 6);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider7 handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider7_Callback(hObject, eventdata, handles)
sliderGeneralCallback(hObject, eventdata, handles, 7);
return

function valueGeneralCallback(hObject, eventdata, handles, index)
%handles.parameterMaxValuesArray = cell(handles.nrsliders, 1);
%handles.parameterMinValuesArray = cell(handles.nrsliders, 1);
%handles.parameterValuesArray = cell(handles.nrsliders, 1);
%handles.parameterSlidersArray = cell(handles.nrsliders, 1);
%handles.parameterDropdownArray = cell(handles.nrsliders, 1);

% get slider settings
v = ceil(get(hObject,'Value'));                       %%%
% get min max settings
max1 = str2double(get(handles.parameterMaxValuesArray{index},'String'));                %%%
min1 = str2double(get(handles.parameterMinValuesArray{index},'String'));                %%%

% get paramvalues
value = str2double(get(handles.parameterValuesArray{index}, 'String'));

if value < min1, min1 = value; end
if value > max1, max1 = value; end

set(handles.parameterMaxValuesArray{index},'String', num2str(max1));                %%%
set(handles.parameterMinValuesArray{index},'String', num2str(min1));             %%%
% set slider value to correspond to value
if min1 > 0,
    vector1 = logspace(log(min1)/log(10),log(max1)/log(10),101);
else
    vector1 = [min1:(max1-min1)/100:max1];
end
[dummy,sliderValue] = min(abs(vector1-value));
set(handles.parameterSlidersArray{index},'Value',sliderValue);

% update parameter information structure in handles with new values
m = get(handles.modelselection,'Value');
pindex = get(handles.parameterDropdownArray{index},'Value');                         %%%
handles.parammodel(m).values(pindex) = value;
handles.parammodel(m).max(pindex) = max1;
handles.parammodel(m).min(pindex) = min1;
% sim and plot
handles = doSimAndPlot(handles);
% Update handles structure
guidata(hObject, handles);

return

function value1_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 1);
return

function value2_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 2);
return

function value3_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 3);
return

function value4_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 4);
return

function value5_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 5);
return

function value6_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 6);
return

function value7_Callback(hObject, eventdata, handles)
valueGeneralCallback(hObject, eventdata, handles, 7);
return



% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Slider handling function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function slider_Callback(hObject, eventdata, handles)
% % get slider settings
% v1 = ceil(get(handles.manualslider1,'Value'));
% v2 = ceil(get(handles.manualslider2,'Value'));
% v3 = ceil(get(handles.manualslider3,'Value'));
% v4 = ceil(get(handles.manualslider4,'Value'));
% v5 = ceil(get(handles.manualslider5,'Value'));
% v6 = ceil(get(handles.manualslider6,'Value'));
% v7 = ceil(get(handles.manualslider7,'Value'));
% % get min max settings
% max1 = str2double(get(handles.manualmax1,'String'));
% max2 = str2double(get(handles.manualmax2,'String'));
% max3 = str2double(get(handles.manualmax3,'String'));
% max4 = str2double(get(handles.manualmax4,'String'));
% max5 = str2double(get(handles.manualmax5,'String'));
% max6 = str2double(get(handles.manualmax6,'String'));
% max7 = str2double(get(handles.manualmax7,'String'));
% min1 = str2double(get(handles.manualmin1,'String'));
% min2 = str2double(get(handles.manualmin2,'String'));
% min3 = str2double(get(handles.manualmin3,'String'));
% min4 = str2double(get(handles.manualmin4,'String'));
% min5 = str2double(get(handles.manualmin5,'String'));
% min6 = str2double(get(handles.manualmin6,'String'));
% min7 = str2double(get(handles.manualmin7,'String'));
% % construct vectors
% if min1 > 0,
%     vector1 = logspace(log(min1)/log(10),log(max1)/log(10),101);
% else
%     vector1 = [min1:(max1-min1)/100:max1];
% end
% if min2 > 0,
%     vector2 = logspace(log(min2)/log(10),log(max2)/log(10),101);
% else
%     vector2 = [min2:(max2-min2)/100:max2];
% end
% if min3 > 0,
%     vector3 = logspace(log(min3)/log(10),log(max3)/log(10),101);
% else
%     vector3 = [min3:(max3-min3)/100:max3];
% end
% if min4 > 0,
%     vector4 = logspace(log(min4)/log(10),log(max4)/log(10),101);
% else
%     vector4 = [min4:(max4-min4)/100:max4];
% end
% if min5 > 0,
%     vector5 = logspace(log(min5)/log(10),log(max5)/log(10),101);
% else
%     vector5 = [min5:(max5-min5)/100:max5];
% end
% if min6 > 0,
%     vector6 = logspace(log(min6)/log(10),log(max6)/log(10),101);
% else
%     vector6 = [min6:(max6-min6)/100:max6];
% end
% if min7 > 0,
%     vector7 = logspace(log(min7)/log(10),log(max7)/log(10),101);
% else
%     vector7 = [min7:(max7-min7)/100:max7];
% end
% % get paramvalues
% values = [vector1(v1) vector2(v2) vector3(v3) vector4(v4) vector5(v5) vector6(v6) vector7(v7)];
% % set parameter values
% set(handles.value1,'String',values(1));
% set(handles.value2,'String',values(2));
% set(handles.value3,'String',values(3));
% set(handles.value4,'String',values(4));
% set(handles.value5,'String',values(5));
% set(handles.value6,'String',values(6));
% set(handles.value7,'String',values(7));
% % update parameter information structure in handles with new values
% m = get(handles.modelselection,'Value');
% pindex = [];
% pindex(end+1) = get(handles.manualparam1,'Value');
% pindex(end+1) = get(handles.manualparam2,'Value');
% pindex(end+1) = get(handles.manualparam3,'Value');
% pindex(end+1) = get(handles.manualparam4,'Value');
% pindex(end+1) = get(handles.manualparam5,'Value');
% pindex(end+1) = get(handles.manualparam6,'Value');
% pindex(end+1) = get(handles.manualparam7,'Value');
% handles.parammodel(m).values(pindex) = values;
% % sim and plot
% handles = doSimAndPlot(handles);
% % Update handles structure
% guidata(hObject, handles);
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minmax handling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function minmaxGeneralCallback(hObject, eventdata, handles, index)

% get max min values
max1 = str2double(get(handles.parameterMaxValuesArray{index},'String'));                %%%
min1 = str2double(get(handles.parameterMinValuesArray{index},'String'));                %%%
% get slider setting
v1 = ceil(get(handles.parameterSlidersArray{index},'Value'));
% get current value
value1 = str2double(get(handles.parameterValuesArray{index},'String'));
% check and adjust bounds if needed
if value1 < min1, value1 = min1; end
if value1 > max1, value1 = max1; end
% set slider value to correspond to value
if min1 > 0,
    vector1 = logspace(log(min1)/log(10),log(max1)/log(10),101);
else
    vector1 = [min1:(max1-min1)/100:max1];
end
[dummy,sliderValue] = min(abs(vector1-value1));
set(handles.parameterSlidersArray{index},'Value',sliderValue);
% update values field
value1 = vector1(sliderValue);
set(handles.parameterValuesArray{index},'String',value1);
% update parameter information structure in handles with new value
m = get(handles.modelselection,'Value');
pindex = get(handles.parameterDropdownArray{index},'Value');
handles.parammodel(m).max(pindex) = max1;
handles.parammodel(m).min(pindex) = min1;
handles.parammodel(m).values(pindex) = value1;
% sim and plot
handles = doSimAndPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

function minmax1_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 1);
return

function minmax2_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 2);
return

function minmax3_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 3);
return

function minmax4_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 4);
return

function minmax5_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 5);
return

function minmax6_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 6);
return

function minmax7_Callback(hObject, eventdata, handles)
minmaxGeneralCallback(hObject, eventdata, handles, 7);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETER HANDLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function manualParamGeneralCallback(hObject, eventdata, handles, index)

m = get(handles.modelselection,'Value');
pindex = get(handles.parameterDropdownArray{index},'Value');
paramvalue = handles.parammodel(m).values(pindex);
% set current parameter values
set(handles.parameterValuesArray{index},'String',paramvalue);
% set max and min values 
max = handles.parammodel(m).max(pindex);
min = handles.parammodel(m).min(pindex);
set(handles.parameterMaxValuesArray{index},'String',max);
set(handles.parameterMinValuesArray{index},'String',min);
% set slider min max and value
set(handles.parameterSlidersArray{index},'Max',101);
set(handles.parameterSlidersArray{index},'Min',1);
% set slider value to correspond to value
if min > 0,
    vector = logspace(log(min)/log(10),log(max)/log(10),101);
else
    vector = [min:(max-min)/100:max];
end
[dummy,pindex] = minfunct(abs(vector-paramvalue));
set(handles.parameterSlidersArray{index},'Value',pindex);
% set slider steps
set(handles.parameterSlidersArray{index},'SliderStep',[0.01 0.01]);
% Update handles structure
guidata(hObject, handles);
return

function manualparam1_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 1);
return

function [a,b] = minfunct(X)
[a,b] = min(X);
return

function manualparam2_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 2);
return

function manualparam3_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 3);
return

function manualparam4_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 4);
return

function manualparam5_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 5);
return

function manualparam6_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 6);
return

function manualparam7_Callback(hObject, eventdata, handles)
manualParamGeneralCallback(hObject, eventdata, handles, 7);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE PLOTDATA (includes MEXfile data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plotdata, work] = initializePlotdata(projectstruct,modelindex,handles)
% GET ALL MEASUREMENT INFORMATION
infostruct = [];
for k=1:length(projectstruct.models(modelindex)),
    model = projectstruct.models{modelindex(k)};
    experiments = projectstruct.experiments(handles.experimentindices);
    displayFlag = 0;
    expmeasinfo = getexpmeasinfoIQM(model,modelindex(k),experiments,handles.experimentindices,displayFlag);
    infostruct(k).modelstruct = IQMstruct(projectstruct.models{modelindex(k)});
    infostruct(k).modelindex = modelindex(k);
    infostruct(k).expinfostruct = expmeasinfo;
end

work = expmeasinfo;


% ADD SIMULATION DATA TO THE STRUCTURE
plotdata = [];
plotdata.project = projectstruct.name;
plotdata.notes = projectstruct.notes;
plotdata.model = [];
% run through all models
for m=1:length(projectstruct.models(modelindex)),
    % model data
    modelstruct = infostruct(m).modelstruct;
    plotdata.model(m).name = modelstruct.name;      
    plotdata.model(m).notes = modelstruct.notes;
    allmeasuredcomponents = {};
    for e=1:length(infostruct(m).expinfostruct),
        % experiment data
        plotdata.model(m).experiment(e).name = infostruct(m).expinfostruct(e).experimentname;
        plotdata.model(m).experiment(e).mexmodel = infostruct(m).expinfostruct(e).model;
        plotdata.model(m).experiment(e).mexfullpath = infostruct(m).expinfostruct(e).mexfullpath;
        timevector = infostruct(m).expinfostruct(e).timevector;
        timestart = timevector(1);
        timeend = timevector(end);
        timevectorsim = [timestart:(timeend-timestart)/1000:timeend];
        plotdata.model(m).experiment(e).timevector = timevectorsim;
        expstatenames = infostruct(m).expinfostruct(e).statenames;          
        expvariablenames = infostruct(m).expinfostruct(e).variablenames; 
        plotdata.model(m).experiment(e).componentnames = {expstatenames{:} expvariablenames{:}};  
        % simulate to get the state and variable values
        mexmodel = infostruct(m).expinfostruct(e).model;
        ic = infostruct(m).expinfostruct(e).initialconditions;
        plotdata.model(m).experiment(e).initialconditions = ic;
        try
            simdata = feval(mexmodel,timevectorsim,ic,[],handles.integratoroptions);
        catch
            simdata.statevalues = NaN(length(timevectorsim),length(ic));
            simdata.variablevalues = NaN(length(timevectorsim),length(ic));
%            disp(sprintf('Integrator problems for experiment %d.',e));
        end
        % collect all states that are measured
        stateindices = infostruct(m).expinfostruct(e).stateindices;
        variableindices = infostruct(m).expinfostruct(e).variableindices;
        plotdata.model(m).experiment(e).stateindices = stateindices;
        plotdata.model(m).experiment(e).variableindices = variableindices;
        statevalues = simdata.statevalues(:,stateindices);
        variablevalues = simdata.variablevalues(:,variableindices);
        % add simulated state trajectories
        plotdata.model(m).experiment(e).componentvalues = [statevalues variablevalues];
        for meas=1:length(infostruct(m).expinfostruct(e).measurement),
            % measurement data
            plotdata.model(m).experiment(e).measurement(meas).name = infostruct(m).expinfostruct(e).measurement(meas).name;
            timevectormeas = timevector(infostruct(m).expinfostruct(e).measurement(meas).timevectorindices);
            plotdata.model(m).experiment(e).measurement(meas).timevector = timevectormeas;
            % reorder the measurements
            measstatenames = infostruct(m).expinfostruct(e).measurement(meas).statenames;
            measvariablenames = infostruct(m).expinfostruct(e).measurement(meas).variablenames;
            % states 
            for k=1:length(expstatenames),
                index = strmatchIQM(expstatenames{k},measstatenames,'exact');
                if ~isempty(index),
                    plotdata.model(m).experiment(e).measurement(meas).componentnames{k} = measstatenames{index};
                    plotdata.model(m).experiment(e).measurement(meas).componentvalues(:,k) = infostruct(m).expinfostruct(e).measurement(meas).statereferences(:,index);
                    plotdata.model(m).experiment(e).measurement(meas).maxvalues(:,k) = infostruct(m).expinfostruct(e).measurement(meas).statemaxvalues(:,index);
                    plotdata.model(m).experiment(e).measurement(meas).minvalues(:,k) = infostruct(m).expinfostruct(e).measurement(meas).stateminvalues(:,index);
                else
                    plotdata.model(m).experiment(e).measurement(meas).componentnames{k} = 'not available';
                    plotdata.model(m).experiment(e).measurement(meas).componnentvalues(:,k) = NaN(length(timevectormeas),1);
                    plotdata.model(m).experiment(e).measurement(meas).maxvalues(:,k) = NaN(length(timevectormeas),1);
                    plotdata.model(m).experiment(e).measurement(meas).minvalues(:,k) = NaN(length(timevectormeas),1);
                end                    
            end
            offset = length(expstatenames);
            % variables
            for k=1:length(expvariablenames),
                index = strmatchIQM(expvariablenames{k},measvariablenames,'exact');
                if ~isempty(index),
                    plotdata.model(m).experiment(e).measurement(meas).componentnames{offset+k} = measvariablenames{index};
                    plotdata.model(m).experiment(e).measurement(meas).componentvalues(:,offset+k) = infostruct(m).expinfostruct(e).measurement(meas).variablereferences(:,index);
                    plotdata.model(m).experiment(e).measurement(meas).maxvalues(:,offset+k) = infostruct(m).expinfostruct(e).measurement(meas).variablemaxvalues(:,index);
                    plotdata.model(m).experiment(e).measurement(meas).minvalues(:,offset+k) = infostruct(m).expinfostruct(e).measurement(meas).variableminvalues(:,index);
                else
                    plotdata.model(m).experiment(e).measurement(meas).componentnames{offset+k} = 'not available';
                    plotdata.model(m).experiment(e).measurement(meas).componentvalues(:,offset+k) = NaN(length(timevectormeas),1);
                    plotdata.model(m).experiment(e).measurement(meas).maxvalues(:,offset+k) = NaN(length(timevectormeas),1);
                    plotdata.model(m).experiment(e).measurement(meas).minvalues(:,offset+k) = NaN(length(timevectormeas),1);
                end                    
            end      
            allmeasuredcomponents = {allmeasuredcomponents{:} plotdata.model(m).experiment(e).measurement(meas).componentnames{:}};
        end
    end
    allmeasuredcomponents = unique(allmeasuredcomponents);
    plotdata.model(m).allmeascomponents = allmeasuredcomponents;
end 

guidata(handles.SBplot, handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIVIAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in zoombutton.
function legendbutton_Callback(hObject, eventdata, handles)
% toogle the legends in the figure
if handles.legendflag == 0,
    handles.legendflag = 1;
else
    handles.legendflag = 0;
end
% plot
doPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

function residuals = calculateResiduals(OPTIONS, work, ic, pv)

try
    alpha = 0.05;
    modelstruct = IQMstruct(work.IQMmodel);
    for f = 1 : numel(modelstruct.parameters)
        fields{f} = modelstruct.parameters(f).name;
    end
    
    for exp=1:length(work),
         mexmodel = work(exp).model;
         timevector = work(exp).timevector;
         simdata = IQMPsimulate(mexmodel, timevector,ic,fields,pv,OPTIONS);
         work(exp).statevalues = simdata.statevalues;
         work(exp).variablevalues = simdata.variablevalues;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the residuals for each experiment and each measurement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    residuals = [];
    for exp=1:length(work),
        experiment = work(exp);
        residuals(exp).name = ['Residuals for: ' experiment.experimentname];
        AllRxx = [];
        for meas=1:length(experiment.measurement),
            measurement = experiment.measurement(meas);
            timevectorindices = measurement.timevectorindices;
            timevector = experiment.timevector(timevectorindices);
            timescaling = measurement.timescaling(:);
            % states
            stateindices = measurement.stateindices;
            if ~isempty(stateindices),
                simstates = experiment.statevalues(timevectorindices,stateindices);
                stateresiduals = measurement.statereferences - simstates;
                statenames = measurement.statenames;
            else
                stateresiduals = [];
                statenames = {};
            end
            % variables
            variableindices = experiment.measurement(meas).variableindices;
            if ~isempty(variableindices),
                simvariables = experiment.variablevalues(timevectorindices,variableindices);
                variableresiduals = measurement.variablereferences - simvariables;
                variablenames = measurement.variablenames;
            else
                variableresiduals = [];
                variablenames = {};
            end
            values = [stateresiduals variableresiduals];
            % apply the time scaling  
            values = values.*timescaling(:,ones(1,size(values,2)));
            % determine autocorrelations of residuals
            allRxx = [];
            allLag = [];
            indexNaNlater = [];
            for k=1:size(values,2),
                % get values to resample and analyze
                x1 = values(:,k);
                t1 = timevector;
                % remove NaNs
                nanindex = find(isnan(x1));
                x1(nanindex) = [];
                t1(nanindex) = [];
                % determine the sampling time (half the min size for finer resolution)
                deltaT = min(t1(2:end)-t1(1:end-1))/2;
                if ~isempty(x1),
                    [x2,t2] = resampleIQM(t1,x1,deltaT,'linear');
                    x2 = x2-mean(x2); % make mean free
                    [Rxx,lag] = xcorrIQM(x2,x2,length(x2)-1,'coeff');
                else
                    lag = [];
                    Rxx = [];
                    indexNaNlater = [indexNaNlater k];
                end
                allLag = [allLag lag(:)];
                allRxx = [allRxx Rxx(:)];
            end
            % handle indexNaNlater (insert columns of NaN values)
            ncoltot = size(allRxx,2) + length(indexNaNlater);
            distribindices = setdiff([1:ncoltot],indexNaNlater);
            AllLag = allLag(:,1);
            AllRxx(:,distribindices) = allRxx;
            AllRxx(:,indexNaNlater) = NaN*ones(size(AllLag,1),length(indexNaNlater));
            % add to residuals
            residuals(exp).measurement(meas).name = ['Residuals for: ' measurement.name];
            residuals(exp).measurement(meas).components = {statenames{:} variablenames{:}};
            residuals(exp).measurement(meas).time = timevector;        
            residuals(exp).measurement(meas).timescalingFlag = 0;
            residuals(exp).measurement(meas).residuals = values;
            % do a statistical test to check if residuals are white
            H = []; pValue = [];
            for k=1:size(values,2),
                vector = values(:,k);
                vector(find(isnan(vector))) = [];
                if length(vector) >= 3,
                    [H(k), pValue(k)] = swtestIQM(vector,alpha,0);
                else
                    H(k) = NaN;
                    pValue(k) = NaN;
                end
            end
            residuals(exp).measurement(meas).pValue = pValue;
            residuals(exp).measurement(meas).alpha = alpha;
            residuals(exp).measurement(meas).H = H;
            % add also the autocorrelation results
            residuals(exp).measurement(meas).autocorrelation_Rxx = AllRxx;
            residuals(exp).measurement(meas).autocorrelation_lag = AllLag;
        end
    end
catch Mexception
    ''
end
return

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
% set all parameter selection boxes with parameter names (model dependent)
initializeParameterLists(handles);
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in resetall.
function resetall_Callback(hObject, eventdata, handles)
m = get(handles.modelselection,'Value');
% reset all values for current model to starting values
handles.parammodel(m).values = handles.parammodel(m).startvalues;
handles.parammodel(m).max = handles.parammodel(m).startmax;
handles.parammodel(m).min = handles.parammodel(m).startmin;
% reinitialize parameter list
initializeParameterLists(handles,'resetall');
% simulate and plot it
handles = doSimAndPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in resetcurrent.
function resetcurrent_Callback(hObject, eventdata, handles)
m = get(handles.modelselection,'Value');
% reset all values for current model and current parameter selections to starting values
% get current indices
pindices = [];
pindices(end+1) = get(handles.manualparam1,'Value');
pindices(end+1) = get(handles.manualparam2,'Value');
pindices(end+1) = get(handles.manualparam3,'Value');
pindices(end+1) = get(handles.manualparam4,'Value');
pindices(end+1) = get(handles.manualparam5,'Value');
pindices(end+1) = get(handles.manualparam6,'Value');
pindices(end+1) = get(handles.manualparam7,'Value');
% reset values for current parameters
handles.parammodel(m).values(pindices) = handles.parammodel(m).startvalues(pindices);
handles.parammodel(m).max(pindices) = handles.parammodel(m).startmax(pindices);
handles.parammodel(m).min(pindices) = handles.parammodel(m).startmin(pindices);
% reinitialize parameter list
initializeParameterLists(handles,'resetcurrent');
% simulate and plot it
handles = doSimAndPlot(handles);
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in experimentselection.
function experimentselection_Callback(hObject, eventdata, handles)
try
    doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
% check if only one experiment selected ... then show the simulate button
% otherwise hide it
eall = get(handles.experimentselection,'Value');
if length(eall) > 1,
    set(handles.simexpbutton,'Enable','off');
    set(handles.simexpbutton,'BackgroundColor',[0.8 0.8 0.8]);
else
    set(handles.simexpbutton,'Enable','on');
    set(handles.simexpbutton,'BackgroundColor',[1 0.745 0.49]);
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on selection change in experimentselection.
function componentselection_Callback(hObject, eventdata, handles)
try
    doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
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

% --- Simulate current model with current experiment and display in IQMplot window
function simulatebutton_Callback(hObject, eventdata, handles)
m = get(handles.modelselection,'Value');
e = get(handles.experimentselection,'Value');
name = handles.plotdata.model(m).experiment(e).name;

    mexmodel = handles.plotdata.model(m).experiment(e).mexmodel;
    timevector = handles.plotdata.model(m).experiment(e).timevector;
    ic = handles.plotdata.model(m).experiment(e).initialconditions;
    % take away the things that are no parameters
    noParamIndices = strmatchIQM('No Parameter',handles.parammodel(m).names,'exact');
    paramIndices = setdiff([1:length(handles.parammodel(m).names)],noParamIndices);
    paramnames = handles.parammodel(m).names(paramIndices);
    paramvalues = handles.parammodel(m).values(paramIndices);

simdata = IQMPsimulate(mexmodel,timevector,ic,paramnames,paramvalues,handles.integratoroptions);
% remove '_' from name
name = strrep(name,'_',' ');
plotdatastruct = createdatastruct2IQMplotIQM(simdata,name);
IQMplot(plotdatastruct);
% Update handles structure
guidata(hObject, handles);
return

% --- Set axis bounds
function bounds_Callback(hObject, eventdata, handles)
try 
    handles.boundsaxis = eval(get(handles.bounds,'String'));
catch
    handles.boundsaxis = NaN;
    set(handles.bounds,'String','Nominal');
end
if length(handles.boundsaxis) ~= 4,
    handles.boundsaxis = NaN;
    set(handles.bounds,'String','Nominal');
end    
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

function plottingStyle_Callback(hObject, eventdata, handles)
    doPlot(handles);
return 


function simulateAndCompare3D_Callback(hObject, eventdata, handles)
m = get(handles.modelselection,'Value');
eall = get(handles.experimentselection,'Value');
counter = 1;
plotdatastruct = cell(numel(eall), 1);
for e = eall
    name = handles.plotdata.model(m).experiment(e).name;

        mexmodel = handles.plotdata.model(m).experiment(e).mexmodel;
        timevector = handles.plotdata.model(m).experiment(e).timevector;
        ic = handles.plotdata.model(m).experiment(e).initialconditions;
        % take away the things that are no parameters
        noParamIndices = strmatchIQM('No Parameter',handles.parammodel(m).names,'exact');
        paramIndices = setdiff([1:length(handles.parammodel(m).names)],noParamIndices);
        paramnames = handles.parammodel(m).names(paramIndices);
        paramvalues = handles.parammodel(m).values(paramIndices);

    simdata = IQMPsimulate(mexmodel,timevector,ic,paramnames,paramvalues,handles.integratoroptions);
    % remove '_' from name
    name = strrep(name,'_',' ');
    plotdatastruct{counter} = createdatastruct2IQMplotIQM(simdata,name);
    counter = counter + 1;
end

measuredatastruct = cell(numel(eall), 1);
counter = 1;
for e = eall
    measuredatastruct{counter} = handles.plotdata.model(m).experiment(e).measurement;
    counter = counter + 1;
end

IQMplot3D(plotdatastruct, measuredatastruct);


% Update handles structure
guidata(hObject, handles);

return

function plot_derivatives_Callback(hObject, eventdata, handles)
    m = get(handles.modelselection,'Value');
    eall = get(handles.experimentselection,'Value');
    callnames = get(handles.componentselection,'String');
    callselected = get(handles.componentselection,'Value');
    plotdata = handles.plotdata;
    simTable = cell(numel(plotdata.model(m).experiment(1).timevector), numel(eall)*numel(callselected) + 1);
    simTable(:, 1) = num2cell(plotdata.model(m).experiment(1).timevector);
   
    for sp = 1:numel(eall),
        e = eall(sp);

        edata = plotdata.model(m).experiment(e);
        % general information
        titletext = regexprep(edata.name,'_',' ');
        xlabeltext = 'Time';
        % simulated data information
        sim_timevector = edata.timevector;
        sim_componentnames = edata.componentnames;
        sim_componentvalues = edata.componentvalues;
        for k=1:length(edata.componentnames),
            if isempty(strmatchIQM(edata.componentnames{k},callnames(callselected),'exact')),
                sim_componentvalues(:,k) = NaN;
            end
        end
        
        for kplotsim=1:numel(callselected)
            simTable(:, (sp-1)*numel(callselected)+kplotsim + 1) = num2cell(sim_componentvalues(:,callselected(kplotsim)));
            
        end
    end
    timevectors = cell(0, 1);
    sumTimes = 0;
    startIndices = [];
    for sp = 1 : numel(eall)
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
        for meas=1:length(edata.measurement)
            timevectors{end + 1} = edata.measurement(meas).timevector;
            startIndices(end + 1) = sumTimes + 1;
            sumTimes = sumTimes + numel(edata.measurement(meas).timevector);
        end
    end
    % check if times are equal
    allAreEqual = true;
    for i = 2 : numel(timevectors)
        if ~isequal(timevectors{i}, timevectors{i - 1})
           allAreEqual = false;
           break;
        end
    end
    if allAreEqual % timevecotrs will be shared between all the measurements, equivalent to slow kinetics in Aparecium
        measurementTable = cell(numel(timevectors{1}), numel(timevectors) + 1);
        measurementTable(:, 1) = num2cell(timevectors{1});
        startIndices(:) = 1;
    else % timevectors are different and thus each measurement set will have a separete time vector, equivalent to fast kinetics in Aparecium
        measurementTable = cell(sumTimes, numel(timevectors) + 1);
    end
    
    counter = 1;
    for sp = 1 : numel(eall)
        % measured data information
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
            
        for meas=1:length(edata.measurement), % currently works only
            meas_timevector = edata.measurement(meas).timevector;
            meas_componentnames = edata.measurement(meas).componentnames;
            meas_componentvalues = edata.measurement(meas).componentvalues;
            meas_maxvalues = edata.measurement(meas).maxvalues;
            meas_minvalues = edata.measurement(meas).minvalues;
            for k=1:length(edata.measurement(meas).componentnames),
                if isempty(strmatchIQM(edata.measurement(meas).componentnames{k},callnames(callselected),'exact')),
                    meas_componentvalues(:,k) = NaN;
                    %meas_maxvalues(:,k) = NaN;
                    %meas_minvalues(:,k) = NaN;
                end
            end
            for kplotsim=1:numel(callselected)
                measurementTable(startIndices(counter):startIndices(counter) + numel(meas_timevector) - 1, 1) = num2cell(meas_timevector);
                measurementTable(startIndices(counter):startIndices(counter) + numel(meas_timevector) - 1, counter + 1) = num2cell(meas_componentvalues(:, callselected(kplotsim)));
            end
            
%             for kplotsim=1:numel(callselected)
%                 subplot(NROW,NCOL,kplotsim,'Parent',handles.plotpanel);
%                 feval(handles.dataPlotType,meas_timevector,meas_componentvalues(:,callselected(kplotsim)),'*','linewidth',2,'color',colorvector{mod(sp-1,7)+1}); hold on;
%             end

            counter = counter + 1;
        end
    end

% 
%     if ~isnan(handles.boundsaxis),
%         axis(handles.boundsaxis)
%     end
% 
%     end
% 
%     for k=1:numel(callselected)
%         subplot(NROW,NCOL,k,'Parent',handles.plotpanel);
%         expNames = cell(numel(eall)*2, 1);
%         for sp = 1:numel(eall),
%             e = eall(sp);    
%             edata = plotdata.model(m).experiment(e);
%             expNames{sp*2} = edata.name;
%             expNames{sp*2-1} = [edata.name, ' fit'];
%         end
%         if handles.legendflag == 1,
%             hlhlx = legend(expNames);
%             set(hlhlx,'Interpreter','none');
%         end
%         hlhlx = title(sim_componentnames{callselected(k)});
%         set(hlhlx,'Interpreter','none');
%         hlhlx = xlabel(xlabeltext);
%         set(hlhlx,'Interpreter','none');    
%         hold off;
%     end
g = figure('Name','Simulated data');
createTable(g, {}, simTable, 'Buttons', 'on', 'Visible', 'on')
%uitable('data', simTable)
g = figure('Name','Experimental data');
createTable(g, {}, measurementTable, 'Buttons', 'on', 'Visible', 'on')
figure
hold on
simTable = cell2num(simTable);
measurementTable = cell2num(measurementTable);
plot(simTable(2:end, 1), smooth(diff(smooth(simTable(:, 1), simTable(:, 2), 0.05, 'moving'))./diff(simTable(:, 1)), 'moving'))
plot(measurementTable(2:end, 1), smooth(diff(smooth(measurementTable(:, 1), measurementTable(:, 2), 0.05, 'moving'))./diff(measurementTable(:, 1)), 'moving'))

return



function generateTable_Callback(hObject, eventdata, handles)
    m = get(handles.modelselection,'Value');
    eall = get(handles.experimentselection,'Value');
    callnames = get(handles.componentselection,'String');
    callselected = get(handles.componentselection,'Value');
    plotdata = handles.plotdata;
    simTable = cell(numel(plotdata.model(m).experiment(1).timevector), numel(eall)*numel(callselected) + 1);
    simTable(:, 1) = num2cell(plotdata.model(m).experiment(1).timevector);
   
    for sp = 1:numel(eall),
        e = eall(sp);

        edata = plotdata.model(m).experiment(e);
        % general information
        titletext = regexprep(edata.name,'_',' ');
        xlabeltext = 'Time';
        % simulated data information
        sim_timevector = edata.timevector;
        sim_componentnames = edata.componentnames;
        sim_componentvalues = edata.componentvalues;
        for k=1:length(edata.componentnames),
            if isempty(strmatchIQM(edata.componentnames{k},callnames(callselected),'exact')),
                sim_componentvalues(:,k) = NaN;
            end
        end
        
        for kplotsim=1:numel(callselected)
            simTable(:, (sp-1)*numel(callselected)+kplotsim + 1) = num2cell(sim_componentvalues(:,callselected(kplotsim)));
            
        end
    end
    timevectors = cell(0, 1);
    sumTimes = 0;
    startIndices = [];
    for sp = 1 : numel(eall)
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
        for meas=1:length(edata.measurement)
            timevectors{end + 1} = edata.measurement(meas).timevector;
            startIndices(end + 1) = sumTimes + 1;
            sumTimes = sumTimes + numel(edata.measurement(meas).timevector);
        end
    end
    % check if times are equal
    allAreEqual = true;
    for i = 2 : numel(timevectors)
        if ~isequal(timevectors{i}, timevectors{i - 1})
           allAreEqual = false;
           break;
        end
    end
    if allAreEqual % timevecotrs will be shared between all the measurements, equivalent to slow kinetics in Aparecium
        measurementTable = cell(numel(timevectors{1}), numel(timevectors) + 1);
        measurementTable(:, 1) = num2cell(timevectors{1});
        startIndices(:) = 1;
    else % timevectors are different and thus each measurement set will have a separete time vector, equivalent to fast kinetics in Aparecium
        measurementTable = cell(sumTimes, numel(timevectors) + 1);
    end
    
    counter = 1;
    for sp = 1 : numel(eall)
        % measured data information
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
            
        for meas=1:length(edata.measurement), % currently works only
            meas_timevector = edata.measurement(meas).timevector;
            meas_componentnames = edata.measurement(meas).componentnames;
            meas_componentvalues = edata.measurement(meas).componentvalues;
            meas_maxvalues = edata.measurement(meas).maxvalues;
            meas_minvalues = edata.measurement(meas).minvalues;
            for k=1:length(edata.measurement(meas).componentnames),
                if isempty(strmatchIQM(edata.measurement(meas).componentnames{k},callnames(callselected),'exact')),
                    meas_componentvalues(:,k) = NaN;
                    %meas_maxvalues(:,k) = NaN;
                    %meas_minvalues(:,k) = NaN;
                end
            end
            for kplotsim=1:numel(callselected)
                measurementTable(startIndices(counter):startIndices(counter) + numel(meas_timevector) - 1, 1) = num2cell(meas_timevector);
                measurementTable(startIndices(counter):startIndices(counter) + numel(meas_timevector) - 1, counter + 1) = num2cell(meas_componentvalues(:, callselected(kplotsim)));
            end
            
%             for kplotsim=1:numel(callselected)
%                 subplot(NROW,NCOL,kplotsim,'Parent',handles.plotpanel);
%                 feval(handles.dataPlotType,meas_timevector,meas_componentvalues(:,callselected(kplotsim)),'*','linewidth',2,'color',colorvector{mod(sp-1,7)+1}); hold on;
%             end

            counter = counter + 1;
        end
    end

% 
%     if ~isnan(handles.boundsaxis),
%         axis(handles.boundsaxis)
%     end
% 
%     end
% 
%     for k=1:numel(callselected)
%         subplot(NROW,NCOL,k,'Parent',handles.plotpanel);
%         expNames = cell(numel(eall)*2, 1);
%         for sp = 1:numel(eall),
%             e = eall(sp);    
%             edata = plotdata.model(m).experiment(e);
%             expNames{sp*2} = edata.name;
%             expNames{sp*2-1} = [edata.name, ' fit'];
%         end
%         if handles.legendflag == 1,
%             hlhlx = legend(expNames);
%             set(hlhlx,'Interpreter','none');
%         end
%         hlhlx = title(sim_componentnames{callselected(k)});
%         set(hlhlx,'Interpreter','none');
%         hlhlx = xlabel(xlabeltext);
%         set(hlhlx,'Interpreter','none');    
%         hold off;
%     end
g = figure('Name','Simulated data');
createTable(g, {}, simTable, 'Buttons', 'on', 'Visible', 'on')
%uitable('data', simTable)
g = figure('Name','Experimental data');
createTable(g, {}, measurementTable, 'Buttons', 'on', 'Visible', 'on')
%uitable('data', measurementTable)
return

function toggleBrush_Callback(hObject, eventdata, handles)
    brush(handles.SBplot);
return 

function setWeightValue_Callback(hObject, eventdata, handles)
    % workaround to get brushed data
    % get brushed data
    figure_handles = findall(gcf); Handles = figure_handles( arrayfun(@(H) isprop(H, 'BrushData'), figure_handles) );
    info = [num2cell(Handles),arrayfun(@(H) find(H.BrushData), Handles, 'uniform', 0)];
    
    % get brushed data line names
    dataNames = cell(size(info, 1), 1);
    for i = 1 : size(info, 1)
        dataNames{i} = info{i, 1}.DisplayName;
    end
    
    % get experiment names
    m = get(handles.modelselection,'Value');
    plotdata = handles.plotdata;
    eall = get(handles.experimentselection,'Value');
    datasetNames = cell(numel(eall), 1);
    for sp = 1 : numel(eall)
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
        datasetNames{sp} = edata.name;
    end
    
    intersectingNames = intersect(dataNames, datasetNames);
    
    % match experiment names with brushed data names
    correctIndices = [];
    datasetIndices = [];
    datasetTitles = cell(1, 0);
    datasetTitleIndices = [];
    correspondingExperimentalIndices = [];
    for brushableObject = 1 : size(info, 1)
        index = findStringFromCellArray(get(handles.experimentselection, 'String'), info{brushableObject, 1}.DisplayName);
        if isequal(index, -1)
            
        else
            correctIndices(end + 1) = brushableObject;
            datasetIndices(end + 1) = index;
            objectParent = get(info{brushableObject, 1}, 'Parent');
            parentTitle = get(objectParent.Title);
            datasetTitles{end + 1} = parentTitle.String;
            datasetTitleIndices(end + 1) = findStringFromCellArray(handles.plotdata.model(1).allmeascomponents, parentTitle.String);
        end
    end
    
    weightValue = str2num(get(handles.weight, 'String'));
    
    %now indices must exists for all experiments
    for i = 1 : numel(handles.experimentindices)
        handles.estimation.experiments.measurementindices{i, 1} = 1;
    end

    
    for i = 1 : numel(correctIndices)
        if numel(handles.estimation.experiments.measurementweight) < datasetIndices(i) || isempty(handles.estimation.experiments.measurementweight{datasetIndices(i)})
            edata = plotdata.model(m).experiment(datasetIndices(i));
            nrOfMeasurements = numel(edata.measurement.timevector); 
            handles.estimation.experiments.measurementweight{datasetIndices(i)}{1} = ones(nrOfMeasurements, numel(handles.plotdata.model(1).allmeascomponents));           
        end 
        handles.estimation.experiments.measurementweight{datasetIndices(i)}{1}(info{correctIndices(i), 2}, datasetTitleIndices(i)) = weightValue;
    end
    guidata(hObject, handles);
return

function MergeExperiment_Callback(hObject, eventdata, handles)
% This function allows merging selected experimental points to other
% experiments - should be used for cases where some area of measurements is
% duplicates but conditions change after some event.

    % get the selection
    figure_handles = findall(gcf); Handles = figure_handles( arrayfun(@(H) isprop(H, 'BrushData'), figure_handles) );
    info = [num2cell(Handles),arrayfun(@(H) find(H.BrushData), Handles, 'uniform', 0)];
    
    % get brushed data line names
    dataNames = cell(size(info, 1), 1);
    for i = 1 : size(info, 1)
        dataNames{i} = info{i, 1}.DisplayName;
    end
    
    % get experiment names
    m = get(handles.modelselection,'Value');
    plotdata = handles.plotdata;
    eall = get(handles.experimentselection,'Value');
    datasetNames = cell(numel(eall), 1);
    for sp = 1 : numel(eall)
        e = eall(sp);
        edata = plotdata.model(m).experiment(e);
        datasetNames{sp} = edata.name;
    end
    
    intersectingNames = intersect(dataNames, datasetNames);
    
    % match experiment names with brushed data names
    correctIndices = [];
    datasetIndices = [];
    datasetTitles = cell(1, 0);
    datasetTitleIndices = [];
    correspondingExperimentalIndices = [];
    for brushableObject = 1 : size(info, 1)
        index = findStringFromCellArray(get(handles.experimentselection, 'String'), info{brushableObject, 1}.DisplayName);
        if isequal(index, -1)
            
        else
            correctIndices(end + 1) = brushableObject;
            datasetIndices(end + 1) = index;
            objectParent = get(info{brushableObject, 1}, 'Parent');
            parentTitle = get(objectParent.Title);
            datasetTitles{end + 1} = parentTitle.String;
            datasetTitleIndices(end + 1) = findStringFromCellArray(handles.plotdata.model(1).allmeascomponents, parentTitle.String);
        end
    end
    
    % collect the largest overlap between dimensions
    indices = info{correctIndices(1), 2};
    for i = 2 : numel(correctIndices)       
        indices = intersect(indices, info{correctIndices(i), 2});
    end
    
    % calculate mean x values
    
    timeMatrix = zeros(numel(indices), numel(correctIndices));
    
    for i = 1 : numel(correctIndices)
        timeMatrix(:, i) = handles.plotdata.model.experiment(datasetIndices(i)).measurement.timevector(indices);
    end
    
    timeVec = mean(timeMatrix, 2);
    
    for i = 1 : numel(correctIndices)
        handles.plotdata.model.experiment(datasetIndices(i)).measurement.timevector(indices) = timeVec;
    end
    
    valueMatrix = zeros(numel(indices), numel(correctIndices));
    
    for i = 1 : numel(correctIndices)
        valueMatrix(:, i) = handles.plotdata.model.experiment(datasetIndices(i)).measurement.componentvalues(indices, datasetTitleIndices(i));
    end
    
    valueVec = mean(valueMatrix, 2);
    
    for i = 1 : numel(correctIndices)
        handles.plotdata.model.experiment(datasetIndices(i)).measurement.componentvalues(indices, datasetTitleIndices(i)) = valueVec;
    end

    guidata(hObject, handles)
return

function RegressionEstimation_Callback(hObject, eventdata, handles)
    mousePointCoords = ginput(2);
    % plot the mouse point coordinates on the figure
    hold on;
    plot(mousePointCoords(:,1), mousePointCoords(:,2),'--o','Color', [0, 0, 0],'MarkerSize',8);
    c = [[1; 1]  mousePointCoords(:,1)]\mousePointCoords(:,2);                        % Calculate Parameter Vector
    slope_m = c(2);
    intercept_b = c(1);
    set(handles.uitable1, 'data', [slope_m, intercept_b]);
return