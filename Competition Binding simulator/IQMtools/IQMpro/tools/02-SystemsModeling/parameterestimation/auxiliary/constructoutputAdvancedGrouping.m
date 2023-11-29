%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the output argument of the parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = constructoutputAdvancedGrouping(projectstruct,estimation,parametersAdvancedGrouping,initialconditions,experiments,workestimation,Popt,ICopt,FVALopt)
% construct optimized project 
% 1) add optimized parameters which are shared by all experiments to the
% model
% 2) add other optimized parameters and initial conditions to the experiment descriptions
projectoptstruct = projectstruct;

% find how many different models must be generated from this fit
% convert to table to access table functionality
t = cell2table(estimation.parametersadvanced.parameterGroupingTable);
groups = findgroups(t);
uniqueGroups = unique(groups);

differentiatingTableIndices = zeros(0, 1);
for columnIndex = 1 : size(t, 2)
    if numel(unique(t{:, columnIndex})) > 1
        differentiatingTableIndices = [differentiatingTableIndices, columnIndex];
    end
end

for uniqueGroupIndex = 1 : numel(uniqueGroups)
    % find correct parameter names
    
    groupExample = t{find(groups == uniqueGroups(uniqueGroupIndex), 1), :};
    LIA = ismember(parametersAdvancedGrouping.paramGroup,groupExample);
    projectoptstruct.models{end + 1} = IQMparameters(projectoptstruct.models{estimation.modelindex},parametersAdvancedGrouping.names(LIA),Popt(LIA));
    differentiatingString = '';
    for differentiatingIndex = 1 : numel(differentiatingTableIndices)
        groupNames = unique(t{:, differentiatingTableIndices(differentiatingIndex)});
        memberIndex = ismember(unique(t{:, differentiatingTableIndices(differentiatingIndex)}), parametersAdvancedGrouping.paramGroup(LIA));
        differentiatingString = [differentiatingString, '_', groupNames{memberIndex}];
    end
    
    m = struct(projectoptstruct.models{end});
    m.name = [m.name, differentiatingString];
    projectoptstruct.models{end} = IQMmodel(m);
    
    icnames = initialconditions.names;
    for k=1:length(experiments),
        experiment = experiments(k).experiment;
        experimentstruct = struct(experiment);
        allpresenticnames = {experimentstruct.paramicsettings.name}; % gets icnames and parameternames mixed! but thats ok ...
        % optimized initialconditions
        for k2=1:length(icnames),
            % check if initial condition already defined in the experiment
            index = strmatchIQM(icnames{k2},allpresenticnames,'exact');
            if ~isempty(index),
                % overwrite
                experimentstruct.paramicsettings(index).name = icnames{k2};
                experimentstruct.paramicsettings(index).formula = sprintf('%g',ICopt(workestimation(k).indicesinICvector(k2)));
                experimentstruct.paramicsettings(index).notes = 'optimized value';
                experimentstruct.paramicsettings(index).icflag = 1;
            else
                % append
                experimentstruct.paramicsettings(end+1).name = icnames{k2};
                experimentstruct.paramicsettings(end).formula = sprintf('%g',ICopt(workestimation(k).indicesinICvector(k2)));
                experimentstruct.paramicsettings(end).notes = 'optimized value';
                experimentstruct.paramicsettings(end).icflag = 1;
            end
        end
   
        % get index of experiment in the project
        indexexp = estimation.experiments.indices(k);
        projectoptstruct.experiments(indexexp).experiment = IQMexperiment(experimentstruct);
    end
end
projectopt = IQMprojectSB(projectoptstruct);
% output is a structure
output = [];
output.parameters = parametersAdvancedGrouping.names;
output.Popt = Popt;
output.parameterslocal = cell(0, 1);
output.PLOCALopt = [];
output.icnames = initialconditions.names;
if ~isempty(initialconditions.names),
    output.ICopt = reshape(ICopt',length(initialconditions.names),length(ICopt)/length(initialconditions.names))';
else 
    output.ICopt = [];
end
output.FVALopt = FVALopt;
output.projectopt = projectopt;
output.estimation = estimation;
return

