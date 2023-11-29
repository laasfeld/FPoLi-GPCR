%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD PARAMETER INDICES AND NOMINAL VALUES TO WORKESTIMATION STRUCTURE
% Could be different indices for different experiments ... so its done here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [workestimation, parametersAdvancedGrouping] = addparameterindices(workestimation,parameters)
[allparameters,nominalparamvalues] = IQMparameters(workestimation(1).IQMmodel);
parametersNominal = [];
parametersNames = [];
parametersAdvancedLowerBounds = [];
parametersAdvancedHigherBounds = [];
paramGroup = [];
parameterIndexRanges = cell(size(parameters.parameterGroupingTable, 2), 1);
for parameter = 1 : size(parameters.parameterGroupingTable, 2)
    parameterIndexRanges{parameter}(1) = numel(parametersNominal)+1;
    %parametersNominal = [parametersNominal, ones(1, numel(unique(parameters.parameterGroupingTable(:, parameter)))) * nominalparamvalues(parameter)];
    additionalParamNames = cell(1, numel(unique(parameters.parameterGroupingTable(:, parameter))));
    for index = 1 : numel(unique(parameters.parameterGroupingTable(:, parameter)))
        additionalParamNames{index} = allparameters{parameter};
        parametersNominal = [parametersNominal, parameters.nominalValue{parameter}(index)];
    end
    parametersNames = [parametersNames, additionalParamNames];
    parameterIndexRanges{parameter}(2) = numel(parametersNominal);
    
    uniqueParamGroups = unique(parameters.parameterGroupingTable(:, parameter));
    possibleParameters = parameters.paramGroupNames{parameter};
    for group = 1 : numel(uniqueParamGroups)
        groupIndex = find(strcmp(possibleParameters, uniqueParamGroups{group}));
        parametersAdvancedLowerBounds = [parametersAdvancedLowerBounds, parameters.lowbounds{parameter}(groupIndex) ];
        parametersAdvancedHigherBounds = [parametersAdvancedHigherBounds,  parameters.highbounds{parameter}(groupIndex) ];
        numel(parametersAdvancedHigherBounds)
        paramGroup = [paramGroup, uniqueParamGroups(group)];
    end
    
end

for k=1:length(workestimation),
    [allparameters,nominalparamvalues] = IQMparameters(workestimation(k).IQMmodel);
    paramindices = getnamevectorindices(allparameters,parameters.names);
    workestimation(k).paramindices = [];
    workestimation(k).paramindiceslocal = [];
    workestimation(k).paramindicesAdvanced = paramindices; 
    workestimation(k).paramnominal = nominalparamvalues;
end

for k=1:length(workestimation)
    indicesInAdvancedParamVector = [];
    for parameter = 1 : size(parameters.parameterGroupingTable, 2)
        uniqueParams = unique(parameters.parameterGroupingTable(:, parameter));
        thisParam = parameters.parameterGroupingTable(k, parameter);
        resultIndex = find(strcmp(uniqueParams, thisParam));
        indicesInAdvancedParamVector = [indicesInAdvancedParamVector, parameterIndexRanges{parameter}(1)+resultIndex - 1];        
    end
    workestimation(k).indicesInAdvancedParamVector = indicesInAdvancedParamVector;
end
parametersAdvancedGrouping.names = parametersNames;
parametersAdvancedGrouping.pliv = parametersNominal;
parametersAdvancedGrouping.pllowerbounds = parametersAdvancedLowerBounds;
parametersAdvancedGrouping.plhigherbounds = parametersAdvancedHigherBounds;
parametersAdvancedGrouping.paramGroup = paramGroup;

return
