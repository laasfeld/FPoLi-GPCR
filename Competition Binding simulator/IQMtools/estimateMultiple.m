%%
filenames = cell(1000, 1);
estimations = cell(0,0);
% handles.UserData.project, project
% getEstimationStructure(handles), estimations{end+1}
%%
for estimationIndex = 1 : numel(estimations)
    output = IQMparameterestimation(project,estimations{estimationIndex},1);
    ps = struct(output.projectopt);
    model = ps.models{1};
    filename = [path, filenames{estimationIndex}];
    IQMcreateTEXTfile(model,filename);
end
    