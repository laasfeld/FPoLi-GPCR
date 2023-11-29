function [outputArg1,outputArg2] = runSingleSimulation(kf_index, Kd_index, run, first_pipetting_to_meas_start, pipeting_time_deltas, Cycle_time, measurement_time_delay_deltas, edited_model, addPipetingNoise, addFANoise)
save_file_name = fullfile('C:\Users\laasfeld\Aparecium paper simulations\fits_no_corrections_added_noise', ['run_', num2str(run),'_', num2str(kf_index), '_', num2str(Kd_index), '.mat']);
if exist(save_file_name, 'file')
   return 
end
proj = IQMprojectSB('C:\Users\laasfeld\Aparecium paper simulations\project_shell.iqmp');

proj_struct = struct(proj);
proj_struct.experiments = [];
proj_struct.models = {edited_model};
proj_struct.estimations = [];
simulations = load(fullfile('C:\Users\laasfeld\Aparecium paper simulations\simulations', [num2str(kf_index), '_', num2str(Kd_index), '.mat']));
experiments = load(fullfile('C:\Users\laasfeld\Aparecium paper simulations\simulations', ['exp_', num2str(kf_index), '_', num2str(Kd_index), '.mat']));

for well_index = 1 : 64

    proj_struct.experiments(well_index).name = ['well ', num2str(well_index)];
    proj_struct.experiments(well_index).notes = '';
    
    if addPipetingNoise
        for treatmentIndex = 1 : numel(experiments.experiment_list{well_index}.paramicsettings)
            ic_value = str2double(experiments.experiment_list{well_index}.paramicsettings(treatmentIndex).formula);
            experiments.experiment_list{well_index}.paramicsettings(treatmentIndex).formula = num2str(ic_value + ic_value * randn() * 0.03);
        end
    end
    
    proj_struct.experiments(well_index).experiment = IQMexperiment(experiments.experiment_list{well_index});

    measurementStruct.name = num2str(well_index);
    measurementStruct.notes = 'Automatically generated note';

    well_total_timeshift = (first_pipetting_to_meas_start - pipeting_time_deltas(well_index)) + measurement_time_delay_deltas(well_index);
    real_measurement_times = well_total_timeshift : Cycle_time : 10000;           
    time = real_measurement_times - well_total_timeshift; %- well_total_timeshift + %measurement_time_delay_deltas(well_index) + (first_pipetting_to_meas_start - max(pipeting_time_deltas));
    measurements = simulations.modeling_result_list{well_index}(round(real_measurement_times) + 1); % simulation started from 0 seconds, index starts from 1

    measurement_data_struct.name = 'FA';
    measurement_data_struct.notes = 'automatically generated notes';
    
    if addFANoise
        measurements = measurements + measurements.*randn(size(measurements))*0.01*sqrt( 3/str2double(experiments.experiment_list{well_index}.paramicsettings(2).formula) );
    end
    
    measurement_data_struct.values = measurements;
    measurement_data_struct.maxvalues = measurements;
    measurement_data_struct.minvalues = measurements;


    measurementStruct.time = time';
    measurementStruct.data = measurement_data_struct;

    proj_struct.experiments(well_index).measurements = {IQMmeasurement(measurementStruct)};
end

estimation.modelindex = 1;
estimation.experiments.indices = 1:64;
estimation.experiments.measurementindices = cell(0,0);
estimation.experiments.weight = ones(1, 64);
estimation.experiments.measurementweight = cell(0,0);        
estimation_parameters.names = {'kf1'; 'Kd1'; 'kf2'; 'Ki2'; 'Afree'; 'Arl'; 'Anbv'; 'n'; 'NBVcommon'};

% set all the lower and upper bounds to be 10 fold higher or lower
% compared to the original value
lowbounds = zeros(0, 0);
highbounds = zeros(0, 0);
modelexp = struct(simulations.modelexp_list{1});
for parameter_index = 1 : numel(estimation_parameters.names)
    name_index = index_of_name(modelexp, estimation_parameters.names{parameter_index});
    if strcmp(estimation_parameters.names{parameter_index}, 'kf1')
        lowbounds(parameter_index) = 10^-5;
        highbounds(parameter_index) = 100;
    elseif strcmp(estimation_parameters.names{parameter_index}, 'Anbv') || strcmp(estimation_parameters.names{parameter_index}, 'Arl') || strcmp(estimation_parameters.names{parameter_index}, 'Afree')
        lowbounds(parameter_index) = 0;
        highbounds(parameter_index) = 0.4;
    else
        lowbounds(parameter_index) = modelexp.parameters(name_index).value / 10;
        highbounds(parameter_index) = modelexp.parameters(name_index).value * 10;
    end
end
estimation_parameters.lowbounds = lowbounds;
estimation_parameters.highbounds = highbounds;

estimation.parameters = estimation_parameters;

estimation.parameterslocal.names = cell(0,0);
estimation.parameterslocal.lowbounds = [];
estimation.parameterslocal.highbounds = [];

estimation.initialconditions.names = cell(0,0);
estimation.initialconditions.lowbounds = [];
estimation.initialconditions.highbounds = [];

estimation.optimization.method = 'simannealingIQM';

optimization_options.tempstart = 100;
optimization_options.tempend = 0.1;
optimization_options.tempfactor = 0.1;
optimization_options.maxitertemp = 100;
optimization_options.maxitertemp0 = 500;
optimization_options.maxtime = 500;

estimation.optimization.options = optimization_options;
estimation.integratior.options.abstol = 10^-6;
estimation.integratior.options.reltol = 10^-6;
estimation.integratior.options.minstep = 0;
estimation.integratior.options.maxstep = Inf;
estimation.integratior.options.maxnumsteps = 1000;

estimation.initialconditionsFlag = 1;
estimation.displayFlag = 2;
estimation.scalingFlag = 1;
estimation.timescalingFlag = 2;

proj_struct.estimations = {estimation};
output = IQMparameterestimation(IQMprojectSB(proj_struct),estimation,1);
save(save_file_name, 'output', 'estimation');

''
end

