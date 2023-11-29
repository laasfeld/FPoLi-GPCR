function [modeling_result_list, experiment_list] = ModelSingleCompetitionExperiment(receptorVolume, FL_concentration, competitive_conc_start, dilution_ratio, kf1, Kd1, kf2, Kd2, Rstock, proj, mod, timeStart, timeStep, timeEnd, addNoise)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Receptor_Volumes = receptorVolume; % in µL
FL_concentrations = FL_concentration;
Competitive_ligand_concentrations = 10.^(competitive_conc_start : log10(dilution_ratio) : competitive_conc_start + 7*log10(dilution_ratio)); %nM

integratoroptions = struct();
integratoroptions.abstol = 1.0000e-06;
integratoroptions.reltol = 1.0000e-06;
integratoroptions.minstep = 0;
integratoroptions.maxstep = Inf;
integratoroptions.maxnumsteps =  1000;


editable_experiment = IQMexperiment();
                
% Set up the model with specific biochemical parameters
mod.parameters(index_of_name(mod, 'kf1')).value = kf1;
mod.parameters(index_of_name(mod, 'Kd1')).value = Kd1;
mod.parameters(index_of_name(mod, 'Arl')).value = 0.15;%0.3;
mod.parameters(index_of_name(mod, 'Afree')).value = 0.04;
mod.parameters(index_of_name(mod, 'Anbv')).value = 0.15;
mod.parameters(index_of_name(mod, 'n')).value = Rstock/100;
mod.parameters(index_of_name(mod, 'ultracorrect')).value = 1;
mod.parameters(index_of_name(mod, 'QYRL')).value = 1;
mod.parameters(index_of_name(mod, 'QYNBVL')).value = 1;
mod.parameters(index_of_name(mod, 'kf3')).value = 0;
mod.parameters(index_of_name(mod, 'kr3')).value = 0;

mod.parameters(index_of_name(mod, 'kf2')).value = kf2;
mod.parameters(index_of_name(mod, 'Ki2')).value = Kd2;

mod.parameters(index_of_name(mod, 'NBVcommon')).value = 10000;

mod.parameters(index_of_name(mod, 'kf4')).value = 0.16;
mod.parameters(index_of_name(mod, 'kr4')).value = 300000;

model_obj = IQMmodel(mod);

experiment_list = cell(16, 1);
well_order_list = cell(16, 1);
modelexp_list = cell(16, 1);
modeling_result_list = cell(16, 1);
experiment_counter = 1;
global mexmodel
for Competitive_ligand_concentration_index = 1 : numel(Competitive_ligand_concentrations)
    for duplicate_index = 1 : 2
        well_order = (Competitive_ligand_concentration_index - 1)*2 + duplicate_index;

        disp(['Simulating - ', num2str(Competitive_ligand_concentration_index), ' ', num2str(duplicate_index), ' ', num2str(experiment_counter)]);
        experiment_struct = struct(editable_experiment);
        experiment_struct.paramicsettings(1).name = 'V';
        experiment_struct.paramicsettings(2).name = 'L';
        experiment_struct.paramicsettings(3).name = 'C';

        experiment_struct.paramicsettings(1).icflag = 1;
        experiment_struct.paramicsettings(2).icflag = 1;
        experiment_struct.paramicsettings(3).icflag = 1;
        if addNoise
            experiment_struct.paramicsettings(1).formula = num2str(Receptor_Volumes(1)*(1+0.02*randn()));
            experiment_struct.paramicsettings(2).formula = num2str(FL_concentrations(1)*(1+0.02*randn()));
            experiment_struct.paramicsettings(3).formula = num2str(Competitive_ligand_concentrations(Competitive_ligand_concentration_index)*Kd2*(1+0.02*randn()));
        else
            experiment_struct.paramicsettings(1).formula = num2str(Receptor_Volumes(1));
            experiment_struct.paramicsettings(2).formula = num2str(FL_concentrations(1));
            experiment_struct.paramicsettings(3).formula = num2str(Competitive_ligand_concentrations(Competitive_ligand_concentration_index)*Kd2);
        end

        experiment_list{experiment_counter} = experiment_struct;
        modelexp_list{experiment_counter} = IQMmergemodexp(model_obj, IQMexperiment(experiment_struct));
        if ~exist('mexmodel', 'var') || isempty(mexmodel)
            [mexmodel, mexmodelFullpath] = IQMmakeTempMEXmodel(modelexp_list{experiment_counter});
        end
        global mexmodel

        %struct(simulations.modelexp_list{well_order}).states(2).initialCondition
        well_order_list{experiment_counter} = well_order;

        %% old version
        %modeling_result_list{experiment_counter} = IQMsimulate(modelexp_list{experiment_counter}, 0:5:54000).variablevalues(:,5);
        %%
        %% mex version
        icCell = cell(numel(IQMstruct(modelexp_list{experiment_counter}).states), 1);
        [icCell{:}] = IQMstruct(modelexp_list{experiment_counter}).states(:).initialCondition;
        ic = cell2mat(icCell);

        pvCell = cell(numel(IQMstruct(modelexp_list{experiment_counter}).parameters), 1);
        [pvCell{:}] = IQMstruct(modelexp_list{experiment_counter}).parameters(:).value;
        pv = cell2mat(pvCell);

        simdata = feval(mexmodel,[timeStart: timeStep: timeEnd],ic,pv,integratoroptions);
        if addNoise
            modeling_result_list{experiment_counter} = simdata.variablevalues(:,5)+0.0015*randn(size(simdata.variablevalues(:,5)));
        else
            modeling_result_list{experiment_counter} = simdata.variablevalues(:,5);           
        end
        
        experiment_counter = experiment_counter + 1;
    end
end

end

