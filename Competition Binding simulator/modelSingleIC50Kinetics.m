function timeIC50CellArray = modelSingleIC50Kinetics(modeling_result_list, experiment_list, timeStart, timeStep, timeEnd, samplingTimes, plateReadingMode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    competitorConcentrations = zeros(16,1);
    for competitorWellIndex = 1 : 16
        competitorConcentrations(competitorWellIndex) = str2num(experiment_list{competitorWellIndex}.paramicsettings(3).formula);
    end
    
    timeIC50CellArray = cell(numel(timeStart : timeStep :  timeEnd), 1);
    timeIC50CellArrayCounter = 1;
    for cycleStart = timeStart : timeStep :  timeEnd
        timeSamplingEstimation = zeros(16, 1);
        timeSamplingCounter = 1;
        correspondingTimeSamplingindex = zeros(16, 1);
        sampleTimes = zeros(16, 1);
        FAvalues = zeros(16, 1);
        for concentrationIndex = 0 : 7
            for replicate = 0 : 1
                if strcmp(plateReadingMode, 'Instantaneous')
                    timeSamplingEstimation(timeSamplingCounter) = cycleStart;%+ ((concentrationIndex*12+replicate)*1.5)%+144-((concentrationIndex*12+replicate)*1.5);  
                elseif strcmp(plateReadingMode, 'Top down dilution')
                    timeSamplingEstimation(timeSamplingCounter) = cycleStart + ((concentrationIndex+replicate*8)*1.5);%+ ((concentrationIndex*12+replicate)*1.5)%+144-((concentrationIndex*12+replicate)*1.5);  
                elseif strcmp(plateReadingMode, 'Inversed dilution')
                    timeSamplingEstimation(timeSamplingCounter) = cycleStart + 24-((concentrationIndex+replicate*8)*1.5);%+ ((concentrationIndex*12+replicate)*1.5)%+144-((concentrationIndex*12+replicate)*1.5);  
                end
                %

                [value, index] = min(abs(samplingTimes - timeSamplingEstimation(timeSamplingCounter)));
                correspondingTimeSamplingindex(timeSamplingCounter) = index;
                FAvalues(timeSamplingCounter) = modeling_result_list{timeSamplingCounter}(index);
                sampleTimes(timeSamplingCounter) = samplingTimes(index);
                timeSamplingCounter = timeSamplingCounter + 1;
            end
        end
           
        %plot(log10(competitorConcentrations+0.000000001), FAvalues, 'o');
        
        %ylim([0 0.3])
        %drawnow();
        %pause(0.1);
        cycleStart
        timeIC50CellArray{timeIC50CellArrayCounter} = [sampleTimes, log10(competitorConcentrations+0.000000001), FAvalues, ];
        timeIC50CellArrayCounter = timeIC50CellArrayCounter + 1;
    end
end

