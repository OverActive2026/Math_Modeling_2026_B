function test_q4_online_precooling_cases23_results()
%TEST_Q4_ONLINE_PRECOOLING_CASES23_RESULTS Validate saved case 2/3 solves.
% Run run_q4_online_precooling_cases23 before calling this test.
    codeDir = fileparts(mfilename('fullpath'));
    for caseNumber = [2 3]
        resultFile = fullfile(codeDir,'results_control',sprintf( ...
            'q4_online_case%d_tuning_melt.mat',caseNumber));
        if ~isfile(resultFile)
            error('q4_cases23_test:MissingResult', ...
                'Run run_q4_online_precooling_cases23 first.');
        end
        saved = load(resultFile,'online','reference','initialState', ...
            'durationMin','simOpt');
        assert(saved.durationMin == 20*(caseNumber-1));
        assert(saved.simOpt.kMelt == 1);
        for r = {saved.online,saved.reference}
            trajectory = r{1};
            assert(trajectory.success);
            assert(all(trajectory.finalTemperatureC > 0));
            assert(max(abs(trajectory.temperatureC(1,:)- ...
                saved.initialState.cellAverageTemperatureC(:).')) < 1e-9);
            assert(trajectory.minimumCellVoltageV >= ...
                trajectory.options.minimumVoltageV-1e-7);
            assert(trajectory.maximumIceVolumeFraction < ...
                trajectory.options.maximumIceVolumeFraction);
            assert(all(trajectory.powerDensityWcm2(:) >= 0 & ...
                trajectory.powerDensityWcm2(:) <= 1));
            durationS = trajectory.intervalEndS- ...
                trajectory.intervalStartS;
            reconstructedJ = trajectory.options.cellAreaM2*1e4* ...
                sum(sum(trajectory.powerDensityWcm2.*durationS));
            assert(abs(reconstructedJ- ...
                trajectory.totalAuxiliaryEnergyJ) < 1e-7);
        end
        if caseNumber == 3
            % Warm inner cells must lose previously accumulated ice.
            icePeak = max(saved.online.iceVolumeFraction,[],1);
            iceFinal = saved.online.iceVolumeFraction(end,:);
            assert(all(icePeak(2:4)-iceFinal(2:4) > 0.01));
        end
    end
    fprintf('test_q4_online_precooling_cases23_results PASS\n');
end
