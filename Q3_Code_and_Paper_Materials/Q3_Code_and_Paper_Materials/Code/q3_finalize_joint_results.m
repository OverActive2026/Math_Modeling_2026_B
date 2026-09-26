function summary=q3_finalize_joint_results(reuseUnchanged)
% Verify the best physical records and write the question's comparison table.
    if nargin<1, reuseUnchanged=false; end
    codeDir=fileparts(mfilename('fullpath'));
    outDir=fullfile(codeDir,'results');
    q3_export_joint_history();
    a=load(fullfile(outDir,'q3_min_energy_preheat.mat'),'run');
    pre=a.run.best;
    d=load(fullfile(outDir,'q3_joint_de_checkpoint.mat'),'run');
    assert(strcmp(d.run.best.status,'success'),'No feasible DE solution.');
    co=d.run.best;
    targetedPath=fullfile(outDir,'q3_targeted_checkpoint.mat');
    targetedEvaluations=0;
    if exist(targetedPath,'file')
        targeted=load(targetedPath,'run');
        targetedEvaluations=targeted.run.evaluations;
        if targeted.run.bestPre.energyJ<pre.energyJ
            p=targeted.run.bestPre;
            pre=struct('q',p.q,'thRequestedS',p.th,'energyJ',p.energyJ);
        end
        if targeted.run.bestCo.energyJ<co.energyJ
            co=targeted.run.bestCo;
        end
    end
    refinements=[dir(fullfile(outDir,'q3_small_allocation_*.mat')); ...
        dir(fullfile(outDir,'q3_end_margin_*.mat')); ...
        dir(fullfile(outDir,'q3_preheat_timing_*.mat')); ...
        dir(fullfile(outDir,'q3_budget_case_*.mat')); ...
        dir(fullfile(outDir,'q3_multistart_best.mat')); ...
        dir(fullfile(outDir,'q3_multistart_stop_*.mat'))];
    for k=1:numel(refinements)
        point=load(fullfile(outDir,refinements(k).name),'c');
        if strcmp(point.c.mode,'preheat')
            if ~strcmp(point.c.status,'success') && ~isempty(point.c.result)
                fprintf('Preheat timing diagnostic: th=%.9f Tmin=%.9f reason=%s\n', ...
                    point.c.th,min(point.c.result.TavgC(end,:)),point.c.reason);
            end
            if strcmp(point.c.status,'success') && point.c.energyJ<pre.energyJ
                pre=struct('q',point.c.q,'thRequestedS',point.c.th, ...
                    'energyJ',point.c.energyJ);
            end
        elseif strcmp(point.c.status,'success') && point.c.energyJ<co.energyJ
            co=point.c;
        end
    end
    cfg=q3_optimization_config('smoke'); opt=cfg.simOpt;
    opt.RelTol=1e-6; opt.AbsTol=1e-10; opt.MaxStep=.0125;
    verified=cell(3,1);
    modes={'preheat','coheat','coheat'};
    inputQ={pre.q,co.q,ones(1,5)};
    % Remove the redundant requested-heating tail beyond first success.
    % A small time margin keeps the planned shutoff just after the event;
    % the shortened common duration is independently simulated below.
    coRequestedTh=min(co.th,co.startupTimeS+1e-5);
    inputTh=[pre.thRequestedS,coRequestedTh,40];
    for i=1:3
        fprintf('Strict verification %d/3: %s\n',i,modes{i});
        checkpointPath=fullfile(outDir,sprintf('q3_strict_verification_%02d.mat',i));
        reuse=false;
        if reuseUnchanged && exist(checkpointPath,'file')
            old=load(checkpointPath,'checkpoint','opt');
            snapshot=dir(checkpointPath);
            physics={'pemfc_stack5_simulate_q3.m','pemfc_setup_ice.m', ...
                'thermal_temperature_state_ice.m','water_ice_state_ice.m', ...
                'gas_transport_state_ice.m','q3_auxiliary_heating.m', ...
                'q2_current_strategy.m','q3_evaluate_energy6.m'};
            fresh=true;
            for f=1:numel(physics)
                source=dir(fullfile(codeDir,physics{f}));
                fresh=fresh && source.datenum<=snapshot.datenum;
            end
            reuse=fresh && isequaln(old.opt,opt) && ...
                strcmp(old.checkpoint.mode,modes{i}) && ...
                strcmp(old.checkpoint.status,'success') && ...
                isequal(old.checkpoint.q,inputQ{i}) && old.checkpoint.th==inputTh(i);
            if reuse, verified{i}=old.checkpoint; end
        end
        if ~reuse
            verified{i}=q3_evaluate_energy6(modes{i},inputQ{i},inputTh(i),opt);
        else
            fprintf('Reused unchanged completed strict verification.\n');
        end
        checkpoint=verified{i}; %#ok<NASGU>
        save(fullfile(outDir,sprintf('q3_strict_verification_%02d.mat',i)), ...
            'checkpoint','opt','-v7.3');
    end
    for i=1:3
        assert(strcmp(verified{i}.status,'success'), ...
            'Independent recheck failed: %s',verified{i}.reason);
    end
    allocationRecheck=verified{2};
    if verified{3}.energyJ < verified{2}.energyJ || ...
            (verified{3}.energyJ==verified{2}.energyJ && ...
            verified{3}.startupTimeS<verified{2}.startupTimeS)
        verified{2}=verified{3};
    end
    summary=table;
    summary.mode=["preheat";"coheat"];
    Q=vertcat(verified{1}.q,verified{2}.q);
    for k=1:5, summary.(sprintf('q%d',k))=Q(:,k); end
    summary.thRequestedS=[verified{1}.th;verified{2}.th];
    R={verified{1}.result,verified{2}.result};
    summary.thUsedS=[R{1}.actualHeatingTimeS;R{2}.actualHeatingTimeS];
    energies=vertcat(R{1}.auxiliaryEnergyJCell,R{2}.auxiliaryEnergyJCell);
    for k=1:5, summary.(sprintf('E%dJ',k))=energies(:,k); end
    summary.ETotalJ=sum(energies,2);
    summary.startupTimeS=[R{1}.successTimeS;R{2}.successTimeS];
    summary.minimumVoltageV=[R{1}.minimumCellVoltageV;R{2}.minimumCellVoltageV];
    summary.maximumIceVolumeFraction= ...
        [R{1}.maximumIceVolumeFraction;R{2}.maximumIceVolumeFraction];
    summary.endCellMaxIce=[max(R{1}.maxIceVolumeFractionCell(:,[1 5]),[],'all'); ...
        max(R{2}.maxIceVolumeFractionCell(:,[1 5]),[],'all')];
    summary.centerCellMaxIce=[max(R{1}.maxIceVolumeFractionCell(:,3)); ...
        max(R{2}.maxIceVolumeFractionCell(:,3))];
    summary.chargeCcm2=[R{1}.chargeUsedCcm2;R{2}.chargeUsedCcm2];
    summary.success=true(2,1);
    summary.resultStatus=repmat("best_verified_finite_search",2,1);
    writetable(summary,fullfile(outDir,'q3_final_energy_table.csv'),'Encoding','UTF-8');
    save(fullfile(outDir,'q3_final_energy_table.mat'),'summary','verified','cfg', ...
        'allocationRecheck','-v7.3');
    fid=fopen(fullfile(outDir,'q3_energy_results_report.md'),'w','n','UTF-8');
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid,'# 问题3辅助能耗优化结果\n\n');
    fprintf(fid,'以下是当前搜索找到并经完整物理模型复算的最低能耗可行候选，尚无全局最优性证明。\n\n');
    fprintf(fid,'初温与环境均为−30 ℃，gammaIce=3.5，kFreeze=0.4，采用原粗网格。五片功率各不超过1 W/cm²。\n\n');
    fprintf(fid,'|方式|五片功率 W/cm²|实际加热 s|启动 s|各片能耗 J|总能耗 J|端片最大冰体积分数|\n');
    fprintf(fid,'|---|---|---:|---:|---|---:|---:|\n');
    names={'纯预热','恒功率协同'};
    for i=1:2
        fprintf(fid,'|%s|%s|%.6f|%.6f|%s|%.3f|%.6f|\n',names{i}, ...
            mat2str(Q(i,:),5),summary.thUsedS(i),summary.startupTimeS(i), ...
            mat2str(energies(i,:),7),summary.ETotalJ(i),summary.endCellMaxIce(i));
    end
    fprintf(fid,'\n纯预热以五片均超过0 ℃作为启动时刻。协同从零时刻按题定斜坡加载；停热时刻独立搜索，停热后仍观察到启动、低压等物理失败或20 C/cm²电荷上限。\n');
    fprintf(fid,'\n协同搜索采用六变量差分进化，已累计评估%d次，当前变异代号为%d，当前代待处理序号为%d。日志与状态保存在q3_joint_de_checkpoint.mat。\n', ...
        d.run.evaluations,d.run.generation,d.run.pendingIndex);
    fprintf(fid,'\n最终记录严格按最低辅助能耗选择，能耗相等时比较启动时间；差分进化的0.1 J容差只用于种群替换，不允许把全局最低能耗记录逐次抬高。未另设题外启动期限或要求各片能耗必须相等/不等。\n');
    fprintf(fid,'\n另完成%d次定向物理评价，检查此前数值未判定候选、多片联合降功率、共同提前停热、稀疏长加热和均匀低功率长加热，并围绕新可行点逐片搜索。明细见q3_targeted_history.csv。有限次数不构成收敛或全局最优证明。\n',targetedEvaluations);
    budgetPath=fullfile(outDir,'q3_budget_checkpoint.mat');
    if exist(budgetPath,'file')
        budgetRun=load(budgetPath,'run');
        fprintf(fid,'\n能量预算补充搜索另完成%d项候选，覆盖0.4、0.5、0.6、0.7、0.9均匀功率及两端较低功率分配，并缩小端片功率步长。详见q3_budget_history.csv。较优局部方案可在代际边界进入差分进化种群；相关导入记录保存在eliteImports字段。\n',budgetRun.run.evaluations);
    end
    fprintf(fid,'\n物理敏感度局部精修也进入候选比较。两种协同方案在相同更紧容差下复算：全满功率能耗%.6f J，选定分配能耗%.6f J，相差%.6f J。\n', ...
        verified{3}.energyJ,verified{2}.energyJ,verified{3}.energyJ-verified{2}.energyJ);
    fprintf(fid,'\n模型保留队友的独立端板与共享双极板热容假设，未重新标定冰参数。温度、电压和冰轨迹保存在同名MAT文件中。\n');
    fprintf(fid,'\n数值事件沿用现有模型的0.01 K过零裕量，因此表内时间对应五片均达到约0.01 ℃，没有追加带载阶段。复核仍为原空间网格，仅将RelTol收紧至1e-6、AbsTol至1e-10、MaxStep至0.0125 s。\n');
    fprintf(fid,'\n纯预热停热时刻另按保存的首次过温事件精修，旧脚本的0.01秒额外时间缓冲已缩小。连续加热协同候选超过成功时刻的冗余计划加热尾段也被缩短并独立复算；实际能耗始终计到成功或关热先到者。\n');
    fprintf(fid,'\n未判定的求解器失败不能用其部分积分能耗参与可行解排名；小幅分配收益也不能解释为具有全局最优性或工程参数不确定性下的保证。\n');
    multistartPath=fullfile(outDir,'q3_multistart_result.mat');
    if exist(multistartPath,'file')
        multi=load(multistartPath,'result');
        fprintf(fid,'\n本轮另用两个不同初始种群进行六变量差分进化，新增%d个变异候选评价（不含复用的16个初始记录，变异候选可能命中缓存）。每组8个个体、2个完整变异代；明细见q3_multistart_history.csv。该搜索规模仍不能证明全局最优。\n',multi.result.newTrialCount);
    end
    cutoffPath=fullfile(outDir,'q3_multistart_cutoff.mat');
    if exist(cutoffPath,'file')
        cutoff=load(cutoffPath,'run');
        fprintf(fid,'\n另对%d个不同可行功率分配检验较低能量预算下的共同停热时刻，完成%d次候选评价，停热后保留原电流曲线和完整成功条件。明细见q3_multistart_cutoff_history.csv。\n', ...
            size(cutoff.run.Q,1),cutoff.run.evaluations);
    end
    disp(summary);
    if ispc
        script=fullfile(codeDir,'export_q3_table4.ps1');
        [status,message]=system(sprintf( ...
            'powershell -NoProfile -ExecutionPolicy Bypass -File "%s"',script));
        assert(status==0,'Table 4 export failed: %s',message);
        fprintf('%s',message);
    end
end
