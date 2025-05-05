
function compareVariabilityParams_FULL()
    % Housekeeping -------------------------------------------------------
    clear; clc; close all;
    fprintf('\n[ compareVariabilityParams_FULL ]  started at %s\n', datestr(now));

    % Result container ---------------------------------------------------
    AllParamResults = struct();

    % Sweep definitions --------------------------------------------------
    sweepDefs = { ...
        {'vthsigma',   0:0.008:0.10, 'V_{th} variability (σ/μ)'}, ... % 14 pts
        {'tausigma',   0:0.020:0.50, 'τ variability (σ/μ)'     }, ... % 36 pts
        {'refpersigma',0:0.020:0.50, 'RefPer variability (σ/μ)'}  };   % 36 pts

    % Loop over parameters ----------------------------------------------
    for p = 1:numel(sweepDefs)
        pName = sweepDefs{p}{1};
        xVals = sweepDefs{p}{2};
        label = sweepDefs{p}{3};

        fprintf('\n--- Sweep: %s  (%d values) ---\n', pName, numel(xVals));

        % Run reservoir for each σ/μ
        [acc, mf, sf, firstDrop] = runParamSweep(pName, xVals);

        % Store in struct
        AllParamResults.(pName) = struct('values',   xVals, ...
                                        'accuracy', acc,   ...
                                        'meanFreq', mf,    ...
                                        'stdFreq',  sf);

        % Generate diagnostic plots
        makeAllPlots(label, xVals, acc, mf, sf, firstDrop);
    end

    save('AllParamResults_v5_pptfigs.mat','AllParamResults','-v7');
    saveAllFigures('fig_sweeps');
    fprintf('\n✔ Saved results -> AllParamResults_FULL.mat');
    fprintf('\n✔ Exported figures -> fig_sweeps/\n');
end

function [accArr, meanF, stdF, dropIdx] = runParamSweep(paramName, grid)
    n = numel(grid);
    accArr = zeros(1,n); meanF = accArr; stdF = accArr;

    for i = 1:n
        sigma_by_mu = grid(i);
        fprintf('  %s = %.3f  -->  ', paramName, sigma_by_mu);
        [acc, mf, sf] = runModelOneParam(paramName, sigma_by_mu);
        accArr(i) = acc; meanF(i) = mf; stdF(i) = sf;
        fprintf('Acc %.1f %% | 〈f〉 %.2f Hz | σ_f %.2f Hz\n', acc, mf, sf);
    end

    dropIdx = find(accArr < 90, 1, 'first');
end

% ======================================================================
%  MAIN PLOTTING ROUTINE  (calls 5 sub‑plots)
% ======================================================================
function makeAllPlots(lbl, x, acc, mf, sf, dropIdx)
    plotAccuracy(lbl,x,acc,dropIdx);
    plotFiringStats(lbl,x,mf,sf);
    plotAccuracyDerivative(lbl,x,acc);
    plotCV(lbl,x,mf,sf);
    plotLogAccuracy(lbl,x,acc);
end

% ----------------------------------------------------------------------
function plotAccuracy(lbl,x,acc,idx)
    figure('Name',[lbl ' – Accuracy'],'NumberTitle','off');
    plot(x,acc,'-o','LineWidth',1.4); hold on; grid on;
    yline(90,'r--','90 %','LineWidth',1.2);
    if ~isempty(idx)
        xline(x(idx),'k--',sprintf('↓ %.3f',x(idx)),'LineWidth',1);
    end
    xlabel(lbl); ylabel('Accuracy (%)');
    title(['Accuracy vs ', lbl]);
end

% ----------------------------------------------------------------------
function plotFiringStats(lbl,x,mf,sf)
    figure('Name',[lbl ' – Firing stats'],'NumberTitle','off');
    plot(x,mf,'-o','LineWidth',1.3); hold on;
    plot(x,sf,'-s','LineWidth',1.3);
    legend('\langle f \rangle','\sigma_f','Location','best'); grid on;
    xlabel(lbl); ylabel('Hz');
    title({['Firing‑rate statistics vs ', lbl], ...
           '\fontsize{8}〈f〉 = (Σ spikes)/(N·T),   σ_f = std(f_i)'});
end

% ----------------------------------------------------------------------
function plotAccuracyDerivative(lbl,x,acc)
    dAcc = diff(acc) ./ diff(x);
    figure('Name',[lbl ' – Accuracy & dAcc'],'NumberTitle','off');
    yyaxis left;  plot(x,acc,'-o','LineWidth',1.2); ylabel('Accuracy (%)');
    yyaxis right; plot(x(2:end),dAcc,'--s','LineWidth',1.2); ylabel('ΔAcc/Δσ');
    grid on; xlabel(lbl); title(['Accuracy & derivative vs ', lbl]);
end

% ----------------------------------------------------------------------
function plotCV(lbl,x,mf,sf)
    CV = sf ./ mf;
    figure('Name',[lbl ' – CV'],'NumberTitle','off');
    plot(x,CV,'-^','LineWidth',1.3); grid on;
    xlabel(lbl); ylabel('CV = σ_f / 〈f〉');
    title(['Coefficient of variation vs ', lbl]);
end

% ----------------------------------------------------------------------
function plotLogAccuracy(lbl,x,acc)
    figure('Name',[lbl ' – log10'],'NumberTitle','off');
    semilogx(x(x>0), acc(x>0), '-o','LineWidth',1.3); grid on;
    xlabel('log_{10}(σ/μ)'); ylabel('Accuracy (%)');
    title(['Accuracy vs log_{10} variability — ', lbl]);
end


function saveAllFigures(folder)
    if ~exist(folder,'dir'), mkdir(folder); end
    figs = findall(0,'Type','figure');
    for f = reshape(figs,1,[])
        name = get(f,'Name'); if isempty(name), name = ['Fig' num2str(f.Number)]; end
        name = regexprep(name,'[\\/:*?"<>| ]','_');
        saveas(f, fullfile(folder,[name '.png']));
    end
end

% ======================================================================
%  RESERVOIR SIMULATION 
% ======================================================================
function [meanAcc, meanFreq, stdFreq] = runModelOneParam(paramName, sigma_by_mu)
    % DATA ---------------------------------------------------------------
    load('preprocessing.mat','DATA');
    load('connections.mat','Gin','Gres');
    load('trained_weights.mat','trained_weights');
    Nin  = size(Gin,2); Nres = size(Gin,1);

    % CONSTANTS ----------------------------------------------------------
    dt=1e-3; tau_mean=64e-3; Vrst=0;
    RefPer_mean=2e-3; RPS_mean=ceil(RefPer_mean/dt); Vth_mean=20;
    tau1=8e-3; tau2=4e-3; ds=1e-3; alpha_Gin=8; alpha_Gres=1;

    % Variability --------------------------------------------------------
    useVarTh=0; useVarTau=0; useVarRef=0;
    Vth_sigma=0; tau_sigma=0; Ref_sigma=0;
    switch lower(paramName)
        case 'vthsigma',    useVarTh=1;  Vth_sigma = sigma_by_mu*Vth_mean;
        case 'tausigma',    useVarTau=1; tau_sigma = sigma_by_mu*tau_mean;
        case 'refpersigma', useVarRef=1; Ref_sigma = sigma_by_mu*RefPer_mean;
    end

    % Synaptic impulse ---------------------------------------------------
    I0=1/(tau1-tau2); TS=0:dt:ds+6*(tau1+tau2);
    Iwave=I0*(exp(-(TS-ds)/tau1)-exp(-(TS-ds)/tau2));
    Iwave(1:max(ceil(ds/dt),1))=0;
    Gnet=[alpha_Gin*Gin, alpha_Gres*Gres'];

    % Static per-neuron params ------------------------------------------
    if useVarTau,  rng(2); tauVec=max(tau_mean+tau_sigma*randn(Nres,1),1e-6);
    else, tauVec=tau_mean*ones(Nres,1); end
    if useVarRef,  rng(3); RefVec=max(RefPer_mean+Ref_sigma*randn(Nres,1),0); RPSVec=ceil(RefVec/dt);
    else, RefVec=RefPer_mean*ones(Nres,1); RPSVec=RPS_mean*ones(Nres,1); end

    % Simulation ---------------------------------------------------------
    parfor s=1:numel(DATA)
        S_in = DATA(s).S; T = size(S_in,2);
        if useVarTh, rng(1); VthVec=max(Vth_mean+Vth_sigma*randn(Nres,1),0.1);
        else, VthVec=Vth_mean*ones(Nres,1); end
        V=zeros(Nres,1); Ibuf=zeros(Nin+Nres,numel(Iwave)); RP=zeros(Nres,1);
        res_prev=zeros(Nres,1); RES=zeros(Nres,T);
        for t=1:T
            spikes=[S_in(:,t); res_prev];
            Ibuf=circshift(Ibuf,-1,2); Ibuf(:,end)=0; Ibuf=Ibuf+spikes*Iwave; Itot=Ibuf(:,1);
            RP=max(RP-1,0);
            V=V.*(1-dt./tauVec) + dt*(Gnet*Itot); V(RP>0)=0;
            res_prev=V>VthVec;
            V(res_prev)=Vrst; V(V<0)=Vrst; RP(res_prev)=RPSVec(res_prev);
            RES(:,t)=res_prev;
        end
        DATA(s).RES = RES;
    end

    % Readout ------------------------------------------------------------
    rZ=[]; for s=1:numel(DATA), rZ(:,s)=sum(DATA(s).RES,2); end
    folds=5; stride=numel(DATA)/folds; ACC=zeros(1,folds);
    for k=1:folds
        test=(k-1)*stride+(1:stride);
        W=trained_weights(k).W; [~,pred]=max(W*rZ,[],1); true=[DATA.type]+1;
        ACC(k)=100*sum(pred(test)==true(test))/numel(test);
    end
    meanAcc=mean(ACC);

    % Firing stats -------------------------------------------------------
    allSpk=cat(2,DATA.RES); [~,freq]=calcMeanDevSpkFreq(allSpk,dt);
    meanFreq=mean(freq); stdFreq=std(freq);
end

% ======================================================================
function [meanDev,freq] = calcMeanDevSpkFreq(resSpikes,dt)
    total=sum(resSpikes,2); freq=total/(size(resSpikes,2)*dt);
    meanDev=mean(abs(freq-mean(freq)));
end
