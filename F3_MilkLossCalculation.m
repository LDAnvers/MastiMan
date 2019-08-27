function [OutData] = F3_MilkLossCalculation(data, lme, Level, CaseData, fig)
% Data contains the milkings data with model predictions and residuals
% INPUTS:   data                contains milkings data of this farm,
%                               merged, but not sorted out (only raw data).
%                               Not sorted in one column for QMY data
%           lme                 contains trained linear mixed model of farm
%           Level               predictions at quarter or at udder level
%           CaseData            AnimalID, lactation, DIM detection, days to
%                               calculate losses for, Quarter (if level=1)
%           figure              plot figure of model and losses: upper
%                               panel= data + model; lower panel =
%                               residuals of 1 or 4 quarters
%
% OUTPUTS:  OutData             Milkings data, containing prediction and
%                               residuals for all quarters or the whole
%                               udder during X days, as indicated in
%                               CaseData
%           h                   figure handle
%
% STEP 1: Check inputs in CaseData
% STEP 2: Prepare FIX, Xp, Zp, prepare input data, e.g. delete meas 
%         MI>24/<6 and MY < 0.25 and calculate exppMI, logMI, etc
% STEP 3: Calculate losses using bayesian prediction, using Resi and bi,
%         input data/if DIM detection < 5 use the marginal model, 
%         otherwise use the data of the specific case and quarters.
% STEP 3: plot results if figure == 1

%% STEP 1: Check inputs in CaseData



% data = InData;              % raw merged datasets
% CaseData = [626 3 90 2; 42 3 160 4];
% lme = lme_TMY;
% Level = 2;
% fig = 1;




%% Predict reference MY for mastitis cases
% select appropriate data / cows
% put data and cows in the right format: calculate pMI, Cat, Qp, LP, ...
% Use the estimates of the mixed model to predict the MY during the
% mastitis
% Formula: Wood_predicted_i = Xp*FIX + Zp*bi
%       Wood_predicted_i = the predicted quarter milk yield at time i
%       Xp  =   Fixed effects predictors at timepoint i
%           =   [1 DIM_i logMI_i log DIM_i logMI_i.*DIM_i]
%       FIX =   Fixed effects estimates of the mixed model (beta) for that
%               cow in that category (depending on QP and lac)
%       Zp  =   Random effects predictors at timepoint i
%           =   [1 DIM_i logMI_i logDIM_i]
%       bi  =   Random effect estimate of the mixed model for that cow 
%               in that category (depending on QP and Lac)
%           =   D*Zi'*Wi*Resi
%               D = variance-covariance matrix of the mixed model
%               Zi = [ones(ni,1) DIM logMI logDIM] of all the data included
%                       for that cow/quarter lactation, ni = length ref set
%               Resi = logMY - Xi*FIX = residuals of fixed effects
%               Wi = (Zi*D*Zi'+Epsi)^(-1);
%               Epsi = Vare*eye(ni); Vare = residual variance

% prepare output data
if Level == 1
    data.QMYp1(:,1) = NaN;         % fill in QMY for output
    data.QMYp2(:,1) = NaN;         % fill in QMY for output
    data.QMYp3(:,1) = NaN;         % fill in QMY for output
    data.QMYp4(:,1) = NaN;         % fill in QMY for output
else
    data.TMYp(:,1) = NaN;         % fill in TMY for output
end


% pre-detection: 5 days
[D, Vare] = covarianceParameters(lme); D = D{:,:};    % obtain D, the Covariancematrix and the residual variance
beta = lme.Coefficients.Estimate(:,:);                % obtain beta, the estimates of the fixed effects


for i = 1:length(CaseData(:,1))             % for all cases in CaseData
    if Level == 1                           % if we do calculations at Ql
        ind = find(data.BA == CaseData(i,1) & data.Lac == CaseData(i,2));
        subset = array2table([data.BA(ind) data.Lac(ind) data.TMY(ind) ...
            data.MYLF(ind) data.MYRF(ind) data.MYLH(ind) data.MYRH(ind)...
            data.DIM(ind) data.MI(ind)],...
            'VariableNames',{'BA','Lac','TMY','MYLF','MYRF','MYLH','MYRH','DIM','MI'});
        
        subset.pMI(1,1) = NaN;                    % previous milking interval
        subset.pMI(2:end,1) = subset.MI(1:end-1); % previous milking interval
        subset.logMI(:,1) = log(subset.MI);       % milking interval
        subset.logDIM(:,1) = log(subset.DIM);     % logDIM
        
        for j = 1:4                               % calculate per quarter
            
            idx = find(subset.DIM < CaseData(i,3)-5 & isnan(subset.logMI)==0 & ...
                subset.MI > 4 & subset.pMI > 4 & subset.MI < 24 & subset.pMI < 24 & ...
                subset{:,3+j} > 0.25);   % find casedata before detection date
            
            % Prepare reference data for calculation of bi
            n = length(idx);
            Xref = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)  subset.logMI(idx).*subset.DIM(idx)];
            Zref = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)];
            
            % Set up fixed effects coefficients
            if subset.Lac(1) == 1        % cat = 1; QP = 1; L = 1
                if j < 3
                    FIX = [beta(1)                          beta(2)         beta(3)     beta(6)           beta(7)]';
                else
                    FIX = [beta(1)+beta(5)                  beta(2)         beta(3)     beta(6)+beta(10)  beta(7)]';
                end
            else
                if j < 3
                    FIX = [beta(1)+beta(4)                  beta(2)+beta(8) beta(3)     beta(6)           beta(7)]';
                else                                        % cat = 4; QP = 2; L = 2+
                    FIX = [beta(1)+beta(4)+beta(5)+beta(9)  beta(2)+beta(8) beta(3)     beta(6)+beta(10)  beta(7)]';
                end
            end
            
            if j == 1
                logMY = log(subset.MYLF(idx)*1000);           % log MY first Q
            elseif j==2
                logMY = log(subset.MYRF(idx)*1000);           % log MY second Q
            elseif j==3
                logMY = log(subset.MYLH(idx)*1000);           % log MY third Q
            else
                logMY = log(subset.MYRH(idx)*1000);           % log MY fourth Q
            end
            
            % Calculate residuals of fixed effects part of the model for
            % reference period
            Resi = logMY - Xref*FIX;
            
            % Calculate W for reference period
            Epsi = Vare*eye(n);
            W = (Zref*D*Zref'+Epsi)^(-1);
            
            % Calculate bi
            b = D*Zref'*W*Resi;
            
            % Calculate model prediction
            WOODref = Xref*FIX + Zref*b;        % est logMY
            WOODrefR = exp(WOODref)/1000;       % est milk yield
            
            % Add to subset (later: fill in to data)
            if j == 1
                subset.QMY_1p(:,1) = NaN;
                subset.QMY_1p(idx,1) = WOODrefR;         % log MY first Q
            elseif j==2
                subset.QMY_2p(:,1) = NaN;
                subset.QMY_2p(idx,1) = WOODrefR;           % log MY second Q
            elseif j==3
                subset.QMY_3p(:,1) = NaN;
                subset.QMY_3p(idx,1) = WOODrefR;          % log MY third Q
            else
                subset.QMY_4p(:,1) = NaN;
                subset.QMY_4p(idx,1) = WOODrefR;           % log MY fourth Q
            end
            
            % Estimate logMY during mastitis: find meas in ref period
            idx = find(subset.DIM > CaseData(i,3)-5 & subset.DIM < CaseData(i,3)+21 & isnan(subset.logMI)==0) ;% find casedata after detection date-5
            
            % Prepare input of mastitis period
            n = length(idx);
            Xmas = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)  subset.logMI(idx).*subset.DIM(idx)];
            Zmas = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)];
            
            % Predict milk yield in reference period
            WOODmas = Xmas*FIX + Zmas*b;
            WOODmasR = exp(WOODmas)/1000;
            
            % Add to subset (later: fill in to data)
            if j == 1
                subset.QMY_1p(idx,1) = WOODmasR;         % log MY first Q
            elseif j==2
                subset.QMY_2p(idx,1) = WOODmasR;         % log MY second Q
            elseif j==3
                subset.QMY_3p(idx,1) = WOODmasR;         % log MY third Q
            else
                subset.QMY_4p(idx,1) = WOODmasR;         % log MY fourth Q
            end
        end
        
        data.QMYp1(ind,1) = subset.QMY_1p;         % fill in QMY for output
        data.QMYp2(ind,1) = subset.QMY_2p;         % fill in QMY for output
        data.QMYp3(ind,1) = subset.QMY_3p;         % fill in QMY for output
        data.QMYp4(ind,1) = subset.QMY_4p;         % fill in QMY for output
        h = [];
        
        if fig == 1
            h = figure; hold on;
%             COL1 = [0.4 0.698 1; 0.4 0.698 1; 0.4 0.698 1; 0.4 0.698 1]; COL1(CaseData(i,4),:) = [75/255 0 130/255];
            COL2 = [119/255 136/255 153/255;119/255 136/255 153/255; 119/255 136/255 153/255; 119/255 136/255 153/255]; COL2(CaseData(i,4),:) = [75/255 0 130/255];
            
            VN1 = find(string(data.Properties.VariableNames) == 'MYLF');    %MY
            VN2 = find(string(data.Properties.VariableNames) == 'QMYp1');    %predicted MY
            
            for j = 1:4
                subplot(2,4,j); hold on;
                plot(data.DIM(ind),data{ind,VN1-1+j},'.-','LineWidth',1,'Color',[0.4 0.698 1],'MarkerSize',12)
                plot(data.DIM(ind),data{ind,VN2-1+j},'.-','LineWidth',1,'Color',COL2(j,:),'MarkerSize',12)
                plot([CaseData(i,3) CaseData(i,3)],[0 10],'Color',[76/255 153/255 0],'LineWidth',2)
                xlim([0 CaseData(i,3)+21]); ylim([0 10]); box on; xlabel('DIM (days)'); ylabel('QMY (kg)'); title(['MY Q ' num2str(j) ' Case ' num2str(i)])
                subplot(2,4,4+j); hold on;
                plot(data.DIM(ind),zeros(length(ind),1),'k:')
                plot(data.DIM(ind),data{ind,VN2-1+j}-data{ind,VN1-1+j},'.-','LineWidth',1,'Color',COL2(j,:),'MarkerSize',12)
                plot([CaseData(i,3) CaseData(i,3)],[-5 5],'Color',[76/255 153/255 0],'LineWidth',2)
                xlim([0 CaseData(i,3)+21]); ylim([-5 5]); box on; xlabel('DIM (days)'); ylabel('Residual QMY (kg)'); title(['Residual MY Q ' num2str(j) ' Case ' num2str(i)])
            end
        end
    elseif Level == 2                   % calculations at udder level
        ind = find(data.BA == CaseData(i,1) & data.Lac == CaseData(i,2));
        subset = array2table([data.BA(ind) data.Lac(ind) data.TMY(ind) ...
            data.DIM(ind) data.MI(ind)],...
            'VariableNames',{'BA','Lac','TMY','DIM','MI'});
        
        subset.pMI(1,1) = NaN;                    % previous milking interval
        subset.pMI(2:end,1) = subset.MI(1:end-1); % previous milking interval
        subset.logMI(:,1) = log(subset.MI);       % milking interval
        subset.logDIM(:,1) = log(subset.DIM);     % logDIM
        
        
        idx = find(subset.DIM < CaseData(i,3)-5 & isnan(subset.logMI)==0 & ...
            subset.MI > 4 & subset.pMI > 4 & subset.MI < 24 & subset.pMI < 24 & ...
            subset.TMY > 0.5);         % find casedata before detection date
        
        % Prepare reference data for calculation of bi
        n = length(idx);
        Xref = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)  subset.logMI(idx).*subset.DIM(idx)];
        Zref = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)];
        
        % Set up fixed effects coefficients
        if subset.Lac(1) == 1        % L = 1
            FIX = [beta(1)         beta(2)          beta(3)   beta(5)  beta(6)]';
        else
            FIX = [beta(1)+beta(4) beta(2)*beta(7)  beta(3)   beta(5)  beta(6)]';
        end
        
        logMY = log(subset.TMY(idx)*1000);           % log MY first Q

        % Calculate residuals of fixed effects part of the model for
        % reference period
        Resi = logMY - Xref*FIX;
        
        % Calculate W for reference period
        Epsi = Vare*eye(n);
        W = (Zref*D*Zref'+Epsi)^(-1);
        
        % Calculate bi
        b = D*Zref'*W*Resi;
        
        % Calculate model prediction
        WOODref = Xref*FIX + Zref*b;        % est logMY
        WOODrefR = exp(WOODref)/1000;       % est milk yield
        
        % add to subset (later: fill in to data)
        subset.TMY_p(:,1) = NaN;
        subset.TMY_p(idx,1) = WOODrefR;

        % Estimate logMY during mastitis: find meas in ref period
        idx = find(subset.DIM > CaseData(i,3)-5 & subset.DIM < CaseData(i,3)+21 & isnan(subset.logMI)==0) ;% find casedata after detection date-5
        
        % Prepare input of mastitis period
        n = length(idx);
        Xmas = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)  subset.logMI(idx).*subset.DIM(idx)];
        Zmas = [ones(n,1) subset.DIM(idx) subset.logMI(idx) subset.logDIM(idx)];
        
        % Predict milk yield in reference period
        WOODmas = Xmas*FIX + Zmas*b;
        WOODmasR = exp(WOODmas)/1000;
        
        % Add to subset (later: fill in to data)
        subset.TMY_p(idx,1) = WOODmasR;         % log MY first Q
    
        data.TMYp(ind,1) = subset.TMY_p;         % fill in TMY for output
        h = [];
    
        if fig == 1
            h = figure; hold on;
                       
            subplot(1,2,1); hold on;
            plot(data.DIM(ind),data.TMY(ind),'.-','LineWidth',1,'Color',[0.4 0.698 1],'MarkerSize',12)
            plot(data.DIM(ind),data.TMYp(ind),'.-','LineWidth',1,'Color',[128/255 0 128/255],'MarkerSize',12)
            plot([CaseData(i,3) CaseData(i,3)],[0 35],'Color',[76/255 153/255 0],'LineWidth',2)
            xlim([0 CaseData(i,3)+21]); ylim([0 35]); box on; xlabel('DIM (days)'); ylabel('MY (kg)'); title(['TMY Case ' num2str(i)])
            subplot(1,2,2); hold on;
            plot(data.DIM(ind),zeros(length(ind),1),'k:')
            plot(data.DIM(ind),data.TMY(ind)-data.TMYp(ind),'.-','LineWidth',1,'Color',[128/255 0 128/255],'MarkerSize',12)
            plot([CaseData(i,3) CaseData(i,3)],[-10 10],'Color',[76/255 153/255 0],'LineWidth',2)
            xlim([0 CaseData(i,3)+21]); ylim([-10 10]); box on; xlabel('DIM (days)'); ylabel('Residual MY (kg)'); title(['Residual TMY Case ' num2str(i)])
        end
    end
end

% Output data
OutData = data;
