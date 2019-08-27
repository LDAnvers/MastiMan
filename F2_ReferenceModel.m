function [lme, GOF, data, MPE, MAE] = F2_ReferenceModel(data,Level,FarmN, ValOut)
% This function will calculate the prior distribution of a specific farm
% for the reference modelling, based on a mixed model
% 
% INPUTS:  data          - dataset to estimate model parameters on,
%                          containing at least the following variables:
%                               * T         = counter / within-data ID
%                                               ~Nominal
%                               * DIM       = days in milk
%                                               ~Numeric
%                               * logDIM    = logarithm of days in milk
%                                               ~Numeric
%                               * QMY / TMY = milk yield
%                                               ~Numeric
%                               * logMY     = log (MY * 1000)
%                               * LP
%          Level         - QMY = 1 or TMY = 2
%          FarmN         - For Summary and GOF
%
%
% OUTPUTS: lme           - model parameters of the mixed model
%                          for QMY(1) or TMY (2)
%          GOF           - Goodness of fit measures of the fitted linear
%                          mixed model
%          data          - final dataset to fit the mixed model, including
%                          the model predictions and the residuals; 
%          ValData       - if ValOut = 1, then data contains also the
%                          results of the crossvalidation using Bayesian
%                          predictions, at day
%                          10,15,20,30,40,65,115,165,215,265 with 2/3 input
%                          data to train the model and 1/3 test set data
%                               
% STEP 0: prepare data so that extreme values (MI/MY) are not affecting
%         model fit (ensure robust estimation of prior distribution)
% STEP 1: prepare modelling: set up models and matrices to generate output
% STEP 2: fit model and calculate model performance + generate outputs
% STEP 3: crossvalidation

% data = D1_TMY;
% Level = 2;
% FarmN = 1;

%% STEP 0: prepare data
% put extra restrictions on data:
%   DIM < 305 (effect preg: 8 months + +/- 65 days to conceive)
%   MI < 24 & > 4h
%   pMI < 24 & > 4h
%   QMY > 0.25 L | TMY > 1 L

GOF = array2table([FarmN length(unique(data.T)) length(data.T) 0 0 0],'VariableNames',{'FarmN','NLac','Nmeas','NmeasSel', 'NmeasProc','NDIM305'});


if Level == 1
    GOF.NDIM305 = length(find(data.DIM > 305)); % milkings deleted > 305 DIM
    ind = find(data.QMY > 0.5 & data.DIM < 305 & isnan(data.QMY)==0 & ...
            data.MI > 4 & data.pMI > 4 & data.MI < 24 & data.pMI < 24 & ...
            isnan(data.logMI) == 0);     % select lactations
       
    data = data(ind,:);                 % Number of measurements selected
    GOF.NmeasSel = length(ind);         % Fill in how many selected
    GOF.NmeasProc = length(ind)/GOF.Nmeas;  % fill in %
else
    GOF.NDIM305 = length(find(data.DIM > 305)); % milkings deleted > 305 DIM
    ind = find(data.TMY > 2 & data.DIM < 305 & isnan(data.TMY)==0 & ...
            data.MI > 4 & data.pMI > 4 & data.MI < 24 & data.pMI < 24 & ...
            isnan(data.logMI) == 0);     % select lactations
       
    data = data(ind,:);                 % Number of measurements selected
    GOF.NmeasSel = length(ind);         % Fill in how many selected
    GOF.NmeasProc = length(ind)/GOF.Nmeas;  % fill in %
end

%% STEP 1: set up model and prepare outputs
% define model depending on Level

if Level == 1
    mod = 'logMY ~ 1+ DIM + logDIM + logMI  + LP*QP + DIM:LP + logDIM:QP + logMI:DIM + (1 + DIM + logDIM + logMI | T)';
else
    mod = 'logMY ~ 1+ DIM + logDIM + logMI  + LP + DIM:LP + logMI:DIM + (1 + DIM + logDIM + logMI | T)';     % TBC that this is the best model for udder level!!!
end


%% STEP 2: fit model and evaluate model performance
% fit model
if Level == 1
    data.LP = categorical(double(data.LP));
    data.QP = categorical(double(data.QP));
    data.T = categorical(double(data.T));

    data.logDIM(data.logDIM < 0) = 0;
    data.DIM(data.DIM < 0) = 0;
    lme = fitlme(data,mod,'FitMethod','ML');  % fit linear mixed model

    GOF.LL = lme.LogLikelihood;
    GOF.AIC = lme.ModelCriterion{1,1};
    GOF.Rsq = lme.Rsquared.Ordinary;

%     D = covarianceParameters(lme); D = D{:};
%     F = lme.Coefficients.Estimate(:,:);  % IC / DIM / logMI / LP / logDIM / DIM:logMI / DIM:QP / LP:QP:logDIM

    data.predMY(:,1) = exp(predict(lme))/1000;      % fill in QMY in data
    data.resMY(:,1) = data.QMY - data.predMY;       % fill in residual QMY
    data.resMYproc = abs(data.resMY)./data.QMY*100; % procentual error

    GOF.RMSE = sum(data.resMY.^2)./length(data.resMY); % mean RMSE
    GOF.MPE = mean(data.resMYproc);                 % mean predicitve error
else
    data.T = categorical(double(data.T));
    data.LP = categorical(double(data.LP));
    
    data.logDIM(data.logDIM < 0) = 0;
    data.DIM(data.DIM < 0) = 0;
    lme = fitlme(data,mod,'FitMethod','ML');  % fit linear mixed model

    GOF.LL = lme.LogLikelihood;
    GOF.AIC = lme.ModelCriterion{1,1};
    GOF.Rsq = lme.Rsquared.Ordinary;

%     D = covarianceParameters(lme); D = D{:};
%     F = lme.Coefficients.Estimate(:,:);  % IC / DIM / logMI / LP / logDIM / DIM:logMI / DIM:LP

    data.predMY(:,1) = exp(predict(lme))/1000;      % fill in QMY in data
    data.resMY(:,1) = data.TMY - data.predMY;       % fill in residual QMY
    data.resMYproc = abs(data.resMY)./data.TMY*100; % procentual error

    GOF.RMSE = sum(data.resMY.^2)./length(data.resMY); % mean RMSE
    GOF.MPE = mean(data.resMYproc);                 % mean predicitve error
    
end


% plot figures to check
% for i = randi(max(double(data.T),1,5))
%    figure; hold on;
%    plot(data.DIM(double(data.T) == i), data.TMY(double(data.T)==i),'k.-','MarkerSize',12)
%    plot(data.DIM(double(data.T) == i), data.predMY(double(data.T)==i),'r.-','MarkerSize',15)
%     
% end

%% STEP 3: Cross-validation for the input dataset
% In this step, we'll do a 10x 1/3-2/3 cross validation to check model
% performance on part of the data, and using Bayes techniques to predict
% the milk yield 50 days ahead.

if ValOut == 1
    
    % Set up the datasets for the cross-validation
    % first reset random number generator and generate random nr of datapoints
    rng(0);
    T = (1:max(double(data.T)))';
    testN = zeros(round(1/3*max(T)),10);                            % test
    trainN = zeros(max(T)-round(1/3*max(T)),10);                    % training
    for i = 1:10
        testN(:,i) = sort(randperm(max(T),round(1/3*max(T))))';     % T for testset
        trainN(:,i) = T(ismember(T, testN(:,i))==0);                % T for training set
    end
    
    % Select subsets from 'data' to fit and test model
    for i = 1:10
        C = sprintf('set_%d',i);
        
        training.(C) = data(ismember(double(data.T),trainN(:,i)) == 1,:);    % training
        test.(C) = data(ismember(double(data.T),testN(:,i)) == 1,:);    % test
    end
    
    % fit models and save model values
    for i = 1:10
        C = sprintf('set_%d',i);
        ValMod.(C) = fitlme(training.(C),mod,'FitMethod','ML');  % fit linear mixed model
    end
    
    % define the input data (number of days taken in)
    INdata = [10 15 20 30 40 65 115 165 215 265];% days input for prediction
    PREDdata = 0:49;                             % days to predict
    M = length(PREDdata);                        % reference period and prolonged ref period
    
    if Level == 1
        tic
        for i = 1:10
            
            C = sprintf('set_%d',i);                 % define the set we're working in
            
            % define the prior distribution parameters
            [D, Vare] = covarianceParameters(ValMod.(C)); D = D{:,:};    % obtain D, the Covariancematrix and the residual variance
            beta = ValMod.(C).Coefficients.Estimate(:,:);                % obtain beta, the estimates of the fixed effects
            
            QL = sort(unique(double(test.(C).T)));    % define all subjects to predict
            
            for j = 1:length(QL)                      % select the QL to make the predictions for
                for k = 1:length(INdata)              % select the number of days to include in the est of bi
                    
                    ind = find(double(test.(C).T) == QL(j) & test.(C).DIM <= INdata(k));    % find all the data before INdata to train model
                    
                    % prepare input of reference period:  IC / DIM / logMI / LP / logDIM / DIM:logMI / DIM:LP
                    n = length(ind);
                    Xref = [ones(n,1) test.(C).DIM(ind) test.(C).logMI(ind) test.(C).logDIM(ind) test.(C).logMI(ind).*test.(C).DIM(ind)];
                    Zref = [ones(n,1) test.(C).DIM(ind) test.(C).logMI(ind) test.(C).logDIM(ind)];
                    
                    if double(test.(C).Cat(ind(1))) == 1        % cat = 1; QP = 1; L = 1
                        FIX = [beta(1)                          beta(2)         beta(3)     beta(6)           beta(7)]';
                    elseif double(test.(C).Cat(ind(1))) == 2    % cat = 2; QP = 2; L = 1
                        FIX = [beta(1)+beta(5)                  beta(2)         beta(3)     beta(6)+beta(10)  beta(7)]';
                    elseif double(test.(C).Cat(ind(1))) == 3    % cat = 3; QP = 1; L = 2+
                        FIX = [beta(1)+beta(4)                  beta(2)+beta(8) beta(3)     beta(6)           beta(7)]';
                    else                                        % cat = 4; QP = 2; L = 2+
                        FIX = [beta(1)+beta(4)+beta(5)+beta(9)  beta(2)+beta(8) beta(3)     beta(6)+beta(10)  beta(7)]';
                    end
                    
                    % calculate residuals of fixed effects part of the model for
                    % reference period
                    Resi = test.(C).logMY(ind) - Xref*FIX;
                    
                    % Calculate W for reference period
                    Epsi = Vare*eye(n);
                    W = (Zref*D*Zref'+Epsi)^(-1);
                    
                    % Calculate bi
                    b = D*Zref'*W*Resi;
                    
                    % predict INdata(k)+ PREDdata
                    for l = 1:M
                        
                        % find all data >= INdata(k)+PREDdata(l) and < Indata(k) +
                        % PREDdata(l)
                        idx = find(double(test.(C).T) == QL(j) & test.(C).DIM >= INdata(k)+PREDdata(l) & test.(C).DIM < INdata(k)+PREDdata(l)+1 & isnan(test.(C).logMI)==0);
                        
                        % find data we want to predict
                        n = length(idx);
                        Xref = [ones(n,1) test.(C).DIM(idx) test.(C).logMI(idx) test.(C).logDIM(idx) test.(C).logMI(idx).*test.(C).DIM(idx)];
                        Zref = [ones(n,1) test.(C).DIM(idx) test.(C).logMI(idx) test.(C).logDIM(idx)];
                        
                        % predict logMYg
                        MYpred =  exp(Xref*FIX + Zref*b)/1000;      %  MY
                        
                        %                plot(test.(C).DIM(idx),test.(C).QMY(idx),'.-'); hold on;
                        %                plot(test.(C).DIM(idx),MYpred,'.-');
                        
                        MPE.(C)(j,M*k-M+l) = mean(abs((test.(C).QMY(idx)- MYpred)./test.(C).QMY(idx)));    % mean proportional error per day
                        MAE.(C)(j,M*k-M+l) = mean(abs((test.(C).QMY(idx)- MYpred)));
                    end
                end
            end
        end
        clear ans b B beta C D Epsi FIX i idx ind INdata j k l M MYpred n QL Resi Vare W Xref Zref
        toc
    else        % level == 2
        tic
        for i = 1:10
            
            C = sprintf('set_%d',i);                 % define the set we're working in
            
            % define the prior distribution parameters
            [D, Vare] = covarianceParameters(ValMod.(C)); D = D{:,:};    % obtain D, the Covariancematrix and the residual variance
            beta = ValMod.(C).Coefficients.Estimate(:,:);                % obtain beta, the estimates of the fixed effects
            
            QL = sort(unique(double(test.(C).T)));    % define all subjects to predict
            
            for j = 1:length(QL)                      % select the QL to make the predictions for
                for k = 1:length(INdata)              % select the number of days to include in the est of bi
                    
                    ind = find(double(test.(C).T) == QL(j) & test.(C).DIM <= INdata(k));    % find all the data before INdata to train model
                    
                    % prepare input of reference period:  IC / DIM / logMI / LP / logDIM / DIM:logMI / DIM:LP
                    n = length(ind);
                    Xref = [ones(n,1) test.(C).DIM(ind) test.(C).logMI(ind) test.(C).logDIM(ind) test.(C).logMI(ind).*test.(C).DIM(ind)];
                    Zref = [ones(n,1) test.(C).DIM(ind) test.(C).logMI(ind) test.(C).logDIM(ind) ];
                    
                    if test.(C).Lac(ind(1)) == 1        % L = 1
                        FIX = [beta(1)         beta(2)          beta(3)   beta(5)  beta(6)]';
                    else 
                        FIX = [beta(1)+beta(4) beta(2)*beta(7)  beta(3)   beta(5)  beta(6)]';
                    end
                    
                    % calculate residuals of fixed effects part of the model for
                    % reference period
                    Resi = test.(C).logMY(ind) - Xref*FIX;
                    
                    % Calculate W for reference period
                    Epsi = Vare*eye(n);
                    W = (Zref*D*Zref'+Epsi)^(-1);
                    
                    % Calculate bi
                    b = D*Zref'*W*Resi;
                    
                    % predict INdata(k)+ PREDdata
                    for l = 1:M
                        
                        % find all data >= INdata(k)+PREDdata(l) and < Indata(k) + PREDdata(l)
                        idx = find(double(test.(C).T) == QL(j) & test.(C).DIM >= INdata(k)+PREDdata(l) & test.(C).DIM < INdata(k)+PREDdata(l)+1 & isnan(test.(C).logMI)==0);
                        
                        % find data we want to predict
                        n = length(idx);
                        Xref = [ones(n,1) test.(C).DIM(idx) test.(C).logMI(idx) test.(C).logDIM(idx) test.(C).logMI(idx).*test.(C).DIM(idx)];
                        Zref = [ones(n,1) test.(C).DIM(idx) test.(C).logMI(idx) test.(C).logDIM(idx)];
                        
                        % predict logMYg
                        MYpred =  exp(Xref*FIX + Zref*b)/1000;      %  MY
                        
                        %                plot(test.(C).DIM(idx),test.(C).QMY(idx),'.-'); hold on;
                        %                plot(test.(C).DIM(idx),MYpred,'.-');
                        
                        MPE.(C)(j,M*k-M+l) = mean(abs((test.(C).TMY(idx)- MYpred)./test.(C).TMY(idx)));    % mean proportional error per day
                        MAE.(C)(j,M*k-M+l) = mean(abs((test.(C).TMY(idx)- MYpred)));
                    end
                end
            end
        end
        clear ans b B beta C D Epsi FIX i idx ind INdata j k l M MYpred n QL Resi Vare W Xref Zref
        toc
    end
    
    % MPE and MAE consist of 10 tables, each containing n rows with prediction errors
    % 1:50 contain the MPE of day 1 - 50 estimated with 10 days of data
    % 51:100 contain the MPE of day 1 - 50 estimated with 15 days of data
    % 101:150 contain the MPE of day 1 - 50 estimated with 15 days of data.
    
    % We first want to average over all rows / testset = average all subjects
    for i = 1:10
        C = sprintf('set_%d',i);
        
        MPE.means(i,:) = nanmean(MPE.(C));
        MPE.std(i,:) = nanstd(MPE.(C));
    
        MAE.means(i,:) = nanmean(MAE.(C));
        MAE.std(i,:) = nanstd(MAE.(C));
    end
    
    COL = [0 0 1; 0 0 0.5 ;0 0.5 1; 0 1 1; 0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; ...
        0.5 0 0.5; 1 0 1; 1 0 0.5];
    h = figure('OuterPosition',[50 100 1000 500]); hold on
    for i = 1:10
        subplot(1,2,1); hold on; box on
        b(i,1:2) = regress(mean(MPE.means(:,50*i-49:50*i))'.*100,[ones(50,1) (1:50)']);
        plot(1:50,mean(MPE.means(:,50*i-49:50*i)).*100,'.', 'MarkerSize',15,'Color',COL(i,:))
        plot(1:50, b(i,1)+b(i,2)*(1:50),'Color',COL(i,:), 'LineWidth', 2)
        
        subplot(1,2,2); hold on; box on
        b(i,1:2) = regress(mean(MAE.means(:,50*i-49:50*i))',[ones(50,1) (1:50)']);
        plot(1:50,mean(MAE.means(:,50*i-49:50*i)),'.', 'MarkerSize',15,'Color',COL(i,:))
        plot(1:50, b(i,1)+b(i,2)*(1:50),'Color',COL(i,:), 'LineWidth', 2)
    end
    legend({'10','','15','','20','','30','','40','','65','','115','','165','','215','','265',''},'AutoUpdate', 'off','FontName', 'Times New Roman')
    
    subplot(1,2,1);
    set(gca,'FontName','Times New Roman','FontSize',14)
    xlabel('Days predicted ahead')
    ylabel('Mean proportional prediction error (%)')
    xlim([0 51])
    
    subplot(1,2,2);
    set(gca,'FontName','Times New Roman','FontSize',14)
    xlabel('Days predicted ahead')
    ylabel('Mean absolute prediction error (kg)')
    xlim([0 51])
    
    saveas(h, ['C:\Users\u0084712\Box Sync\Documents\IWT-LA mastitis\Research\S1_Reference_modelling\Results\F' num2str(FarmN) '_L' num2str(Level) '_Result_CV.fig'])
end




