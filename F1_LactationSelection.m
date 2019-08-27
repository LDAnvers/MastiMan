function [dataset, Summary, Delsum] = F1_LactationSelection(InData, OutData, FarmN)
% This script aims to construct a dataset ready for training a mixed model,
% and produces a general overview of the retained and discarded lactations
%
% INPUTS:   InData      - dataset provided as a table, containing 
%                               BA, Lac, LFMY/RFMY/LHMY/RHMY, DIM, MI, EC,
%                                           if Outdata = 1
%                               BA, Lac, TMY, DIM, MI, LFEC/RFEC/LHEC/RHEC, Dest 
%                                           if OutData = 2
%           OutData     - whether output should be QMY (1) or TMY (2) data
%           Nfarm       - the farm number to be added to the output dataset
%
% OUTPUTS:  dataset     - a new dataset containing:
%                               FarmN  Farm number 
%                               BA     Basic Animal Identifyer
%                               Lac    Lactation number
%                               DIM    Days in milk of the milking
%                               MI     Milking interval
%                               LP     Lactation parity NOMINAL
%                               Q  (1) Quarter number
%                               QP (1) Quarter place NOMINAL      
%                               Cat (1)Category LP / QP NOMINAL
%                               T      Teller NOMINAL (for rand)
%                               TMY (2)  total milk yield
%                               QMY (1)  yield first quarter
%                               logTMY log total milk yield
%                               logQMY log quarter milk yield
%                               LogDIM Linearized days in milk (log)
%            Summary      - summary of the  data + why deleted
%
% STEP 0: Make sure all lactations are properly separated by plotting all
%         data per lactation - in separate script!
%
% STEP 1: Check requirements of the input arguments to start the script
%      STEP 1a: check validity of OutData
%      STEP 1b: check validity of InData - table format
%      STEP 1c: check validity of InData - variable names and completeness
%
% STEP 2: Selection procedure
%      STEP 2a: selection on start, length and completeness of the lactation
%      STEP 2b: selection on milk yield (QMY / TMY)
%
% STEP 3: Summary and output definition


%% STEP 1: check requirements for the OutData and InData
% STEP 1a: check validity of OutData
if OutData < 1 || OutData >= 3                 % nonvalid input argument 'OutData'
    warning('Please define valid OutData')
    return
end

% STEP 1b: check validity of InData - table format
if istable(InData) == 0
    warning('InData should be a table')
    return
end

% STEP 1c: check validity of InData - Variable names and completeness
VN = InData.Properties.VariableNames;
C = 0;                                          % counter for N° lacking vars
if sum(ismember(VN,'BA')) < 1;      warning('No BA variable');   C=C+1;   end
if sum(ismember(VN,'Lac')) < 1;     warning('No Lac variable');  C=C+1;   end
if sum(ismember(VN,'DIM')) < 1;     warning('No DIM variable');  C=C+1;   end
if sum(ismember(VN,'MI')) < 1;      warning('No MI variable');   C=C+1;   end
% if sum(ismember(VN,'ECLF')) < 1;    warning('No ECLF variable'); C=C+1;   end
% if sum(ismember(VN,'ECRF')) < 1;    warning('No ECRF variable'); C=C+1;   end
% if sum(ismember(VN,'ECLH')) < 1;    warning('No ECLH variable'); C=C+1;   end
% if sum(ismember(VN,'ECRH')) < 1;    warning('No ECRH variable'); C=C+1;   end
if sum(ismember(VN,'TMY')) < 1;     warning('No TMY variable');  C=C+1;   end

if OutData == 1 || isempty(OutData) == 1        % Defined QMY or default = QMY
    if sum(ismember(VN,'MYLF')) < 1;     warning('No MYLF variable'); C=C+1;   end
    if sum(ismember(VN,'MYRF')) < 1;     warning('No MYLF variable'); C=C+1;   end
    if sum(ismember(VN,'MYLH')) < 1;     warning('No MYLF variable'); C=C+1;   end
    if sum(ismember(VN,'MYRH')) < 1;     warning('No MYLF variable'); C=C+1;   end
end

if C > 0            % quit function if one of the variables are missing
    return
end
clear C

%% STEP 2: Selection procedure
% STEP 2a: selection on start, length and completeness of the lactations
%          this step does not take in QMY, only TMY yields

% Summarize InData
Summary = array2table(unique([InData.BA(:,1) InData.Lac(:,1)],'rows'),'VariableNames',{'BA','Lac'});
Summary.Start(:,1) = NaN;                           % fill in startDIM
Summary.End(:,1) = NaN;                             % fill in endDIM
Summary.Nmeas(:,1) = NaN;                           % fill in N° measurements

Summary(isnan(Summary.Lac) == 1,:)=[];              % delete rows without LacN

for i = 1:length(Summary.BA)                        % scan and summarize data
    ind = find(InData.BA == Summary.BA(i) & InData.Lac == Summary.Lac(i));
    if isempty(ind) == 0
        Summary.Start(i,1) = min(InData.DIM(ind));  % fill in startDIM
        Summary.End(i,1) = max(InData.DIM(ind));    % fill in endDIM
        Summary.Nmeas(i,1) = length(ind);           % fill in n° measurements
    end
end
clear i ind

% Selection 1: based on summary / criterion start < 4 DIM & end > 200 DIM
Summary.Del1(:,1) = 0;                              % fill in deletion 1
Summary.Del1(Summary.Start > 4 | Summary.End < 200,1) = 1; % these are not eligible

% STEP 2b: selection on milk yield (QMY / TMY) - only take into account
%          first 305 days (to fit models on) and milkings with MI > 5h
Wood = @(p,t) p(1).*t.^(p(2)).*exp(-p(3).*t);       % Wood model
opts = optimset('Display','off');
p0 = [0.2 0.2 0.005];                               % initial pars Wood
if OutData == 1 || isempty(OutData) == 1            % Defined QMY or default = QMY
    Summary.Del2a(:,1) = 0;                         % fill in selection 2a
     Summary.Del2b(:,1) = 0;                         % fill in selection 2b
      Summary.Del2c(:,1) = 0;                         % fill in selection 2c
       Summary.Del2d(:,1) = 0;                         % fill in selection 2d
    ColIndex(1) = find(strcmp(VN,'MYLF'));          % index 1 of MYLF
     ColIndex(2) = find(strcmp(VN,'MYRF'));          % index 2 of MYRF
      ColIndex(3) = find(strcmp(VN,'MYLH'));          % index 3 of MYLH
       ColIndex(4) = find(strcmp(VN,'MYRH'));          % index 4 of MYRH
    
    for i = 1:length(Summary.BA(Summary.Del1==0))   
        subset = [Summary.BA(Summary.Del1==0)  Summary.Lac(Summary.Del1==0)]; % non deleted in first round
        ind = find(InData.BA == subset(i,1) & InData.Lac == subset(i,2) & InData.DIM < 305 & InData.MI > 5 & isnan(InData.MI)==0 & isnan(InData.TMY)==0);    % find data cow & lac i
        k = find(Summary.BA == subset(i,1) & Summary.Lac == subset(i,2));   % find index of subset(i) in Summary

        for j = 1:4
            Qdata = [InData.DIM(ind,1),InData{ind,ColIndex(j)}./InData.MI(ind)];  % select Qdata
            
            %====== POSSIBILITY 1 = 2-step approach via Wood model ======%
                p1 = lsqcurvefit(Wood,p0,Qdata(:,1),Qdata(:,2),[],[],opts);      % plot Wood curve
                idxH = find(Qdata(:,2) > 0.8*(Wood(p1,Qdata(:,1))));     % all indices of milkings above 80% 
                
                p2 = lsqcurvefit(Wood,p0,Qdata(idxH,1),Qdata(idxH,2),[],[],opts);   % refit using Wood for only milkings above 80%
                crit = (Qdata(:,2) < 0.8*(Wood(p2,Qdata(:,1))));         % find all lower than 80% of new Wood
                result = (strfind((crit(10:end-10)'),[1 1 1 1 1 1 1 1 1 1 1]));  % Find 10 successive milkings below 80%
                if isempty(result)==0
                    Summary{k,6+j} = 1;             % delete these milking
                else
                    Summary{k,6+j} = 0; 
                end     
                if sum(Qdata(Qdata(:,2)<0.05)) > 0.3*length(Qdata(:,1)); Summary{k,6+j} = 1; end % delete if >30% < 0.05kg
%                 figure(i); subplot(2,2,j);hold on; xlabel('DIM (days)');ylabel('MY (kg/h)'); ylim([0 1.1])
%                 plot(Qdata(:,1),Qdata(:,2),'k.','MarkerSize',12)
%                 plot(Qdata(:,1),Wood(p1,Qdata(:,1)),'r','LineWidth',1.5)
%                 idxL = find(Qdata(:,2) < 0.8*(Wood(p1,Qdata(:,1))));     % all indices of milkings below 80%
%                 plot(Qdata(idxL,1), Qdata(idxL,2),'rx','LineWidth',2,'MarkerSize',7)
%                 plot(Qdata(:,1),Wood(p2,Qdata(:,1)),'b','LineWidth',1.5)

            %====== POSSIBILITY 2 = 80% approach via stepwise smoothing ======%
                % Here should be the MoSAR APPROACH BE IMPLEMENTED
                % This is based on an iterative smoothing process
                % in which the milk lactation is smoothed until 80% increase
                % To be added when discussed with NIC/AHMED/PIERRE
        end
    end
        
elseif OutData == 2                                 % Defined TMY table = output
    Summary.Del2(:,1) = 0;                          % fill in selection 2
    subset = [Summary.BA(Summary.Del1==0)  Summary.Lac(Summary.Del1==0)]; % non deleted in first round
    
    for i = 1:length(subset(:,1))
        ind = find(InData.BA == subset(i,1) & InData.Lac == subset(i,2) & InData.DIM < 305 & InData.MI > 5);    % find data cow & lac i
        k = find(Summary.BA == subset(i,1) & Summary.Lac == subset(i,2));   % find index of subset(i) in Summary
        Data = [InData.DIM(ind) InData.TMY(ind)./InData.MI(ind)];
        
        %====== POSSIBILITY 1 = 2-step approach via Wood model ======%
            p1 = lsqcurvefit(Wood,p0,Data(:,1),Data(:,2),[],[],opts);         % plot Wood curve

            idxH = find(Data(:,2) > 0.8*(Wood(p1,Data(:,1))));     % all indices of milkings above 80%
        
            p2 = lsqcurvefit(Wood,p0,Data(idxH,1),Data(idxH,2),[],[],opts);   % refit using Wood for only milkings above 80%

            crit = (Data(:,2) < 0.8*(Wood(p2,Data(:,1))));         % find all lower than 80% of new Wood
            result = (strfind((crit(10:end-10)'),[1 1 1 1 1 1 1 1 1]));  % Find 9 successive milkings below 80%
            if isempty(result)==0 
                Summary.Del2(k) = 1; % delete these lact (Sum)
            else
                Summary.Del2(k) = 0;
            end     
%               figure(i);hold on; xlabel('DIM (days)');ylabel('MY (kg/h)');
%               plot(Data(:,1),Data(:,2),'ko','MarkerSize',3, 'LineWidth',2,'MarkerFaceColor','k')
%               plot(Data(:,1),Wood(p1,Data(:,1)),'r','LineWidth',1.5)
%               idxL = find(Data(:,2) < 0.8*(Wood(p1,Data(:,1))));     % all indices of milkings below 80%
%               plot(Data(idxL,1), Data(idxL,2),'rx','LineWidth',2,'MarkerSize',7)
%               plot(Data(:,1),Wood(p2,Data(:,1)),'b','LineWidth',1.5)
%               idxL2 = find(Data(:,2) < 0.8*(Wood(p2,Data(:,1))));     % all indices of milkings below 80%
%               plot(Data(idxL2,1), Data(idxL2,2),'bo','LineWidth',2,'MarkerSize',7)
        %====== POSSIBILITY 2 = 80% approach via stepwise smoothing ======%
            % Here should be the MoSAR APPROACH BE IMPLEMENTED
            % This is based on an iterative smoothing process
            % in which the milk lactation is smoothed until 80% increase
            % To be added when discussed with NIC/AHMED/PIERRE
        
    end
end
clear ind idxL idxH i j k Qdata Data crit ColIndex p0 p1 p2 result subset VN Wood ind 
% STEP 2c - Selection on EC or Sep  ( I think it unnecessary)
% to be implemented later if needed


%% STEP 3: Summary and output definition
VN = InData.Properties.VariableNames;               % get variable names

if OutData == 1 || isempty(OutData) == 1            % Defined QMY or default = QMY
    dataset = array2table(zeros(1000000,15),'VariableNames',{'FarmN','BA','Lac','DIM','QMY','MI','pMI','logMI','T','LP','QP','Q','Cat','logDIM','logMY'});
    T = 1;                      % subject counter
    TM = 0;                     % milkings counter
    ColIndex(1) = find(strcmp(VN,'MYLF'));          % index 1 of MYLF
     ColIndex(2) = find(strcmp(VN,'MYRF'));         % index 2 of MYRF
      ColIndex(3) = find(strcmp(VN,'MYLH'));        % index 3 of MYLH
       ColIndex(4) = find(strcmp(VN,'MYRH'));       % index 4 of MYRH
    for i = 1:length(Summary.BA)
        for j = 1:4
            if sum([Summary{i,6} Summary{i,6+j}]) == 0  % only proceed when Qlac is eligible
                ind = find(InData.BA == Summary.BA(i) & InData.Lac == Summary.Lac(i));
                L = length(ind);                    % N° milkings
                dataset.FarmN(TM+1:TM+L) = FarmN;   % fill in Farm number
                dataset.BA(TM+1:TM+L) = Summary.BA(i);      % fill in BA
                dataset.Lac(TM+1:TM+L) = Summary.Lac(i);    % fill in Lac
                dataset.DIM(TM+1:TM+L) = InData.DIM(ind);   % fill in DIM
                dataset.QMY(TM+1:TM+L) = InData{ind,ColIndex(j)};   % fill in QMY
                dataset.MI(TM+1:TM+L) = InData.MI(ind);     % fill in MI
                dataset.T(TM+1:TM+L) = T; T=T+1;            % fill in Tel
                dataset.Q(TM+1:TM+L) = j;                   % fill in quarter
                dataset.pMI(TM+1) = NaN;            % fill in first previous milking interval
                dataset.pMI(TM+2:TM+L) = dataset.MI(TM+1:TM+L-1);     % fill in previous milking interval
                dataset.logMI(TM+1:TM+L) = log(InData.MI(ind));       % fill in MI
                if Summary.Lac(i) == 1 && j < 3     % category 1 front/first
                    dataset.LP(TM+1:TM+L) = 1;      % first lactation
                    dataset.QP(TM+1:TM+L) = 1;      % front quarter
                    dataset.Cat(TM+1:TM+L) = 1;     % Category 1
                elseif Summary.Lac(i) == 1 && j > 2 % category 2 hind/first
                    dataset.LP(TM+1:TM+L) = 1;      % first lactation
                    dataset.QP(TM+1:TM+L) = 2;      % hind quarter
                    dataset.Cat(TM+1:TM+L) = 2;     % Category 2
                elseif Summary.Lac(i) > 1 && j < 3  % category 3 front/higher
                    dataset.LP(TM+1:TM+L) = 2;      % higher lactation
                    dataset.QP(TM+1:TM+L) = 1;      % front quarter
                    dataset.Cat(TM+1:TM+L) = 3;     % Category
                else                                % category 4 hind/higher
                    dataset.LP(TM+1:TM+L) = 2;      % higher lactation
                    dataset.QP(TM+1:TM+L) = 2;      % hind quarter
                    dataset.Cat(TM+1:TM+L) = 4;     % Category
                end
                dataset.logDIM(TM+1:TM+L) = log(InData.DIM(ind)+1);   % fill in log(DIM + 1)
                dataset.logMY(TM+1:TM+L) = log(dataset.QMY(TM+1:TM+L)*1000);   % fill in log(MY*1000)
                TM = TM+L;
            end
        end
    end

    dataset.Q = categorical(dataset.Q);             % set nominal
    dataset.QP = categorical(dataset.QP);           % set nominal
    dataset.LP = categorical(dataset.LP);           % set nominal
    dataset.T = categorical(dataset.T);             % set nominal
    dataset(dataset.BA == 0,:) = [];                % delete empty rows
    
elseif OutData == 2                                 % Defined TMY table = output
    dataset = array2table(zeros(1000000,12),'VariableNames',{'FarmN','BA','Lac','DIM','TMY','MI','pMI','logMI','T','LP','logDIM','logMY'});
    T = 1;                      % subject counter
    TM = 0;                     % milkings counter
    for i = 1:length(Summary.BA)
        if sum([Summary.Del1(i) Summary.Del2(i)]) == 0  % only proceed when Lac is eligible
                ind = find(InData.BA == Summary.BA(i) & InData.Lac == Summary.Lac(i));
                L = length(ind);                    % N° milkings
                dataset.FarmN(TM+1:TM+L) = FarmN;   % fill in Farm number
                dataset.BA(TM+1:TM+L) = Summary.BA(i);      % fill in BA
                dataset.Lac(TM+1:TM+L) = Summary.Lac(i);    % fill in Lac
                dataset.DIM(TM+1:TM+L) = InData.DIM(ind);   % fill in DIM
                dataset.TMY(TM+1:TM+L) = InData.TMY(ind);   % fill in QMY
                dataset.MI(TM+1:TM+L) = InData.MI(ind);     % fill in MI
                dataset.T(TM+1:TM+L) = T; T=T+1;            % fill in Tel
                dataset.logMI(TM+1:TM+L) = log(InData.MI(ind));       % fill in MI
                if Summary.Lac(i) == 1              % category 1 front/first
                    dataset.LP(TM+1:TM+L) = 1;      % first lactation
                else
                    dataset.LP(TM+1:TM+L) = 2;      % second lactation
                end
                dataset.logDIM(TM+1:TM+L) = log(InData.DIM(ind)+1);   % fill in log(DIM + 1)
                dataset.logMY(TM+1:TM+L) = log(dataset.TMY(TM+1:TM+L)*1000);   % fill in log(MY*1000)
                dataset.pMI(TM+1) = NaN;            % fill in first previous milking interval
                dataset.pMI(TM+2:TM+L) = dataset.MI(TM+1:TM+L-1);     % fill in previous milking interval
                TM = TM+L;
        end
    end
    dataset.LP = categorical(dataset.LP);               % set nominal
    dataset.T = categorical(dataset.T);                 % set nominal
    dataset(dataset.BA == 0,:) = [];                % delete empty rows    
end

% SUMMARIZE DATA & DELETION
if OutData == 1 || isempty(OutData) == 1            % Defined QMY or default = QMY
    Summary.NDel1(Summary.Del1==1,1) = 4*Summary.Nmeas(Summary.Del1==1);        % N° of meas deleted incomplete
    Summary.NDel2a(Summary.Del2a==1,1) = Summary.Nmeas(Summary.Del2a==1);       % N° Q1 deleted Wood perturbed
    Summary.NDel2b(Summary.Del2b==1,1) = Summary.Nmeas(Summary.Del2b==1);       % N° Q2 deleted Wood perturbed
    Summary.NDel2c(Summary.Del2c==1,1) = Summary.Nmeas(Summary.Del2c==1);       % N° Q3 deleted Wood perturbed
    Summary.NDel2d(Summary.Del2d==1,1) = Summary.Nmeas(Summary.Del2d==1);       % N° Q4 deleted Wood perturbed
    Summary.NDel2(:,1) = sum(Summary{:,12:15},2);                               % Tot N° deleted

    Delsum = array2table([length(unique(Summary.BA)) length(Summary.BA) sum(Summary.Nmeas)*4 ...   'UniqueCows','UniqueLact','Nmeas'
                      sum(Summary.Del1) sum(Summary.NDel1) sum(Summary.NDel1)/(sum(Summary.Nmeas)*4)*100 ... 'NLacDel1','NmeasDel1','ProcDel1'
                      sum(sum(Summary{:,7:10},2)) sum(sum(Summary{:,7:10},2)>0) sum(Summary.NDel2) sum(Summary.NDel2)/(sum(Summary.Nmeas)*4)*100 ... 'NQLacDel2','NUniqueLacDel2','NmeasDel2','ProcDel2'
                      sum(Summary.NDel2)+sum(Summary.NDel1) (sum(Summary.NDel2)+sum(Summary.NDel1))/(sum(Summary.Nmeas)*4)*100 ... 'TotDelmeas','TotprocDel'
                      length(dataset.BA) length(dataset.BA)/(sum(Summary.Nmeas)*4)*100],...'TotRemMeas','TotprocRem'
                   'VariableNames',{'UniqueCows','UniqueLact','Nmeas','NLacDel1','NmeasDel1','ProcDel1',...
                                    'NLacDel2','NUniqueLacDel2','NmeasDel2','ProcDel2',...
                                    'TotDelmeas','TotprocDel','TotRemMeas','TotProcRem'},...
                   'RowNames',{['QMY_' num2str(FarmN)]});
               
else              
    Summary.NDel1(Summary.Del1==1,1) = Summary.Nmeas(Summary.Del1==1);          % N° deleted incomplete
    Summary.NDel2(Summary.Del2==1,1) = Summary.Nmeas(Summary.Del2==1);          % N° deleted Wood perturbed
    
    Delsum = array2table([length(unique(Summary.BA)) length(Summary.BA) sum(Summary.Nmeas) ... 'UniqueCows','UniqueLact','Nmeas'
                      sum(Summary.Del1) sum(Summary.NDel1) sum(Summary.NDel1)/(sum(Summary.Nmeas))*100 ... 'NLacDel1','NmeasDel1','ProcDel1'
                      sum(Summary.Del2) sum(Summary.Del2) sum(Summary.NDel2) sum(Summary.NDel2)/(sum(Summary.Nmeas))*100 ... 'NQLacDel2','NUniqueLacDel2','NmeasDel2','ProcDel2'
                      sum(Summary.NDel2)+sum(Summary.NDel1) (sum(Summary.NDel2)+sum(Summary.NDel1))/sum(Summary.Nmeas)*100 ... 'TotDelmeas','TotprocDel'
                      length(dataset.BA) length(dataset.BA)/(sum(Summary.Nmeas))*100],... 'TotRemMeas','TotprocRem' 
                   'VariableNames',{'UniqueCows','UniqueLact','Nmeas','NLacDel1','NmeasDel1','ProcDel1',...
                                    'NLacDel2','NUniqueLacDel2','NmeasDel2','ProcDel2',...
                                    'TotDelmeas','TotprocDel','TotRemMeas','TotProcRem'},...
                   'RowNames',{['TMY_' num2str(FarmN)]});
end
