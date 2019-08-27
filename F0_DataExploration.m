function [SUM, DataSel, Summary] = F0_DataExploration(data,type, FarmN, SUM)
% Function for data exploration of the different farm datasets; following
% aspects will be calculated:
%       - Number of unique animals (NCows)
%       - Number of unique lactations (NLac)
%       - Number of unique animals with complete lactations (NCCows)
%       - Number of unique lactations with complete lactations (NCLac)
%       - Average Parity of complete lactations (AParity)
%       - SD of parity of complete lactations (SDParity)
%       - Average fit of Wood model (R2) TMY/h on first 300DIM of CLac
%       - SD of fit of Wood model
%       - Average 
%
%
%
%
%
% INPUTS:   data    = dataset of the farm
%           type    = Delaval or Lely (might have different VarNames)
%
% OUTPUTS:  SUM     = summary of the dataset (see above)
%           h1      = figure handle of the first figure
%           h2      = figure handle of the second figure
%           DataSel = New dataset containing only complete lactations -
%                     added to the existing dataset InData saved as
%                     D1_alldatasets; this contains:
%                       all data of data + MILKING INTERVAL + FARMN + type

SUM.FarmN(FarmN,1) = FarmN;
SUM.Type(FarmN,1) = type;
SUM.Start(FarmN,1) = min(data.EndTime);
SUM.End(FarmN,1) = max(data.EndTime);
SUM.NCows(FarmN,1) = length(unique(data.BA));
SUM.NLac(FarmN,1) = size(unique([data.BA data.Lac],'rows'),1);
 
% summarize data for selection
Summary = array2table(unique([data.BA(:,1) data.Lac(:,1)],'rows'),'VariableNames',{'BA','Lac'});
Summary.FarmN(:,1) = FarmN;
Summary.Start(:,1) = NaN;                           % fill in startDIM
Summary.End(:,1) = NaN;                             % fill in endDIM
for i = 1:length(Summary.BA)                        % scan and summarize data
    ind = find(data.BA == Summary.BA(i) & data.Lac == Summary.Lac(i));
    if isempty(ind) == 0
        Summary.Start(i,1) = min(data.DIM(ind));  % fill in startDIM
        Summary.End(i,1) = max(data.DIM(ind));    % fill in endDIM
    end
end
clear i ind

% select eligble cows
ind = find(Summary.Start<4 & Summary.End >200);

SUM.NCCows(FarmN,1) = length(unique(Summary.BA(ind)));
SUM.NCLac(FarmN,1) = length(ind);
SUM.AParity(FarmN,1) = mean(Summary.Lac(ind));
SUM.SDParity(FarmN,1) = std(Summary.Lac(ind));

sel = ismember([data.BA data.Lac],[Summary.BA(ind,1) Summary.Lac(ind,1)],'rows');
DataSel = data(sel,:);
DataSel.MI(:,1) = (datenum(DataSel.EndTime) - datenum(DataSel.PEndTime))*24;
DataSel.MI(DataSel.MI>100,1) = NaN; DataSel(DataSel.MI ==0,:) = [];

SUM.AMI(FarmN,1) = nanmean(DataSel.MI);
SUM.SDMI(FarmN,1) = nanstd(DataSel.MI);
SUM.AMYh(FarmN,1) = nanmean((DataSel.TMY./DataSel.MI).*24);
SUM.SDMYh(FarmN,1) = nanstd((DataSel.TMY./DataSel.MI).*24);
SUM.AMYd100(FarmN,1) = nanmean((DataSel.TMY(DataSel.DIM<100)./DataSel.MI(DataSel.DIM<100)).*24);
SUM.SDMYd100(FarmN,1) = nanstd((DataSel.TMY(DataSel.DIM<100)./DataSel.MI(DataSel.DIM<100)).*24);


Wood = @(p,t) p(1).*t.^(p(2)).*exp(-p(3).*t);       % Wood model
opts = optimset('Display','off');
p0 = [1 0.2 0.004];                               % initial pars Wood

Summary.MY305(:,1) = NaN;
Summary.R2(:,1) = NaN;
Summary.PeakYieldh(:,1) = NaN;
Summary.PeakDIM(:,1) = NaN;
Summary.RMSE(:,1) = NaN;
Summary.p1(:,1) = NaN;
Summary.p2(:,1) = NaN;
Summary.p3(:,1) = NaN;
tic
for i =  1:length(ind)
    
    idx = find(DataSel.BA == Summary.BA(ind(i)) & DataSel.Lac == Summary.Lac(ind(i)) & DataSel.DIM < 306 & isnan(DataSel.TMY)==0 & isnan(DataSel.MI)==0 & DataSel.MI > 4 & DataSel.MI < 24);
    
    DIM = DataSel.DIM(idx);
    TMY = DataSel.TMY(idx)./DataSel.MI(idx);
    
    Summary.MY305(ind(i)) = nansum(TMY);
    
    p = lsqcurvefit(Wood,p0,DIM,TMY,[],[],opts);      % plot Wood curve
%     figure;
%     plot(DIM, TMY,'ob','MarkerSize',4,'MarkerFaceColor','b')
%     hold on
    WOOD = Wood(p, DIM);
%     plot(DIM,WOOD,'r','LineWidth',2)
%     
    Summary.R2(ind(i)) = 1-(sum((TMY-WOOD).*(TMY-WOOD))/sum((TMY-mean(TMY)).*(TMY-mean(TMY)))); % r2
    Summary.RMSE(ind(i)) = mean(sqrt((TMY-WOOD).*(TMY-WOOD)));
    [Summary.PeakYieldh(ind(i)),m] = max(WOOD);
    Summary.PeakDIM(ind(i)) = DIM(m);
    Summary.p1(ind(i)) = p(1);
    Summary.p2(ind(i)) = p(2);
    Summary.p3(ind(i)) = p(3);
end
toc


SUM.AR2Wood(FarmN,1) = nanmean(Summary.R2);
SUM.SDR2Wood(FarmN,1) = nanstd(Summary.R2);
SUM.ARMSEWood(FarmN,1) = nanmean(Summary.RMSE);
SUM.SDRMSEWood(FarmN,1) = nanstd(Summary.RMSE);
SUM.APeakYield(FarmN,1) = nanmean(Summary.PeakYieldh)*24;
SUM.SDPeakYield(FarmN,1) = nanstd(Summary.PeakYieldh)*24;
SUM.APeakYield1(FarmN,1) = nanmean(Summary.PeakYieldh(Summary.Lac ==1))*24;
SUM.SDPeakYield1(FarmN,1) = nanstd(Summary.PeakYieldh(Summary.Lac ==1))*24;
SUM.APeakYield2(FarmN,1) = nanmean(Summary.PeakYieldh(Summary.Lac ~=1))*24;
SUM.SDPeakYield2(FarmN,1) = nanstd(Summary.PeakYieldh(Summary.Lac ~=1))*24;

DataSel.T(:,1) = 1:length(DataSel.BA);
Farm = array2table([(1:length(DataSel.BA))' ones(length(DataSel.BA),1)*FarmN],'VariableNames',{'T','FarmID'});
DataSel = innerjoin(Farm,DataSel,'Keys','T');

SJ = array2table(unique([DataSel.BA DataSel.Lac],'rows'),'VariableNames',{'BA','Lac'});
SJ.SJ(:,1) = 1:length(SJ.BA);
DataSel = outerjoin(DataSel,SJ,'Keys',{'BA','Lac'},'MergeKeys',1);

Summary = Summary(isnan(Summary.p1)==0,:);








