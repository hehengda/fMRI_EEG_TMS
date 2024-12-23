clear all; clc; 
close all;    

addpath('functions')

SUBJECTS={'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8','sub9','sub10',...
'sub11','sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19',...
'sub20','sub21','sub22','sub23','sub24','sub25','sub26','sub27','sub28'}; % 28 subs

response_pre = zeros(400,length(SUBJECTS));

%% ROI labels
atlas_lookup = 'data/Schaefer2018_400Parcels_17Networks_labels_wRGB.txt';
fileID1 = fopen(atlas_lookup,'r') ;
formatSpec = '%f%s%f%f%f%f' ;
atlas_lookup_list = textscan(fileID1 , formatSpec);
fclose(fileID1) ;
roi_names = atlas_lookup_list{2};
roi_RGB(:,1) = atlas_lookup_list{3};
roi_RGB(:,2) = atlas_lookup_list{4};
roi_RGB(:,3) = atlas_lookup_list{5};

roi_net_name = cell(length(roi_names),1);
for j = 1:length(roi_names)
    roi_name_split = strsplit(roi_names{j},'_');
    roi_net_name{j} = roi_name_split{2};
end

net_labels = zeros(400,1);
net_unique = unique(roi_net_name);
for i = 1:length(net_unique)
    net_idx_cell = strfind(roi_net_name,net_unique{i});
    net_idx = find(not(cellfun('isempty',net_idx_cell)));

    net_labels(net_idx) = i;
end

roi_net_name_LHRH = cell(length(roi_names),1);
for j = 1:length(roi_names)
    roi_name_split = strsplit(roi_names{j},'_');
    roi_net_name_LHRH{j} = [roi_name_split{1} '-' roi_name_split{2}];
end

net_labels_LHRH = zeros(400,1);
net_unique_LHRH = unique(roi_net_name_LHRH);
for i = 1:length(net_unique_LHRH)
    net_idx_cell = strfind(roi_net_name_LHRH,net_unique_LHRH{i});
    net_idx = find(not(cellfun('isempty',net_idx_cell)));
    net_labels_LHRH(net_idx) = i;
end

%%
for sub_i = 1:length(SUBJECTS)
    sub = SUBJECTS{sub_i};
    disp(sub)
    
    %% pre treatment sesssion 1
    sub_dir = ['Subjects/' sub '/ses-fmri01'];

    roi_labels = zeros(400,1);
    load([sub_dir,'/ROIs_timeseries_SPM.mat']);

    TMS_regressor_FLOBS = textread([sub_dir '/firstlevel/FLOBS.feat/design.txt']);
    TMS_regressor_FLOBS = TMS_regressor_FLOBS(:,1:3)';

    count = 1;
    % LH
    for label = 1001:1200
        roi_labels(count) = label;
        b = glmfit(TMS_regressor_FLOBS',roi_ts(count,:));
        response_pre(count,sub_i) = sign(b(2)).*sqrt(sum(b(2:4).^2));
        
        count = count + 1;
    end

    % RH
    for label = 2001:2200
        roi_labels(count) = label;
        b = glmfit(TMS_regressor_FLOBS',roi_ts(count,:));
        response_pre(count,sub_i) = sign(b(2)).*sqrt(sum(b(2:4).^2));
        
        count = count + 1;
    end
end

%%
load data/hengda_subinfo.mat

sub_HAMD_post_pre = zeros(1,length(SUBJECTS));
sub_HAMD_pre = zeros(1,length(SUBJECTS));
sync = zeros(1,length(SUBJECTS));
idlist = cell(1,length(SUBJECTS));

for sub = 1:length(SUBJECTS)
    id = SUBJECTS{sub};
    
    Index = find(ismember(sub_info.NDAR,id(4:end)));
    
    if isempty(Index)
        sub_HAMD_post_pre(sub) = NaN;
        sub_HAMD_pre(sub) = NaN;
        sync(sub) = NaN;
        idlist{sub} = NaN;
    else
        sub_HAMD_post_pre(sub) = (sub_info.HAMD(1,Index) - sub_info.HAMD(8,Index))./sub_info.HAMD(1,Index);
        sub_HAMD_pre(sub) = sub_info.HAMD(1,Index);
        findunsync = strfind(sub_info.sync(Index),'unsync');
        sync(sub) = cellfun('isempty',findunsync);
        idlist{sub} = id(end-2:end);
    end

end

%%
numofNETs = 34;
response_pre_net = zeros(numofNETs,length(SUBJECTS));
for sub_i = 1:length(SUBJECTS)
    for i = 1:numofNETs
        select_i = find(net_labels_LHRH == i);
        response_pre_net(i,sub_i) = sum(response_pre(select_i,sub_i));
    end
end

%% spider PLOt response

net_num = 34;

reorder = [linspace(net_num,net_num/2+1,net_num/2),linspace(1,net_num/2,net_num/2)];
input = response_pre_net(reorder,:)';

mid_plot = median(input,1);
low_plot = quantile(input,0.25);
high_plot = quantile(input,0.75);

plot_input = [low_plot;mid_plot;high_plot] ;
plot_label = net_unique_LHRH(reorder);

lim_low = (1-0.1*sign(min(plot_input(:))))*min(plot_input(:));
lim_high = 1.1*max(plot_input(:));
interval = 6;
unit_len = (lim_high-lim_low)/interval;

a = plot_input(:,1)+unit_len;
b = plot_input(:,end)+unit_len;
mid_align = 2*a.*b./(a+b)*cos(2*pi/(net_num+1))-unit_len;
plot_input = [mid_align plot_input];
plot_label = [{''}; plot_label];

figure,spider_plot_HH(plot_input,'AxesOffset', 1,...
    'AxesLimits',[lim_low*ones(1,net_num+1);lim_high*ones(1,net_num+1)],...
    'AxesLabels', plot_label,...
    'AxesDisplay', 'one',...
    'AxesInterval', interval,...
    'AxesLabelsOffset', 0.05,...
    'AxesPrecision',2,...
    'marker','.',...
    'axesfontsize',22,...
    'labelfontsize',22,...
    'axescolor',[0.7 0.7 0.7],...
    'linetransparency',0.8,...
    'color',[0.4 0.4 0.9;0.9 0.4 0.4;0.4 0.4 0.9],...
    'AxesLabelsRotate', 'on',...
    'axeslabelsedge',[1 1 1]);

%% before treatment (<->) and HAMD improvement

indx = ~(isnan(sub_HAMD_post_pre))&(sync==1);
R_all = zeros(34,1);
P_all = zeros(34,1);
for i = 1:34
    [R,P] = corrcoef(sub_HAMD_post_pre(indx),response_pre_net(i,indx));
    R_all(i) = R(1,2);
    P_all(i) = P(1,2);
end

net_unique_LHRH(find(P_all<(0.05/34)))
R_all(find(P_all<(0.05/34)))
P_all(find(P_all<(0.05/34)))

%% before treatment (<->) and HAMD improvement

indx = ~(isnan(sub_HAMD_post_pre))&(sync==0);
R_all = zeros(34,1);
P_all = zeros(34,1);
for i = 1:34
    [R,P] = corrcoef(sub_HAMD_post_pre(indx),response_pre_net(i,indx));
    R_all(i) = R(1,2);
    P_all(i) = P(1,2);
end

net_unique_LHRH(find(P_all<(0.05/34)))
R_all(find(P_all<(0.05/34)))
P_all(find(P_all<(0.05/34)))

%% LH-ContB
indx = ~(isnan(sub_HAMD_post_pre));
xbeta = sub_HAMD_post_pre(indx&(sync==1));
ybeta = response_pre_net(2,indx&(sync==1));
figure,
scatter(xbeta,ybeta,'or','LineWidth',2)
ylabel(['Evoked response at pre-treatment scan'])
xlabel(['HAMD percent improvement'])
title([net_unique_LHRH(2) ' sync(red) and unsync(blue)' ])
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--r') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['sync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.55, 0.75, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

hold on
xbeta = sub_HAMD_post_pre(indx&(sync==0));
ybeta = response_pre_net(2,indx&(sync==0));
scatter(xbeta,ybeta,'ob','LineWidth',2)
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--b') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['unsync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.65, 0.3, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

%% RH-ContB
indx = ~(isnan(sub_HAMD_post_pre));
xbeta = sub_HAMD_post_pre(indx&(sync==1));
ybeta = response_pre_net(19,indx&(sync==1));
figure,
scatter(xbeta,ybeta,'or','LineWidth',2)
ylabel(['Evoked response at pre-treatment scan'])
xlabel(['HAMD percent improvement'])
title([net_unique_LHRH(19) ' sync(red) and unsync(blue)' ])
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--r') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['sync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.65, 0.75, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

hold on
xbeta = sub_HAMD_post_pre(indx&(sync==0));
ybeta = response_pre_net(19,indx&(sync==0));
scatter(xbeta,ybeta,'ob','LineWidth',2)
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--b') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['unsync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.65, 0.18, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

%% RH-LimbicB
indx = ~(isnan(sub_HAMD_post_pre));
xbeta = sub_HAMD_post_pre(indx&(sync==1));
ybeta = response_pre_net(27,indx&(sync==1));
figure,
scatter(xbeta,ybeta,'or','LineWidth',2)
ylabel(['Evoked response at pre-treatment scan'])
xlabel(['HAMD percent improvement'])
title([net_unique_LHRH(27) ' sync(red) and unsync(blue)' ])
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--r') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['sync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.25, 0.55, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

hold on
xbeta = sub_HAMD_post_pre(indx&(sync==0));
ybeta = response_pre_net(27,indx&(sync==0));
scatter(xbeta,ybeta,'ob','LineWidth',2)
set(gca,'FontSize',20)
grid on

hold on
mdl = fitlm(xbeta,ybeta);
intercept = table2array(mdl.Coefficients(1,1));
beta = table2array(mdl.Coefficients(2,1));
pValue = table2array(mdl.Coefficients(2,4));
tvalue = table2array(mdl.Coefficients(2,3));

xi = linspace(min(xbeta),max(xbeta)) ;
yi = beta*xi+intercept;
plot(xi,yi,'--b') ;
[R,P] = corrcoef(xbeta,ybeta);

str = {['unsync'],['\beta = ',num2str(beta,'%.4d')],['p = ',num2str(pValue,'%.4d')],['r = ',num2str(R(1,2),'%.4d')]};
a = annotation('textbox', [0.5, 0.3, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

