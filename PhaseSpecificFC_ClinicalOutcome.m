clear all; clc; close all;

%% brain parcellation network labels
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

% 17 networks label
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

% 34 networks label with left and right
net_labels_LHRH = zeros(400,1);
net_unique_LHRH = unique(roi_net_name_LHRH);
for i = 1:length(net_unique_LHRH)
    net_idx_cell = strfind(roi_net_name_LHRH,net_unique_LHRH{i});
    net_idx = find(not(cellfun('isempty',net_idx_cell)));
    net_labels_LHRH(net_idx) = i;
end

% large ROI labels
net_nodes_labels = zeros(400,1);

for j = 1:length(roi_names)
    roi_name_split = strsplit(roi_names{j},'_');
    namesize = size(roi_name_split,2);
    roi_nodes_name_LHRH{j} = [];
    for nn = 1:(namesize-1)
        roi_nodes_name_LHRH{j} = [roi_nodes_name_LHRH{j} roi_name_split{nn}];
    end
end


numofnodes = 106;

nodes_labels_LHRH = zeros(400,1);
nodes_unique_LHRH = unique(roi_nodes_name_LHRH);
nodes_netlabel = zeros(length(nodes_unique_LHRH),1);
nodes_netlabel_LHRH = zeros(length(nodes_unique_LHRH),1);
nodes_roi_RGB = zeros(length(nodes_unique_LHRH),3);
nodes_unique_LHRH_netname = cell(1,numofnodes);
for i = 1:length(nodes_unique_LHRH)
    net_idx_cell = strfind(roi_nodes_name_LHRH,nodes_unique_LHRH{i});
    net_idx = find(not(cellfun('isempty',net_idx_cell)));
    nodes_labels_LHRH(net_idx) = i;
    netlabel_LHRH_i = net_labels_LHRH(net_idx);
    netlabel_i = net_labels(net_idx);
    nodes_netlabel(i) = netlabel_i(1);
    nodes_netlabel_LHRH(i) = netlabel_LHRH_i(1);
    nodes_roi_RGB(i,:) = roi_RGB(net_idx(1),:);
    
    nodes_unique_LHRH_netname{i} = roi_net_name_LHRH(net_idx(1));
end

numofROIs = length(roi_names);
numofNETs = 34;

%% subject PPI functional connectivity test

SUBJECTS={'sub1','sub2','sub3','sub5','sub6','sub7','sub8','sub9','sub10',...
    'sub11','sub25','sub14','sub28','sub15'}; % subjects with phase bins

beta_mat_nodes_pre = zeros(numofnodes,numofnodes,length(SUBJECTS),4);
beta_mat_nodes_post = zeros(numofnodes,numofnodes,length(SUBJECTS),4);
 
for sub_i = 1:length(SUBJECTS)
    sub = SUBJECTS{sub_i};

    sub_dir = ['Subjects/' sub '/ses-fmri01'];

    load([sub_dir,'/PPI_alpha_phase/PPI_task_FC_allROIs.mat'])
    
    % symmetrize and summarize into nodes FC
    for i = 1:numofROIs
        for j = 1:numofROIs
            beta_mat_avgpair_pre(i,j,sub_i,:) = (PPI_FC_mat_bins(i,j,:)+PPI_FC_mat_bins(j,i,:))/2;
        end
    end
    for i = 1:numofnodes
        for j = 1:numofnodes
            if i == j
                beta_mat_nodes_pre(i,j,sub_i,:) = 0;
            else
                select_i = find(nodes_labels_LHRH == i);
                select_j = find(nodes_labels_LHRH == j);
                beta_mat_nodes_pre(i,j,sub_i,:) = sum(sum(beta_mat_avgpair_pre(select_i,select_j,sub_i,:))); 
            end
        end
    end
    
    %% post
    
    sub_dir = ['Subjects/' sub '/ses-fmri02'];

    load([sub_dir,'/PPI_alpha_phase/PPI_task_FC_allROIs.mat'])
    
    % symmetrize and summarize into nodes FC

    for i = 1:numofROIs
        for j = 1:numofROIs
            beta_mat_avgpair_post(i,j,sub_i,:) = (PPI_FC_mat_bins(i,j,:)+PPI_FC_mat_bins(j,i,:))/2;
        end
    end
    for i = 1:numofnodes
        for j = 1:numofnodes
            if i == j
                beta_mat_nodes_post(i,j,sub_i,:) = 0;
            else
                select_i = find(nodes_labels_LHRH == i);
                select_j = find(nodes_labels_LHRH == j);
                beta_mat_nodes_post(i,j,sub_i,:) = sum(sum(beta_mat_avgpair_post(select_i,select_j,sub_i,:))); 
            end
        end
    end

end

%%
load data/hengda_subinfo.mat

sub_HAMD_post_pre = zeros(1,length(SUBJECTS));
sub_HAMD_post_pre_value = zeros(1,length(SUBJECTS));
sub_HAMD_pre = zeros(1,length(SUBJECTS));
sync = zeros(1,length(SUBJECTS));

for sub = 1:length(SUBJECTS)
    %disp(['Subject - ',subjects{sub}])
    id = SUBJECTS{sub};
    
    Index = find(ismember(sub_info.NDAR,id(4:end)));
    
    if isempty(Index)
        sub_HAMD_post_pre(sub) = NaN;
        sub_HAMD_post_pre_value(sub) = NaN;
        sub_HAMD_pre(sub) = NaN;
        sync(sub) = NaN;
    else
        sub_HAMD_post_pre(sub) = (sub_info.HAMD(1,Index) - sub_info.HAMD(8,Index))./sub_info.HAMD(1,Index);
        sub_HAMD_post_pre_value(sub) = (sub_info.HAMD(8,Index) - sub_info.HAMD(1,Index));
        sub_HAMD_pre(sub) = sub_info.HAMD(1,Index);
        findunsync = strfind(sub_info.sync(Index),'unsync');
        sync(sub) = cellfun('isempty',findunsync);
    end

end


%%

HLPLLPphase=load('data/HLP_LLP_phase_allPatients.mat');
sub_HLPphase = zeros(1,length(SUBJECTS));
sub_LLPphase = zeros(1,length(SUBJECTS));

for sub = 1:length(SUBJECTS)
    
    id = SUBJECTS{sub};

    Index = find(ismember(HLPLLPphase.SUBJECTS,id));
    
    sub_HLPphase(sub) = HLPLLPphase.subject_HLPerphase_DefaultA_PFCd(Index);
    sub_LLPphase(sub) = HLPLLPphase.subject_LLPerphase_DefaultA_PFCd(Index);

end

phase = sub_HLPphase;

beta_mat_nodes_post_new = zeros(numofnodes,numofnodes,length(SUBJECTS));
beta_mat_nodes_pre_new = zeros(numofnodes,numofnodes,length(SUBJECTS));

for i = 1:length(SUBJECTS)
    beta_mat_nodes_post_new(:,:,i) = beta_mat_nodes_post(:,:,i,phase(i));
    beta_mat_nodes_pre_new(:,:,i) = beta_mat_nodes_pre(:,:,i,phase(i));
end

idx_sub = ~(isnan(sub_HAMD_post_pre))&(sync==1);
p_group_mat_nodes = zeros(numofnodes,numofnodes);
r_group_mat_nodes = zeros(numofnodes,numofnodes);

for i = 1:numofnodes
    for j = 1:numofnodes

        [R,P] = corrcoef(sub_HAMD_post_pre(idx_sub),squeeze(beta_mat_nodes_post_new(i,j,idx_sub)-beta_mat_nodes_pre_new(i,j,idx_sub)));
        p_group_mat_nodes(i,j) = P(1,2);
        r_group_mat_nodes(i,j) = R(1,2);

    end
end


% FDR test
test = p_group_mat_nodes;
select = triu(test,1);
q = 0.05;
allvox_idx = find(select~=0);
p_values = test(allvox_idx);
[p_v_sort idx] = sort(p_values);
V = length(p_values);
r=0;
for i = 1:V
    if p_v_sort(i)<=(i/V)*q
        r = i;
    else
    end
end
selectedvox_idx=allvox_idx(idx(1:r));

[a b] = ind2sub(size(test),selectedvox_idx);
for i = 1:length(a)
    disp([nodes_unique_LHRH{a(i)} ' - ' nodes_unique_LHRH{b(i)} ' r = ' num2str(r_group_mat_nodes(a(i),b(i))) ';p = ' num2str(p_group_mat_nodes(a(i),b(i)))]);
end

%% scatter plot
% LHDefaultBPFCd - RHLimbicBOFC
% 19 - 84

indx = ~(isnan(sub_HAMD_post_pre));
xbeta = sub_HAMD_post_pre(indx&sync);
ybeta = squeeze(beta_mat_nodes_post_new(19,84,indx&sync)-beta_mat_nodes_pre_new(19,84,indx&sync));

figure,
scatter(xbeta,ybeta,'or','LineWidth',2)
ylabel(['Modulated functional connectivity change (post - pre)'])
xlabel(['HAMD percent improvement'])
title([strrep(nodes_unique_LHRH{19},'_','-') ' & ' strrep(nodes_unique_LHRH{84},'_','-') ' sync(red) and unsync(blue)' ])
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
a = annotation('textbox', [0.6, 0.65, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;

hold on
xbeta = sub_HAMD_post_pre(indx&(~sync));
ybeta = squeeze(beta_mat_nodes_post_new(19,84,indx&(~sync))-beta_mat_nodes_pre_new(19,84,indx&(~sync)));
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
a = annotation('textbox', [0.7, 0.25, 0.1, 0.1], 'String', str,'LineStyle','none');
a.FontSize = 18;



