clear all; clc; close all;

SUBJECTS={'sub1','sub2','sub3','sub4','sub5','sub6','sub7','sub8','sub9','sub10',...
'sub11','sub12','sub13','sub14','sub15','sub16','sub17','sub18','sub19',...
'sub20','sub21','sub22','sub23','sub24','sub25','sub26','sub27','sub28'}; % 28 subs"

numofNETs = 34;
numofROIs = 400;

beta_mat = zeros(numofROIs,numofROIs,length(SUBJECTS));
beta_mat_avgpair = zeros(numofROIs,numofROIs,length(SUBJECTS));
beta_mat_networks = zeros(numofNETs,numofNETs,length(SUBJECTS));

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
   
%
for sub_i = 1:length(SUBJECTS)
    sub = SUBJECTS{sub_i};

    sub_dir = ['Subjects/' sub '/ses-fmri01'];

    load([sub_dir,'/PPI/PPI_task_functional_connectivity.mat'])
    
    beta_mat(:,:,sub_i) = PPI_FC_mat;
    
    for i = 1:numofROIs
        for j = 1:numofROIs
            beta_mat_avgpair(i,j,sub_i) = (PPI_FC_mat(i,j)+PPI_FC_mat(j,i))/2;
        end
    end
    
    % network
    for i = 1:numofNETs
        for j = 1:numofNETs
            
            if i == j
                select_i = find(net_labels_LHRH == i);
                select_j = find(net_labels_LHRH == j);
                beta_mat_networks(i,j,sub_i) = sum(sum(beta_mat_avgpair(select_i,select_j,sub_i)))/2;
            else
                select_i = find(net_labels_LHRH == i);
                select_j = find(net_labels_LHRH == j);
                beta_mat_networks(i,j,sub_i) = sum(sum(beta_mat_avgpair(select_i,select_j,sub_i))); 
            end
        end
    end

end

%% NET level
p_group_mat_net = zeros(numofNETs,numofNETs);
t_group_mat_net = zeros(numofNETs,numofNETs);

    for i = 1:numofNETs
        for j = 1:numofNETs

            [h,p,ci,stats] = ttest(squeeze(beta_mat_networks(i,j,:)));
            p_group_mat_net(i,j) = p;
            t_group_mat_net(i,j) = stats.tstat;
            
        end
    end

% transform T to Z
dof = length(SUBJECTS) - 1;
z_group_mat_net = norminv(tcdf(t_group_mat_net,dof),0,1);

%% visualization net FDR
    
select = triu(p_group_mat_net);
q = 0.05;
allvox_idx = find(select~=0);
p_values = p_group_mat_net(allvox_idx);
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
z_group_mat_net_fdr = zeros(size(z_group_mat_net));
z_group_mat_net_fdr(selectedvox_idx) = z_group_mat_net(selectedvox_idx);
p_group_mat_net_fdr = zeros(size(p_group_mat_net));
p_group_mat_net_fdr(selectedvox_idx) = p_group_mat_net(selectedvox_idx);

input = z_group_mat_net_fdr;

figure,imagesc(input)

edgediffmax = max(max(input));
edgediff_min = min(min(input));

steps = 5000;
mymap = [[linspace(0.5,1,-(edgediff_min)*steps)]',...
[linspace(0.65,1,-(edgediff_min)*steps)]',...
[linspace(0.6,1,-(edgediff_min)*steps)]'];

colormap(mymap), colorbar 

grid on
set(gca,'fontsize',18)
title('FDR')

xticks([1:numofNETs])
xticklabels(net_unique_LHRH)
xtickangle(45)
yticks([1:numofNETs])
yticklabels(net_unique_LHRH)

