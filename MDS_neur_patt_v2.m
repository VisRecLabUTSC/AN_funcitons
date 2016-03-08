
function MDS_neur_patt_v2

%%%load bhv data - constr by constr_mn_patt_vfMRI.m
fold_res='../../bhv_exp/anl_results/grp_res/';
conf_mat_vect_bhv=dlmread([fold_res, 'conf_sym_exp_vect_E2_vfMRI_bysubj.txt']);
%conf_mat_vect_bhv=dlmread([fold_res, 'conf_sym_exp_vect_E2_vfMRI_CP1.txt']);
%%%selext normal mean or CP1, 2
subj_sel=1:8;%1:8;%9, 10
conf_mat_vect_bhv=mean(conf_mat_vect_bhv(:, subj_sel), 2); %conf_mat_vect_bhv(:, 9) or 10

%%%load neur data - constr by constr_svm_corr_ROI_neurconsist.m
fold_res='../anl_mvpa/grp_res/';
load([fold_res, 'ROI_perf_mat_conf_alt_8normals.mat'], 'ROI_perf_mat_conf') %1770 x subj_n x ROI_n

%%%constr var: 1 - cross-set pairs; 0 - within-set pairs
%%%to be used for config (not typ where set eff' are averaged out)
stim_subsets=dlmread(['../anl_mvpa/data_det/id_runassign.txt']);
idpair_list=nchoosek(1:60,2);
set_ind=zeros(1770,1);
for subset_k=1:2
    tmp=ismember(idpair_list, stim_subsets(:,subset_k));
    set_ind=set_ind+single(tmp(:,1)&tmp(:,2))*subset_k;
end
set_ind=1-single(set_ind>0);

%%%rFG, lFG, lIFG
ROI_n=size(ROI_perf_mat_conf, 3);


fold_sols='MDS_sols/';
[jnk1, jnk2]=mkdir(fold_sols);

for ROI_k=0:0 %1:11
    %ROI_n %0 - bhv data
    
     if ROI_k==0
         conf_mat_vect=conf_mat_vect_bhv;
         
     else conf_mat_vect=mean(ROI_perf_mat_conf(:,subj_sel, ROI_k), 2);
         
         %%%regress out set eff from fMRI data
         [b,bint,resid] = regress(conf_mat_vect,set_ind);
         conf_mat_vect=resid;
     end
     
     
     
     conf_mat_sym=squareform(conf_mat_vect);
     
     %%%%%%%%%% compute MDS
 
     %%%make it pos and preserve 0 diag (to run MDS) - needed for neur data (neg d')
     conf_mat_sym=conf_mat_sym - min(conf_mat_sym(:));
     conf_mat_sym=conf_mat_sym.*(1-diag(diag(ones(60))));
    

     [Y,eigs] = cmdscale(conf_mat_sym);
     dims_pos = sum(eigs >0.0001); %nmb of pos dims
     Y=Y(:, 1:dims_pos);
%      sum(Y(:,1:10))
%      std(Y(:,1:10))
     
     
     %%%%check on perc explained var
     %eigs=eigs(1:dims_pos,1);
     %perc_expl=eigs/sum(eigs);
     %perc_expl_cum=cumsum(eigs)/sum(eigs);
     %perc_exp_summ_confMDS=[perc_expl perc_expl_cum]
     
     dim_max=20; %restrict #dims for recon purposes
     Y_mat=NaN(60, dim_max, 60);%60 L1out plus the original (all-in) MDS
     
     for ind_k=1:60
         
         %ind_k=ind_k
         ind_train=setdiff(1:60, ind_k);
         conf_mat_sym_L1out=conf_mat_sym(ind_train, ind_train);
         
         [Y_L1out,eigs_L1out] = cmdscale(conf_mat_sym_L1out);
         
         dims_pos_L1out = sum(eigs_L1out >0.0001); %nmb of pos dims
         dims_pos_L1out=min(dims_pos_L1out, dims_pos);

         Y_L1out=Y_L1out(:, 1:dims_pos_L1out);
         %eigs_L1out=eigs_L1out(1:dims_pos,1);
         
         [fit_val,Y_proj_train,transform] = procrustes(Y_L1out,Y(ind_train, 1:size(Y_L1out, 2)));
         c = transform.c;
         T = transform.T;
         b = transform.b;
         Y_proj = b*Y(:, 1:size(Y_L1out, 2))*T + repmat(c(1,:), 60, 1);
         
%compare orig coord vs projected one side-by-side         
%          Y_view(:,1:2:10)=Y(:,1:5);
%          Y_view(:,2:2:10)=Y_proj(:,1:5);
%          Y_view

        Y_mat(ind_train, :, ind_k)=Y_L1out(:, 1:dim_max);
        Y_mat(ind_k, :, ind_k)=Y_proj(ind_k, 1:dim_max);
         
         
     end
     
     MDS_sol_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_L1out.mat'];
     save(MDS_sol_fl, 'Y_mat')
     
     MDS_orig_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_orig.txt'];
     dlmwrite(MDS_orig_fl, Y)
     
end

         
    