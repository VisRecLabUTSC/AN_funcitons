
% function MDS_neur_patt_v2
%
% %%%load bhv data - constr by constr_mn_patt_vfMRI.m
% fold_res='../../bhv_exp/anl_results/grp_res/';
% conf_mat_vect_bhv=dlmread([fold_res, 'conf_sym_exp_vect_E2_vfMRI_bysubj.txt']);
% %conf_mat_vect_bhv=dlmread([fold_res, 'conf_sym_exp_vect_E2_vfMRI_CP1 .txt']);
% %%%selext normal mean or CP1, 2
% subj_sel=1:8;%1:8;%9, 10
% conf_mat_vect_bhv=mean(conf_mat_vect_bhv(:, subj_sel), 2); %conf_mat_vect_bhv(:, 9) or 10
%
% %%%load neur data - constr by constr_svm_corr_ROI_neurconsist.m
% fold_res='../anl_mvpa/grp_res/';
% load([fold_res, 'ROI_perf_mat_conf_alt_8normals.mat'], 'ROI_perf_mat_conf') %1770 x subj_n x ROI_n
%
% %%%constr var: 1 - cross-set pairs; 0 - within-set pairs
% %%%to be used for config (not typ where set eff' are averaged out)
% stim_subsets=dlmread(['../anl_mvpa/data_det/id_runassign.txt']);
% idpair_list=nchoosek(1:60,2);
% set_ind=zeros(1770,1);
% for subset_k=1:2
%     tmp=ismember(idpair_list, stim_subsets(:,subset_k));
%     set_ind=set_ind+single(tmp(:,1)&tmp(:,2))*subset_k;
% end
% set_ind=1-single(set_ind>0);
%
% %%%rFG, lFG, lIFG
% ROI_n=size(ROI_perf_mat_conf, 3);
%
%
% fold_sols='MDS_sols/';
% [jnk1, jnk2]=mkdir(fold_sols);
%
% for ROI_k=1:11
%     %ROI_n %0 - bhv data
%
%      if ROI_k==0
%          conf_mat_vect=conf_mat_vect_bhv;
%
%      else conf_mat_vect=mean(ROI_perf_mat_conf(:,subj_sel, ROI_k), 2);
%
%          %%%regress out set eff from fMRI data
%          [b,bint,resid] = regress(conf_mat_vect,set_ind);
%          conf_mat_vect=resid;
%      end

constr_prot_on=1;



%%%%%%%%%% format matrix
%%%A full dissimilarity matrix is real and symmetric, has zeros along
%%%the main diagonal and positive elements elsewhere
%      conf_mat_sym=squareform(conf_mat_vect);
conf_mat_sym=conf_mat_sym - min(conf_mat_sym(:));
conf_mat_sym=conf_mat_sym.*(1-diag(diag(ones(60))));

%%%%%%%%%% compute MDS
[Y,eigs] = cmdscale(conf_mat_sym);
dims_pos = sum(eigs >0.0001); %nmb of pos dims
Y=Y(:, 1:dims_pos);
%      sum(Y(:,1:10))
%      std(Y(:,1:10))


%%%%check on perc explained var
eigs=eigs(1:dims_pos,1);
perc_expl=eigs/sum(eigs);
perc_expl_cum=cumsum(eigs)/sum(eigs);
perc_exp_summ_confMDS=[perc_expl perc_expl_cum]

%      dim_max=20; %restrict #dims for recon purposes

%%%%plot MDS results
plot_on=1;
if plot_on
    
    %      if ROI_k==1
    %          Y(:,2)=-Y(:,2);%%%to free upper top (for fig inset) & get some match w/ bhv
    %      end
    
    figure
    
    xlabel('1st dimension')
    ylabel('2nd dimension')
    
    col_vect=[0.8 0 0];
    
    %%%choose to plot raw or zscored solution
    plot(Y(:,1),Y(:,2),'.', 'MarkerEdgeColor',col_vect,'MarkerFaceColor',col_vect, 'MarkerSize',21)
    %plot(zscore(Y(:,1)), zscore(Y(:,2)),'d', 'MarkerEdgeColor',[.6 .1 .1],'MarkerFaceColor',[.6 .1 .1], 'MarkerSize',5)
    
    box off
    
    %%%control (set) aspect ratio & size of plot in inches
    %%%(so reproducible in the future); and figure size (to enclose entire plot)
    
    set(gca, 'Units', 'inches')
    set(gca, 'Position', [0.5 0.5 10 8])
    set(gca,'PlotBoxAspectRatio', [1.25 1 1])
    
    set(gcf, 'Units', 'inches')
    set(gcf, 'Position', [2 2 11 9])
    
    %%%choose axis limits appropriately
    %if ROI_k==0
    %%%bhv plot
    %          axis([-0.45 0.6 -0.6 0.45])
    %          set(gca,'XTick',-0.45:0.15:0.6)
    %          set(gca,'YTick',-0.6:0.15:0.45)
    %elseif ROI_k==1
    %%%rFG plot
    %          axis([-1 1 -1 1])
    %          set(gca,'XTick',-1:0.25:1)
    %          set(gca,'YTick',-1:0.25:1)
    %end
    
    %%%label points with 1-60
    %      nm_mat=1:60;
    
    gimage(ims)
    
end


%% construct 'prototypes' for the first 2 dimensions
ims_new=[];
for i=1:size(ims,1)
    ims_new=cat(4,ims_new,ims{i});
end
if constr_prot_on
    for dim_k=2:2
        
        ind_pos=find(Y(:,dim_k)>0)
        ind_neg=find(Y(:,dim_k)<0)
        error
        Y_pos_scale=Y(ind_pos,dim_k)/sum(Y(ind_pos,dim_k));
        Y_neg_scale=-Y(ind_neg,dim_k)/sum(-Y(ind_neg,dim_k));
        
        
        im_mat_pos_neut=ims_new(:,:,:,ind_pos+60);
        im_mat_neg_neut=ims_new(:, :, :,ind_neg+60);
        im_mat_pos_hap=ims_new(:,:,:,ind_pos);
        im_mat_neg_hap=ims_new(:,:,:,ind_neg);
        
        Y_pos_scale_mat=permute(Y_pos_scale, [2 3 4 1]);
        Y_pos_scale_mat=repmat(Y_pos_scale_mat, [size(ims_new,1) size(ims_new,2) size(ims_new,3) 1]);
        
        Y_neg_scale_mat=permute(Y_neg_scale, [2 3 4 1]);
        Y_neg_scale_mat=repmat(Y_neg_scale_mat, [size(ims_new,1) size(ims_new,2) size(ims_new,3) 1]);
        
        prot_pos_neut=uint8(sum(double(im_mat_pos_neut).*Y_pos_scale_mat, 4));
        prot_neg_neut=uint8(sum(double(im_mat_neg_neut).*Y_neg_scale_mat, 4));
        prot_pos_hap=uint8(sum(double(im_mat_pos_hap).*Y_pos_scale_mat, 4));
        prot_neg_hap=uint8(sum(double(im_mat_neg_hap).*Y_neg_scale_mat, 4));
        
%         imtool(prot_neg_neut)
%         imtool(prot_pos_neut)
%         imtool(prot_neg_hap)
%         imtool(prot_pos_hap)
        
        imwrite(prot_pos_neut,['prot_pos_neut_' 'dim' num2str(dim_k) '.tif']);
        imwrite(prot_neg_neut,['prot_neg_neut_' 'dim' num2str(dim_k) '.tif']);
        imwrite(prot_pos_hap,['prot_pos_hap_' 'dim' num2str(dim_k) '.tif']);
        imwrite(prot_neg_hap,['prot_neg_hap_' 'dim' num2str(dim_k) '.tif']);
        
%         imwrite(prot_neg_neut, [fold_prots, 'ROI', sprintf('%02.0f', ROI_k), '_dim', sprintf('%02.0f', dim_k), '_neg_neut.tif'])
%         imwrite(prot_pos_neut, [fold_prots, 'ROI', sprintf('%02.0f', ROI_k), '_dim', sprintf('%02.0f', dim_k), '_pos_neut.tif'])
%         imwrite(prot_neg_hap, [fold_prots, 'ROI', sprintf('%02.0f', ROI_k), '_dim', sprintf('%02.0f', dim_k), '_neg_hap.tif'])
%         imwrite(prot_pos_hap, [fold_prots, 'ROI', sprintf('%02.0f', ROI_k), '_dim', sprintf('%02.0f', dim_k), '_pos_hap.tif'])
        
        %%%check still same L* means, a* (as it should); checked!!
        %              im_stack(:,:,:,1)=prot_neg_neut; im_stack(:,:,:,2)=prot_pos_neut; im_stack(:,:,:,3)=prot_neg_hap; im_stack(:,:,:,4)=prot_pos_hap;
        %              Lab_mns=NaN(4, 3);
        %              im_stack_Lab=NaN(size(im_stack));
        %              for im_k=1:4
        %                  im_rgb=double(im_stack(:,:,:,im_k));
        %
        %                  cform_srgb2lab = makecform('srgb2lab');
        %                  im_lab_uint8 = applycform(im_rgb, cform_srgb2lab);
        %                  im_lab=lab2double(im_lab_uint8);
        %
        %                  im_stack_Lab(:,:,:, im_k)=im_lab;
        %
        %                  Lab_mns(im_k,:)=permute(sum(sum(im_lab))/(sz1*sz2), [3 1 2]);
        %
        %
    end
    %              Lab_mns=Lab_mns
    %
    %              diff_L_neut_map=im_stack_Lab(:,:,2, 2)-im_stack_Lab(:,:,2, 1);
    %                                 diff_L_neut_map=(diff_L_neut_map+fliplr(diff_L_neut_map))/2;
    %
    %              mx=max(max(abs(diff_L_neut_map)))
    %              diff_L_neut_map=0.5+(diff_L_neut_map/mx)/2;
    %              imtool(diff_L_neut_map)
    %              %imtool(diff_L_neut_map.*single(diff_L_neut_map>0.5))
    %              %imtool(diff_L_neut_map.*single(diff_L_neut_map<0.5))
    %              imtool(diff_L_neut_map>0.5)
else
end


%% Y_mat=NaN(60, dim_max, 60);%60 L1out plus the original (all-in) MDS
% 
% for ind_k=1:60
%     
%     %ind_k=ind_k
%     ind_train=setdiff(1:60, ind_k);
%     conf_mat_sym_L1out=conf_mat_sym(ind_train, ind_train);
%     
%     [Y_L1out,eigs_L1out] = cmdscale(conf_mat_sym_L1out);
%     
%     dims_pos_L1out = sum(eigs_L1out >0.0001); %nmb of pos dims
%     dims_pos_L1out=min(dims_pos_L1out, dims_pos);
%     
%     Y_L1out=Y_L1out(:, 1:dims_pos_L1out);
%     %eigs_L1out=eigs_L1out(1:dims_pos,1);
%     
%     [fit_val,Y_proj_train,transform] = procrustes(Y_L1out,Y(ind_train, 1:size(Y_L1out, 2)));
%     c = transform.c;
%     T = transform.T;
%     b = transform.b;
%     Y_proj = b*Y(:, 1:size(Y_L1out, 2))*T + repmat(c(1,:), 60, 1);
%     
%     %compare orig coord vs projected one side-by-side
%     %          Y_view(:,1:2:10)=Y(:,1:5);
%     %          Y_view(:,2:2:10)=Y_proj(:,1:5);
%     %          Y_view
%     
%     Y_mat(ind_train, :, ind_k)=Y_L1out(:, 1:dim_max);
%     Y_mat(ind_k, :, ind_k)=Y_proj(ind_k, 1:dim_max);
%     
%     
% end
% 
% MDS_sol_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_L1out.mat'];
% save(MDS_sol_fl, 'Y_mat')
% 
% MDS_orig_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_orig.txt'];
% dlmwrite(MDS_orig_fl, Y)
% 
% end
% 
% 
