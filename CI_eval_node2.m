
function CI_eval_node2


cond_col=1;
if cond_col
    cond_nm='col';
else cond_nm='gray';
end

flip=0;
if ~flip
    fold_fl=['../../design_stims/stims/exp_categ1', '_', cond_nm, '.mat'];
    load(fold_fl, 'im_mat')
    flip_nm='';
else fold_fl=['../../design_stims/stims/exp_categ1', '_', cond_nm, '_flip.mat'];
    load(fold_fl, 'im_mat_flip')
    im_mat=im_mat_flip;
    flip_nm='_flip';
end
[sz, im_n]=size(im_mat)
sz_comp=sz/3

perm_n=10000;

%%%get im size and mask ind (needed to reconstr ims of vects)
[sz_im, ones_ind]=reverse_ellipse_mask;
%sz=size(ones_ind, 1)

    
cform_lab2srgb = makecform('lab2srgb'); 

fold_sols='MDS_sols/';

% fold_rnk='rnk_CI_perms/';
% [~, ~]=mkdir(fold_rnk);
fold_pval='pval_CI_perms/';
[~, ~]=mkdir(fold_pval);

tic
for ROI_k=[8 9]%11
    
    %ROI_k=ROI_k
    
    MDS_sol_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_L1out.mat'];
    load(MDS_sol_fl, 'Y_mat')
    
    dim_max=20;%size(Y_mat, 2); % #dims for recon purposes
     

%     rnk_CI_neut_mat=NaN(sz, dim_max, 60);
%     rnk_CI_hap_mat=NaN(sz, dim_max, 60);
    pval_CI_neut_mat=NaN(sz, dim_max, 60);
    pval_CI_hap_mat=NaN(sz, dim_max, 60);
    
    for ind_k=1:60
         
         %ind_k=ind_k
         ind_train=setdiff(1:60, ind_k);
         
         im_orig_neut=im_mat(:, (ind_k*2)-1);
         im_orig_hap=im_mat(:, ind_k*2);
         im_orig_neut_train=im_mat(:, (ind_train*2)-1);
         im_orig_hap_train=im_mat(:, ind_train*2);
         


         Y_curr=Y_mat(:,:, ind_k);
         Y_L1out=Y_curr(ind_train,:);
         
        
         
         for dim_k=1:dim_max
             
             ROI_k=ROI_k
             ind_k=ind_k
             dim_k=dim_k
             
             
             CI_mat_neut=NaN(sz, perm_n);
             CI_mat_hap=NaN(sz, perm_n);
             
             parfor perm_k=1:perm_n
             
                 Y_rnd=Y_L1out(:,dim_k);%(:,1)
                 if perm_k>1

%                      cnd=1;
%                      while cnd
                          ind_rand=randperm(59)';
%                          [rho, pval]=corr(ind_rand, [1:59]');
%                          cnd=pval<0.05;
%                      end
                     Y_rnd=Y_L1out(ind_rand,dim_k);%(:,1)
                      

                 end
             
                 ind_pos=find(Y_rnd(:,1)>0);
                 ind_neg=find(Y_rnd(:,1)<0);


                 im_mat_pos_neut=im_mat(:, (ind_pos*2)-1);
                 im_mat_neg_neut=im_mat(:, (ind_neg*2)-1);
                 im_mat_pos_hap=im_mat(:, ind_pos*2);
                 im_mat_neg_hap=im_mat(:, ind_neg*2);

                 %Y_pos_scale=Y(ind_pos,dim_k)/sum(Y(ind_pos,dim_k));
                 %Y_neg_scale=-Y(ind_neg,dim_k)/sum(-Y(ind_neg,dim_k));
                 Y_pos=Y_rnd(ind_pos,1);
                 Y_neg=-Y_rnd(ind_neg,1);
                          

             
                 Y_pos_mat=repmat(Y_pos', [sz 1]);             
                 Y_neg_mat=repmat(Y_neg', [sz 1]);

                 %%%prots - unscaled here
                 prot_pos_neut=sum(im_mat_pos_neut.*Y_pos_mat, 2);
                 prot_neg_neut=sum(im_mat_neg_neut.*Y_neg_mat, 2);
                 prot_pos_hap=sum(im_mat_pos_hap.*Y_pos_mat, 2);
                 prot_neg_hap=sum(im_mat_neg_hap.*Y_neg_mat, 2);

                 CI_neut=prot_pos_neut-prot_neg_neut;
                 CI_hap=prot_pos_hap-prot_neg_hap;
                 
             

                 CI_mat_neut(:,perm_k)=CI_neut;
                 CI_mat_hap(:,perm_k)=CI_hap;
                 
                 %%%view raw CI
%                  view_on=0;
%                  if view_on && perm_k==1
%                      view_neut=1; %1 neut; 0 hap
%                      for cond_k=1:3
%                         im_vect=zeros(sz_im(:,1)*sz_im(:,2), 1);
%                         im=zeros(sz_im(1), sz_im(2));
% 
%                         if view_neut
%                             dt_comp=CI_neut((cond_k-1)*sz_comp+1:cond_k*sz_comp,1);
%                         else dt_comp=CI_hap((cond_k-1)*sz_comp+1:cond_k*sz_comp,1);
%                         end
% 
%                         mx=max(abs(dt_comp));
%                         im_vect(ones_ind, 1)=0.5+(dt_comp/mx)/2;
% 
%                         im(:,:)=reshape(im_vect, sz_im);
%                         imtool(im)
%                      end
%                  end
                 
                
                
             end
             
             toc
             
             %rnk_vect_neut=comp_rnk(CI_mat_neut, sz, perm_n);
             pval_vect_neut=comp_pval(CI_mat_neut, sz, perm_n);
             %sum(pval_vect_neut<0.05)
             %rnk_vect_hap=comp_rnk(CI_mat_hap, sz, perm_n);
             pval_vect_hap=comp_pval(CI_mat_hap, sz, perm_n);
             %sum(pval_vect_hap<0.05)
             
             %rnk_CI_neut_mat(:, dim_k, ind_k)=rnk_vect_neut;
             %rnk_CI_hap_mat(:, dim_k, ind_k)=rnk_vect_hap;
             pval_CI_neut_mat(:, dim_k, ind_k)=pval_vect_neut;
             pval_CI_hap_mat(:, dim_k, ind_k)=pval_vect_hap;    

%%computes signif based on permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             %q=0.1;
%              for cond_k=1:3
%                  
%                  cond_rng=(cond_k-1)*sz_comp+1:cond_k*sz_comp;
%                  
% %                  [pID,pN, pID_cnt, pN_cnt] = FDR_comp(pval_vect(cond_rng,1), q);
% %                  pID=pID
% %                  pID_cnt=pID_cnt
%                  
%                  pID=0.01;
%                  
% %                  tmp=sort(pval_vect(cond_rng,1));
% %                  tmp_trim=tmp(1:100,1)
% %                  hist(tmp)
% %                  error
%                  
%                  CI_neut_anl=0.5+(CI_mat_neut(cond_rng,1)-0.5)*0; %0.6
%                  
%                  
%                  if size(pID, 1)>0
%                  
%                     pval_vect_thr=double(pval_vect(cond_rng,1)<pID);
%                     sum(pval_vect_thr)
%                                             
% 
%                     CI_mat_neut_pol=double(CI_mat_neut(cond_rng,1)>0);%polarity
%                     CI_neut_anl=CI_neut_anl.*(1-pval_vect_thr)+(pval_vect_thr.*CI_mat_neut_pol);
%                  end
%                  
%                  
%                  %%%%display CI anl
%                  view_on=1;
%                  if view_on
%                      
%                     view_neut=1; %1 neut; 0 hap
%                      
%                     im_vect=zeros(sz_im(:,1)*sz_im(:,2), 1);
%                     im=zeros(sz_im(1), sz_im(2));
% 
%                     if view_neut
%                         dt_comp=CI_neut_anl;
%                     else dt_comp=CI_hap_anl;
%                     end
% 
%                     %mx=max(abs(dt_comp));
%                     %im_vect(ones_ind, 1)=0.5+(dt_comp/mx)/2;
% 
%                     im_vect(ones_ind, 1)=dt_comp;
%                     im(:,:)=reshape(im_vect, sz_im);
%                     imtool(im)
%                      
%                  end
                 
                 
%              end
             
           
             
         end
         
             
             
         
    end
    
%     rnk_CI_neut_mat=single(rnk_CI_neut_mat);
%     rnk_CI_hap_mat=single(rnk_CI_hap_mat);
    pval_CI_neut_mat=single(pval_CI_neut_mat);
    pval_CI_hap_mat=single(pval_CI_hap_mat);
    
%     rnk_neut_fl=[fold_rnk, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIrnk_neut.mat'];
%     save(rnk_neut_fl, 'rnk_CI_neut_mat')
%     rnk_hap_fl=[fold_rnk, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIrnk_hap.mat'];
%     save(rnk_hap_fl, 'rnk_CI_hap_mat') 
    pval_neut_fl=[fold_pval, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIpval_neut.mat'];
    save(pval_neut_fl, 'pval_CI_neut_mat')
    pval_hap_fl=[fold_pval, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIpval_hap.mat'];
    save(pval_hap_fl, 'pval_CI_hap_mat')
    
    
    
end

end


%%computes signif based on permutation test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function rnk_vect=comp_rnk(CI_mat, sz, perm_n)
function pval_vect=comp_pval(CI_mat, sz, perm_n)
             
 tmp_sort=sort(abs(CI_mat), 2, 'descend');
 tmp_act=abs(CI_mat(:,1));
 rnk_vect=NaN(sz,1);
 for val_k=1:sz
    rnk_vect(val_k, 1)=find(tmp_act(val_k,1)==tmp_sort(val_k,:), 1, 'last');
 end
 pval_vect=rnk_vect/perm_n;
 pval_vect=pval_vect/2;%%1-tailed


end
         
         
%%uses stack_stims_md ellipse constr to map back vects to ims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sz_im, ones_ind]=reverse_ellipse_mask
mn_iod=60;%round(mean(iod))%80%32%!!!keep orig 80 if resizing (since L is based on it)
a=2.25;
b=3.4;
L=mn_iod*10.5;%14; %5.6
ell_templ=design_ellipse(a, b, L);
ell_templ=single(ell_templ);
    %%%%%%%%%%%!!!resize 0.4 for modelling stims
     %ell_templ=round(imresize(ell_templ, 0.4));
     
     
ell_templ=round((ell_templ+fliplr(ell_templ))/2);
ell_templ=round((ell_templ+flipud(ell_templ))/2);
%imtool(ell_templ)
%error     
ell_templ_vect=logical(ell_templ(:));
ones_ind=find(ell_templ);
sz=sum(ell_templ_vect);
sz_im=size(ell_templ)        
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
         
%%%%design ellipse mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ell_templ] = design_ell_mask
    
mn_iod=60;
a=2.25;
b=3.4;
L=mn_iod*10.5;%14; %5.6
ell_templ=design_ellipse(a, b, L);
ell_templ=single(ell_templ);
    %%%%%%%%%%%!!!resize 0.4 for modelling stims
     %ell_templ=round(imresize(ell_templ, 0.4));
     
     
ell_templ=round((ell_templ+fliplr(ell_templ))/2);
ell_templ=round((ell_templ+flipud(ell_templ))/2);
ell_templ_vect=ell_templ(:)>0;
ell_templ_pixn=sum(sum(ell_templ))
% imtool(ell_templ)
% size(ell_templ)

%ell_templ_col=uint8(repmat(ell_templ, [1 1 3]));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = conv_im_RGB(im_lab, sz_im, ones_ind, cform_lab2srgb)

sz_comp=size(ones_ind, 1);
im_vect=zeros(sz_im(:,1)*sz_im(:,2), 3);
im=zeros(sz_im(1), sz_im(2), 3);
for cond=1:3
    im_vect(ones_ind, cond)=im_lab((cond-1)*sz_comp+1:cond*sz_comp, 1);
    im(:,:,cond)=reshape(im_vect(:, cond), sz_im);
end

%im=uint8(im);

im_back = applycform(im,cform_lab2srgb);
%imtool(im_back)
imtool(imresize(uint8(im_back*255), 2))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    