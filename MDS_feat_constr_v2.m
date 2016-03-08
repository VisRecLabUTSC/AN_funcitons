
%%%implemented selective use of dimensions (based on perm test with CIs)
%%%will need conv to Lab (from RGB image combination)
function MDS_feat_constr_v2


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




%%%get im size and mask ind (needed to reconstr ims of vects)
[sz_im, ones_ind]=reverse_ellipse_mask;
%sz=size(ones_ind, 1)

    
cform_lab2srgb = makecform('lab2srgb'); 




fold_sols='MDS_sols/';
fold_CIsel='pval_CI_sel/';
recon_fold='recon_res_v2/';
[~,~]=mkdir(recon_fold);

for ROI_k=2:11
    
    ROI_k=ROI_k
    
    MDS_sol_fl=[fold_sols, 'ROI', sprintf('%02.0f', ROI_k), '_MDS_L1out.mat'];
    load(MDS_sol_fl, 'Y_mat')
    
    CIsel_neut_fl=[fold_CIsel, 'ROI', sprintf('%02.0f', ROI_k), '_col_CIsel_neut.txt'];
    CIsel_hap_fl=[fold_CIsel, 'ROI', sprintf('%02.0f', ROI_k), '_col_CIsel_hap.txt'];
    
    CIsel_neut_bin=dlmread(CIsel_neut_fl);%binary matrix
    CIsel_hap_bin=dlmread(CIsel_hap_fl);
    
        %%%find if no dims provide info (by perm test); select the 1st one then
        replace_ind=sum(CIsel_neut_bin, 2)==0;
        CIsel_neut_bin(replace_ind, 1)=1;
        replace_ind=sum(CIsel_hap_bin, 2)==0;
        CIsel_hap_bin(replace_ind, 1)=1;
    
    dim_max=20;%size(Y_mat, 2); % #dims for recon purposes
    
    recon_mat=NaN(sz, im_n);
    
    for ind_k=1:60
         
         ind_k=ind_k
         ind_train=setdiff(1:60, ind_k);
         
         im_orig_neut=im_mat(:, (ind_k*2)-1);
         im_orig_hap=im_mat(:, ind_k*2);
         im_orig_neut_train=im_mat(:, (ind_train*2)-1);
         im_orig_hap_train=im_mat(:, ind_train*2);
         

         %%%find informative dimensions
         CIsel_neut_curr=find(CIsel_neut_bin(ind_k,:));
         CIsel_hap_curr=find(CIsel_hap_bin(ind_k,:));
         
         dim_sel_neut=size(CIsel_neut_curr, 2);
         dim_sel_hap=size(CIsel_hap_curr, 2);
         
         %%%select coefs for training images
         Y_curr=Y_mat(:,:, ind_k);
%          Y_L1out=Y_curr(ind_train,:);
         
         Y_L1out_sel_neut=Y_L1out(:, CIsel_neut_curr);
         Y_L1out_sel_hap=Y_L1out(:, CIsel_hap_curr);
         
         %%%find the dist to origin for each training face for the purpose of generating 'origin' face 
         %%%use only diagnostic dims
         dist_L1out_sel_neut=sqrt(sum(Y_L1out_sel_neut(:, 1:dim_sel_neut).^2, 2));
         dist_L1out_sel_hap=sqrt(sum(Y_L1out_sel_hap(:, 1:dim_sel_hap).^2, 2));
%          dist_L1out_sc=dist_L1out.^(-1);
%          dist_L1out_sc=dist_L1out_sc.*(1/sum(dist_L1out_sc, 1));
         %%%norm to unit for the purpose of generating 'origin' face
         dist_L1out_sel_neut_sc=dist_L1out_sel_neut.*(1/sum(dist_L1out_sel_neut, 1));
         dist_L1out_sel_hap_sc=dist_L1out_sel_hap.*(1/sum(dist_L1out_sel_hap, 1));
         
%          sum(dist_L1out_sel_neut_sc) %norm to 1 check
%          sum(dist_L1out_sel_hap_sc)
%          error
         
         %%%face 'origin' constr here
         im_mn_neut=sum(im_orig_neut_train .* repmat(dist_L1out_sel_neut_sc', [sz 1]), 2);
         im_mn_hap=sum(im_orig_hap_train .* repmat(dist_L1out_sel_hap_sc', [sz 1]), 2);
%          im_mn_neut=mean(cat(2, im_mat_pos_neut, im_mat_neg_neut), 2);
%          im_mn_hap=mean(cat(2, im_mat_pos_hap, im_mat_neg_hap), 2);
         
         
         
%          CI_mat_neut=NaN(sz1, sz2, sz3, dim_max);
%          CI_mat_hap=NaN(sz1, sz2, sz3, dim_max);
         CI_mat_neut=NaN(sz, dim_sel_neut);
         CI_mat_hap=NaN(sz, dim_sel_hap);
         
         %%%neut face constr
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for dim_k=1:dim_sel_neut
                          
             ind_pos=find(Y_L1out_sel_neut(:,dim_k)>0);
             ind_neg=find(Y_L1out_sel_neut(:,dim_k)<0);
             
             im_mat_pos_neut=im_orig_neut_train(:, ind_pos);
             im_mat_neg_neut=im_orig_neut_train(:, ind_neg);
             
             Y_pos=Y_L1out_sel_neut(ind_pos,dim_k);
             Y_neg=-Y_L1out_sel_neut(ind_neg,dim_k);
                                     
             Y_pos_mat=repmat(Y_pos', [sz 1]);             
             Y_neg_mat=repmat(Y_neg', [sz 1]);
             
             %%%prots - unscaled here
             prot_pos_neut=sum(im_mat_pos_neut.*Y_pos_mat, 2);
             prot_neg_neut=sum(im_mat_neg_neut.*Y_neg_mat, 2);
             
             CI_neut=prot_pos_neut-prot_neg_neut;
             
             
             cf=Y_curr(ind_k, dim_k);
             
             CI_mat_neut(:, dim_k)=cf*CI_neut/2;
             
         end
         
                  
         recon_im_neut=im_mn_neut+sum(CI_mat_neut, 2);
         recon_mat(:,(ind_k*2)-1)=recon_im_neut;
         
        % check to assess recon appearance         
        %          conv_im_RGB(im_orig_neut, sz_im, ones_ind, cform_lab2srgb)
        %          conv_im_RGB(im_mn_neut, sz_im, ones_ind, cform_lab2srgb)
        %          conv_im_RGB(recon_im_neut, sz_im, ones_ind, cform_lab2srgb)
        %          error
        
        
         %%%hap face constr
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for dim_k=1:dim_sel_hap
                          
             ind_pos=find(Y_L1out_sel_hap(:,dim_k)>0);
             ind_neg=find(Y_L1out_sel_hap(:,dim_k)<0);
             
             im_mat_pos_hap=im_orig_hap_train(:, ind_pos);
             im_mat_neg_hap=im_orig_hap_train(:, ind_neg);
             
             Y_pos=Y_L1out_sel_hap(ind_pos,dim_k);
             Y_neg=-Y_L1out_sel_hap(ind_neg,dim_k);
                                     
             Y_pos_mat=repmat(Y_pos', [sz 1]);             
             Y_neg_mat=repmat(Y_neg', [sz 1]);
             
             %%%prots - unscaled here
             prot_pos_hap=sum(im_mat_pos_hap.*Y_pos_mat, 2);
             prot_neg_hap=sum(im_mat_neg_hap.*Y_neg_mat, 2);
             
             CI_hap=prot_pos_hap-prot_neg_hap;
             
             
             cf=Y_curr(ind_k, dim_k);
             
             CI_mat_hap(:, dim_k)=cf*CI_hap/2;
             
         end
         
         recon_im_hap=im_mn_hap+sum(CI_mat_hap, 2);
         recon_mat(:,ind_k*2)=recon_im_hap;
         
        %%%check to assess recon appearance         
%          conv_im_RGB(im_orig_hap, sz_im, ones_ind, cform_lab2srgb)
%          conv_im_RGB(im_mn_hap, sz_im, ones_ind, cform_lab2srgb)
%          conv_im_RGB(recon_im_hap, sz_im, ones_ind, cform_lab2srgb)
%          error
         
                    
    end
    
    recon_res_fl=[recon_fold, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_dimselqp10', '_recon', flip_nm, '.txt'];
    dlmwrite(recon_res_fl, recon_mat) 
    
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
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    