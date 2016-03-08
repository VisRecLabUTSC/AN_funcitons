
function CI_sel


q=0.10;
%p_thr=0.01; prop_thr=0.05; %%%just for prelim anl to assess pvals
expr_k=2; %1-neut; 2- hap


cond_col=1;
if cond_col
    cond_nm='col';
else cond_nm='gray';
end

fold_pval='pval_CI_perms/';

fold_sel='pval_CI_sel/';
[~, ~]=mkdir(fold_sel);


for ROI_k=0:11%11
    
    ROI_k=ROI_k
    
    dim_max=20;%size(Y_mat, 2); % #dims for recon purposes

    if expr_k==1
        pval_neut_fl=[fold_pval, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIpval_neut.mat'];
        load(pval_neut_fl, 'pval_CI_neut_mat')
        pval_CI_mat=pval_CI_neut_mat;
    else 
        pval_hap_fl=[fold_pval, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIpval_hap.mat'];
        load(pval_hap_fl, 'pval_CI_hap_mat')
        pval_CI_mat=pval_CI_hap_mat;
    end
    
    
%     figure
%     hold on

    sel_res_pthr_bin=zeros(60, dim_max);    
    sel_res_FDR_bin=zeros(60, dim_max);
    
    for ind_k=1:60
        
        %prop_res=mean(pval_CI_mat(:,:, ind_k)<p_thr);
        %plot(prop_res)
        %sel_res_pthr_bin(ind_k, :) = double(prop_res > prop_thr);
        %sel_res=find(sel_res_bin(ind_k, :));
        
        for dim_k=1:dim_max
        
            [pID,pN, pID_cnt, pN_cnt] = FDR_comp(pval_CI_mat(:,dim_k, ind_k), q);
            
%              pID=pID
%              pID_cnt=pID_cnt

             sel_res_FDR_bin(ind_k, dim_k)=size(pID, 1)>0;
                
                       
        end
        
    end
    
    %sel_res_pthr_bin
    sel_res_FDR_bin
    
    if expr_k==1
        sel_neut_fl=[fold_sel, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIsel_neut.txt'];
        dlmwrite(sel_neut_fl, sel_res_FDR_bin)
    else 
        sel_hap_fl=[fold_sel, 'ROI',  sprintf('%02.0f', ROI_k), '_', cond_nm, '_CIsel_hap.txt'];
        dlmwrite(sel_hap_fl, sel_res_FDR_bin)
    end
    
end

        