
%%%compute d' for all possible pairs os stims (within / across 2 subsets)
%%%saves res separately for each vox like eval_voxROI_MA
function []=searchlight_MA_full(subj, ball_r, ball_nm, rnd_cond, univar_cond)

if nargin==0

    subj=2;
    ball_r=5;%2.2;% %in #vox
    ball_nm='r5';%'md';%'sm';%
    rnd_cond=0;
    univar_cond=false;%true;
    
end
%rnd_cond=rnd_cond
dt_type_vect={'F', 'E', 'FE'};%'FN' 'F' 'N'
subjid=['S', sprintf('%02.0f', subj)];
univar_str='';
if univar_cond
    univar_str='_univar';
end

concat_opt=0; %0 for average across TP, 1 for concat

%read masks (mask_var computed from nondet data)
mask_ind_nm=['data_det/',subjid, '/mask_ind.txt'];
mask_ind=dlmread(mask_ind_nm);
coord_mat=mask_ind(:,2:4);

mask_var_nm=['data_det/',subjid, '/var_mask.txt'];
mask_var=dlmread(mask_var_nm);

%stim_subsets=dlmread(['data_det/id_runassign.txt']);


for dt_type_k=2:2
    dt_type=char(dt_type_vect(dt_type_k))
    
    switch dt_type_k
        case 1            
           
            idpair_list=nchoosek(1:60,2);
            pair_n=1770;
        case 2
            
            idpair_list=nchoosek(1:120,2);
            tmp=rem(idpair_list, 2);
            diff_expr_ind=tmp(:,1)~=tmp(:,2);
            diff_id_ind=(idpair_list(:,1)~=(idpair_list(:,2)-1)) | tmp(:,1)==0;
            diff_ind=diff_expr_ind&diff_id_ind;
            idpair_list=idpair_list(diff_ind, :);
%             size(idpair_list)
%             idpair_list(1:10,:)
%             idpair_list(end-9:end, :)
            pair_n=3540;            
        case 3
            error
           
            idpair_list=[1 2];
            pair_n=1;
    end

    
    for startTP=7:9
        
        fl_nm=['data_det/', subjid, '/', subjid, '_', dt_type, '_TP', num2str(startTP), '_zed.mat'];
        load(fl_nm)
%                                     if rnd_cond==1
%                                         whole_patt_rnd(:,1,:)=cat(3, whole_patt_zed(:,1, ind11), whole_patt_zed(:, 2, ind12));
%                                         whole_patt_rnd(:,2,:)=cat(3, whole_patt_zed(:,1, ind21), whole_patt_zed(:, 2, ind22));
%                                         whole_patt_zed=whole_patt_rnd;
%         %                                 size(whole_patt_zed)
%                                     end


        if startTP==7
            V_mat=whole_patt_zed;
        else V_mat=cat (4, V_mat, whole_patt_zed);
        end 
        clear whole_patt_zed
    end
    
    %%%normalize within voxel X TP to 0-1 (and elim outliers +- 3SD)
    for TP_k=1:3
        V_mat_reshape=reshape(V_mat(:,:,:, TP_k), size(V_mat, 1), size(V_mat, 2)*size(V_mat, 3));
        V_mat_reshape_zed=zscore(V_mat_reshape, 0, 2);
        V_mat_reshape_norm=(V_mat_reshape_zed/6)+0.5; % scale to 0-1 and fix outs next

        V_mat_reshape_posout=V_mat_reshape_norm>1;
        V_mat_reshape_negout=V_mat_reshape_norm<0;
        V_mat_reshape_in=~V_mat_reshape_posout & ~V_mat_reshape_negout;

        V_mat_reshape_norm=V_mat_reshape_norm.*double(V_mat_reshape_in)+double(V_mat_reshape_posout); %replace max outs w/ 1, leave min ones as 0
        V_mat(:,:,:, TP_k)=reshape(V_mat_reshape_norm, size(V_mat, 1), size(V_mat, 2), size(V_mat, 3));

        %check_n_posout=sum(sum(single(V_mat_reshape_posout)))
        %check_n_negout=sum(sum(single(V_mat_reshape_negout)))
        %check_n_in=sum(sum(single(V_mat_reshape_in)))
    end

    %error
    

    tic
    

    vox_n=size(mask_ind, 1);
    
    %%%ind of voxs of interest
    %vox_set=[116872 144899 107385 ... ];%%S01
    
    for vox_k=1:vox_n %16510% %S01:16510:33020:49523

        
%         if floor(coord_mat(vox_k, 1)/3)==(coord_mat(vox_k, 1)/3) && floor(coord_mat(vox_k, 2)/3)==(coord_mat(vox_k, 2)/3) &&...
%         floor(coord_mat(vox_k, 3)/3)==(coord_mat(vox_k,3)/3) &&...
%         coord_mat(vox_k, 3)>9 && coord_mat(vox_k, 3)<20 %%%split by z coord across nodes

          %if ismember(mask_ind(vox_k,1), vox_set)
    
            dist_mat=coord_mat-repmat(coord_mat(vox_k,:), [size(coord_mat, 1) 1]);
            dist_mat=sum(dist_mat.^2, 2);
            dist_mat=dist_mat<=(ball_r^2);
            ball_ind=find(dist_mat .* mask_var);

            ball_sz=size(ball_ind, 1);
            
            if ball_sz>40
                
                vox_k=vox_k
                
                
                vox_ind=mask_ind(vox_k,1);
                vox_nm=sprintf('%06.0f', vox_ind);
                if rnd_cond~=1
                    fold_perf=['results_svmlib', univar_str, '_MA_full/', subjid, '/'];
                    [jnk1, jnk2]=mkdir(fold_perf);
                    perf_fl_nm=[fold_perf, subjid, '_MA_ind', vox_nm, '_', dt_type, '_ball', ball_nm, '_perf_svm.txt'] %MA_sel
                %else perf_fl_nm_root=['results', univar_str, '/', subjid, '_', dt_type, '_ball', ball_nm, '_perf_rnd'];
                end
                dlmwrite(perf_fl_nm, zeros(0, 2), 'precision', '% 4.4f')%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!comment out if appending to file
                
                

                V_sel=V_mat(ball_ind, :, :, :);
                
                if concat_opt
                    V_sel=cat(1, V_sel(:,:,:,1), V_sel(:,:,:,2), V_sel(:,:,:,3));
                else V_sel=mean(V_sel, 4);
                end
                


                if univar_cond
                    V_sel=mean(V_sel, 1);
                end

                V_sel=double(V_sel);
                
            
            
                ap_perf=NaN(pair_n,1);
                dp_perf=NaN(pair_n,1);
    

                        
                parfor pair_k=1:pair_n    %parfor


                        idpair=idpair_list(pair_k,:);

                        [ap, dp]=apply_classf_libsvm_MA_full(V_sel, idpair)
                        %[ap dp MI c_mn c_std]=apply_classf_libsvm_MA_full(V_sel, univar_cond, idpair);
                        ap_perf(pair_k,:)=ap;
                        dp_perf(pair_k,:)=dp;


                end


                dlmwrite(perf_fl_nm, [ap_perf dp_perf],  'precision', '% 4.4f')
            end
            %toc
          
    end 
    toc
        
end

