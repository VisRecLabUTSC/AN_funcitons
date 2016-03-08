%store stims as vects (only cropped part) in mat
function stack_stims

col_cond=1;
if col_cond
    col_nm='col'
else col_nm='gray'
end

categ_stims=1;
fold_fl=['allstims/','norm_females_col/']; % write to this folder?
fl_nm=['allstims/norm_females_col/filenames.txt']; % names of original stims used in behavioural portion of experiment?

%%%%%%%%%%design ellipse mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mn_iod=60%round(mean(iod))%80%32%!!!keep orig 80 if resizing (since L is based on it)
a=2.25;
b=3.4;
L=mn_iod*10.5;%14; %5.6
ell_templ=design_ellipse(a, b, L);
ell_templ=single(ell_templ);
    %%%%%%%%%%%!!!resize 0.4 for modelling stims
     %ell_templ=round(imresize(ell_templ, 0.4));
     
     
ell_templ=round((ell_templ+fliplr(ell_templ))/2);
ell_templ=round((ell_templ+flipud(ell_templ))/2);
sum(sum(ell_templ))
%imtool(ell_templ)
%error     
ell_templ_vect=logical(ell_templ(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz=sum(ell_templ_vect);
if col_cond
    sz=sz*3;
end
    
im_mat=NaN(sz, 120); 
im_mat_flip=NaN(sz, 120);
disp(fl_nm)

fid=fopen(fl_nm);
disp(fid)
for fl_k=1:120
    im_nm=fgetl(fid);
    
    im=imread([fold_fl, im_nm]);
    
    
    if col_cond
        %im=RGBToLab(im);
%         imtool(im)
        im=RGB2Lab(im);
%         im=Lab2RGB(im);
%         imtool(im)
%         error
        im_flip=NaN(size(im));
        for comp=1:3
            im_flip(:,:,comp)=fliplr(im(:,:,comp));
        end
    else       
        %im=rgb2gray(im);
        im=RGB2Lab(im);
        im=im(:,:,1);
        im_flip=fliplr(im);
%         imtool(double(im).*ell_templ/255)
%         imtool(double(im_flip).*ell_templ/255)
%         error
    end
    
    if col_cond
        im_cond=NaN(sz/3, 3);
        im_cond_flip=NaN(sz/3, 3);
        for cond=1:3
            tmp=im(:,:,cond);
            tmp=tmp(:);
%             size(tmp)
            im_cond(:,cond)=tmp(ell_templ_vect);
            
            tmp=im_flip(:,:,cond);
            tmp=tmp(:);
            im_cond_flip(:,cond)=tmp(ell_templ_vect);
        end
        im=im_cond(:);
        im_flip=im_cond_flip(:);
        
    else
        im=im(:);
        im=im(ell_templ_vect);
        im_flip=im_flip(:);
        im_flip=im_flip(ell_templ_vect);
    end
    
%     size(im)
%     im(1:40)
%     error
    
    im_mat(:, fl_k)=im;
    im_mat_flip(:, fl_k)=im_flip;
    
end
    
fclose(fid)

%if col_cond
    im_mat=single(im_mat);
    im_mat_flip=single(im_mat_flip);
%end
save(['stims_mats/exp_fem_', col_nm, '.mat'], 'im_mat') %categ', num2str(categ_stims),'_res04_',
save(['stims_mats/exp_fem_', col_nm, '_flip.mat'], 'im_mat_flip') % categ', num2str(categ_stims),'_res04_',
