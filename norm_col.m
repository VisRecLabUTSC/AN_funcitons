function norm_col(fold)
keepNames=1;
namesfile= [fold,'/','filenames.txt'];
fold_crop=[fold,'/norm_faces'];%_res04 %[fold,'/norm_faces/','crop_faces_a225b34L105iod60'];
fold_col=[fold,'/','norm_faces_col_ALL'];%_res04
[status,message,messageid] = mkdir(fold_col);
%coords= dlmread([fold,'/', 'coord.txt']);%only used to count files

% cform_srgb2lab = makecform('srgb2lab');
% cform_lab2srgb = makecform('lab2srgb');


mn_iod=60;%80;

a=2.25;%2.25%2.5 %2.27 (EEG exp)
b=3.4;%3.4%3.34 %3.25 %3.30(EEG exp)
L=mn_iod*10.5;%10.5 %10.4 %10.6 (EEG exp)
ell_templ=design_ellipse(a, b, L);
ell_templ=single(ell_templ);
%     %%%%%%%%%%!!!resize 0.4 for modelling stims
%      ell_templ=round(imresize(ell_templ, 0.4));
     
ell_templ=round((ell_templ+fliplr(ell_templ))/2);
ell_templ=round((ell_templ+flipud(ell_templ))/2);    
ell_templ=logical(ell_templ);

[sz1, sz2]=size(ell_templ);
zero_vect=zeros(sz1*sz2,1);
%imtool(ell_templ)

%ind=find(ell_templ);
ind=ell_templ>0;
fid=fopen(namesfile);

for k=1:118%size(coords,1)%120
    
    U{k}=fgetl(fid);
    im=imread([fold_crop, '/', U{k}(1:end-4),'.tif']);
    
        
%         im_sum=sum(im, 3);
%         im_bnr=im_sum<250*3;
%     %     imtool(im_bnr)
%         im_bnr=im_bnr(:);
%           ind=find(im_bnr);

    
%obsolete, M.Ruzon code
%    im_lab=RGB2Lab(im); 
%    class(im_lab)
%    im_rgb=Lab2RGB(im_lab);
%    class(im_rgb)
   
   im_lab=rgb2lab(im); %predefined matalb function
   
                                   %%reconversion test
                                 %  im_rgb=lab2rgb(im_lab); %predefined matalb function
                                   %im_rgb=uint8(im_rgb*255);
                                %    
                                %    imtool(im_rgb)
                                % 
                                %    tst=double(im)-double(im_rgb); 
                                %    tst_im=abs(tst)>0;
                                %    for k=1:3
                                %     imtool(tst_im(:,:,k))
                                %    end
                                %    
                                %    tst=tst(:); hist(tst)
                                %    error
                                    %im_lab=applycform(im, cform_srgb2lab);
                                    %error
    im_lab=single(im_lab);

    for compn=1:3
        im_compn=im_lab(:,:, compn);

        im_compn=im_compn(:);
        im_compn=im_compn(ind);

        mn_mat(k, compn)=mean(im_compn);
        std_mat(k, compn)=std(im_compn);
        
    end
   
end

fclose(fid)

% mn_mat
% std_mat

mns=mean(mn_mat)
stds=mean(std_mat)
%error

        %old: params for 80 iod
        %mns=[74.1369    9.0508   14.6381];%means for AR set or:   mean(mn_mat)
        %stds=[9.1358    3.1228    4.1022];%mean(std_mat)
        %params for 60 iod (For MEMORY RECON)
        %mns=[74.1236    9.0272   14.6300];%means for AR set or:   mean(mn_mat)
        %stds=[9.0889    3.1079    4.0664];%mean(std_mat)
%mns=[52.4985    17.6935    22.8290]; %values for new memory reconstruction famous faces
%stds=[13.2175    5.3788    6.2662]; % Found the mean and stds values using the old unfamiliar faces 
mns=[52.4982    17.6790    22.8058];
stds=[13.2581    5.4385    6.2295];

%%%%....
%Ash used below values to color normalize final stimuli for EEG recon.
%project 
% mns=[74.6523    8.5256   13.6698];
% stds=[8.5466    2.9109    3.4650];

% mns=[75.5073    9.5192   15.2534]; % New values to fix 0 pixel values
% stds=[8.2222    2.8105    3.6920];

%mns=[75.8523    10.5256   15.6698];
%stds=[9.5466    3.9109    4.4650];

ind=find(ell_templ);
fid=fopen(namesfile);

for k=1:118%size(coords,1)%120

    U{k}=fgetl(fid);
    im=imread([fold_crop, '/', U{k}(1:end-4),'.tif']);


     

    %im_lab=RGB2Lab(im);
    im_lab=rgb2lab(im);
%     im_lab=applycform(im, cform_srgb2lab);
%     im_lab=double(im_lab);


    for compn=1:3
        im_compn=im_lab(:,:, compn);
        
        im_compn=im_compn(:);
        im_compn=im_compn(ind);


        im_compn=im_compn-mn_mat(k, compn);
        im_compn=im_compn/(std_mat(k, compn)/stds(1,compn));

        im_compn=im_compn+mns(1, compn);

        
                    im_vect=zero_vect;
                    im_vect(ind)=im_compn;
                   
                    im_lab(:,:,compn)=reshape(im_vect, [sz1 sz2]);
        
        
        
        %im_lab(:,:,compn)=im_compn;
    end

    %im_rgb=Lab2RGB(im_lab);
    im_rgb=lab2rgb(im_lab);
    
    %im_rgb=applycform(im_lab, cform_lab2srgb);
    %class(im_rgb)
    
    for compn=1:3
        im_rgb(:,:, compn)=single(im_rgb(:,:, compn)+0).*single(ell_templ)+(1-single(ell_templ))*0;
    end
    
   im_rgb=uint8(im_rgb*255);
%    imtool(im_rgb)

   %imwrite(im_rgb, [fold_col, '/', U{k}(1:end-4),'.tif'])
   
   %%%%stim1...- males, stim2... - females
   %imwrite(im_rgb, [fold_col, '/stim2', sprintf('%03.0f', k),'1.tif']) %only neutral
   if keepNames
       cd(fold_crop)
       files=dir('*.tif');
       imwrite(im_rgb, [fold_col, '/' files(k).name])
   else
       imwrite(im_rgb, [fold_col, '/stim2', sprintf('%03.0f', ceil(k/2)), num2str(2-rem(k,2)),'.tif']) %1 or 2 appended for neut and happy
   end
 end

    
    
    
