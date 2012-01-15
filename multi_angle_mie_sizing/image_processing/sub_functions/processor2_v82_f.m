function processor2_v82_f(filename)
% Process the given processing initialization file, or
% if filename is a number, open a GUI and set
% skip_dialog=0, ask for .ini file, pause between angles
% skip_dialog=1, ask for .ini file, do all angles without pause in between
% skip_dialog=2, do all angles without stopping in between
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

warning('off','Images:initSize:adjustingMag');
close all

if isnumeric(filename)==1
    skip_dialog=filename;
    current_path=pwd;   % gets current path for use in dialogs
    [data_file, current_path, filterindex] = uigetfile('*.ini',...
        'Select .ini file:', current_path ,'MultiSelect', 'Off');
    filename=strcat(current_path,data_file);
else
    skip_dialog=2;
    result=strfind(filename,'/');
    if numel(result)==0
        current_path=strcat(pwd,'/');
        data_file=filename;
    else
        lst_sl=result(size(result,2));
        current_path=filename(1:lst_sl);
        data_file=filename(lst_sl+1:length(filename));
    end
end

[keys,sections,subsections] = inifile(filename,'readall');

key_assignment_string=strcat(keys(:,2),'_',keys(:,3),'=',keys(:,4),';');

% Assign values to the initalization file variables
for ii=1:size(key_assignment_string,1)
    eval(char(key_assignment_string(ii)))
end

% Determine how many data angles are to be processed
num_angles=0;
for ii=1:size(subsections,1)
    if sum(strcmp(subsections{ii,2}(1:5),{'angle' 'Angle'}))>0
        num_angles=num_angles+1;
    end
end
if num_angles==0
    error('Check INI file, must be at least one "Angle" entry')
end

clear filterindex ii job_desc keys pathname sections
clear subsections key_assignment_string

%% prepare ouput directory, file
% information will be stored here during execution (things that were 
% calculated/adjusted during processing
% assures that save dir is available and working correctly before
% calculations begin

if output_save_dir(size(output_save_dir,2))~='/'
    output_save_dir=strcat(output_save_dir,...
        '_',output_processing_info_identifier,'/');
else
    output_save_dir=...
        strcat(output_save_dir(1:(length(output_save_dir)-1)),...
        '_',output_processing_info_identifier,'/');
end

[success message messageid]=mkdir(output_save_dir);
if success~=1
    error(strcat('Failure to create output directory: ',message))
end

out_filename=strcat(output_save_dir,output_processing_info_identifier,...
    '_processing_log.txt');

[out_fid message]=fopen(out_filename,'w');
if out_fid<0
    error(strcat('Failure to create output file: ',message))
end

out_process_filename=strcat(output_save_dir,...
    output_processing_info_identifier,'_',data_file);
[success message messageid]=copyfile(filename, out_process_filename);
if success~=1
    error(strcat('Fail to copy .ini file to output directory: ',message))
end

% output may go to screen (1) and the above prepared file.  Any number of
% FIDS may be added to the array, the same output will go to all files in
% the array.
fids=[1 out_fid];

%% General information
% origin which will be 0,0,0 in the figures, with respect to axis definied
% in the INI file
yo=output_figure_origin(2);
zo=output_figure_origin(3);
mfprintf(fids,['Figure origin [0  y-offset  z-offset] '...
    '[ 0  %f  %f ]  \n'],yo,zo);

% find original image size, corrected 'flat' images will 
% (arbitrarily, could pick anything) be same size as original
ref_img=imread(reference_registration_image_name);
size_info=size(ref_img);
z_pix=size_info(1);
y_pix=size_info(2);

% dist between red dots will take up 0.85 of width of image
dy=test_section_red_dot_2(2)-test_section_red_dot_1(2);
dz=test_section_red_dot_3(3)-test_section_red_dot_1(3);
% pixel resolution will be determined by distance between these red dots
pix_res_y=(0.85*y_pix)/dy;
pix_res_z=(0.85*z_pix)/dz;
if pix_res_y<pix_res_z
    pix_res=pix_res_y;
else
    pix_res=pix_res_z;
end

image_width=y_pix/pix_res;
image_height=z_pix/pix_res;

delta_y=test_section_red_dot_2(2)-test_section_red_dot_1(2);
delta_z=test_section_red_dot_2(3)-test_section_red_dot_1(3);

% lower left dot (dot1) will be 10% up and 10% over from bottom left corner
dot1=[0.075*y_pix                  0.925*z_pix];

dot2=[dot1(1)+pix_res*delta_y    dot1(2)-pix_res*delta_z];

delta_y=test_section_red_dot_3(2)-test_section_red_dot_1(2);
delta_z=test_section_red_dot_3(3)-test_section_red_dot_1(3);
dot3=[dot1(1)+pix_res*delta_y    dot1(2)-pix_res*delta_z];

delta_y=test_section_red_dot_4(2)-test_section_red_dot_1(2);
delta_z=test_section_red_dot_4(3)-test_section_red_dot_1(3);
dot4=[dot1(1)+pix_res*delta_y    dot1(2)-pix_res*delta_z];

% pixel coordinates of the registration dots in the 'flat' image
flat_coords=[dot1;dot2;dot3;dot4];

% check that all red dots will be within 'flat' image
if max(flat_coords,1)>y_pix
    error(['requested Y red dot locations mis-fit with'...
        'default image resolution'])
end
if max(flat_coords,2)>z_pix
    error(['requested Z red dot locations mis-fit with'...
        'default image resolution'])
end

% "Exact" coordinates for center of each pixel
% calcualte y-z coordinates in the image,  x coordinate is 0.
for ii=1:y_pix
    ycorr(ii)=test_section_red_dot_1(2)-dot1(1)/pix_res+...
        (ii-1)/(pix_res)+0.5/pix_res;
end
for ii=1:z_pix
    zcorr(ii)=test_section_red_dot_1(3)-...
        (z_pix-dot1(2))/pix_res+(z_pix-ii)/(pix_res)+0.5/pix_res;
end

% determine the sub-image limits, indexes, and coordinates
y_inds=find(ycorr>=output_limits_y(1) & ycorr<=output_limits_y(2));
z_inds=find(zcorr>=output_limits_z(1) & zcorr<=output_limits_z(2));
ycorr_blk=ycorr(y_inds);
zcorr_blk=zcorr(z_inds);

%% **** Reference Image
mfprintf(fids,'**** Reference Image **** \n');
%% Find the red registration dots in image, and correct perspective

% data image is read,
ref_dat_img=imread(reference_data_image_name);

% If there were image distortions (other than perspective!), correct BOTH
% ref_img and ref_dat_img HERE:
        % ---> no current distortion correctons
% -------------

% dots are located in the reference registration image, coords in pixels
coords=redotlocator_v4(ref_img,0);

% scaling from the ini file
scaled_coords=[test_section_red_dot_1;...
    test_section_red_dot_2;...
    test_section_red_dot_3;...
    test_section_red_dot_4];

% camera position, as given
ref_camera_position=[reference_distance*sin(reference_angle*pi/180) ...
    reference_distance*cos(reference_angle*pi/180) reference_cam_height];

% correction for the offset between registration dots and the image plane
new_coords=find_reg_offset(coords, scaled_coords, ref_camera_position);

% transform applied against the data image
% obviously registration image and data image MUST be identical in every
% way; location, resolution etc.
[ref_img_trans_orig xdata ydata]=correct_perspective(ref_dat_img,...
    flat_coords, new_coords, y_pix, z_pix);

%% Set everyhting below the threshold to zero
ref_img_trans_orig(ref_img_trans_orig<output_threshold)=0;

%% apply camera correction
ref_img_trans=make_image_correction(...
    double(ref_img_trans_orig(:,:,output_color_channel))/255,...
    output_correction_file_name);

%% Histogram of data.  Check against too many saturated pixels.
[counts bins]=imhist(ref_img_trans_orig(:,:,output_color_channel));

if strcmp(output_output_eps_figures,'Yes')==1
    figure
    title('histogram')
    stem(bins(2:size(bins)),counts(2:size(bins)))
end

sat_percent=100*counts(size(counts,1))/sum(counts(2:size(counts,1)));
mfprintf(fids,'Percentage saturated pixels: %f \n', sat_percent);
if sat_percent>output_max_percent_sat_pixels
    mfprintf(fids,['    !!! WARNING !!! Percentage saturated pixels' ...
        'above limit of %f \n'], output_max_percent_sat_pixels);
end
%% Interpret intensities from image into XYZ coordinates.

[ref_int ref_int_blk ref_sup act_msk_y act_msk_z ref_blk ref_avg_dat]=...
   intensities2xyz_v2(ref_img_trans, reference_distance,...
   reference_data_image_exposure_time,...
   pix_res, z_pix, y_inds, z_inds, ycorr, zcorr, output_super_pixels_yz,...
   output_super_pixels_size, output_super_pixel_shape, fids);

%% Generate output data and figures
% raw data file
save_ref_data_name=strcat(output_save_dir,...
    output_processing_info_identifier,'_reference_data.mat');
save(save_ref_data_name, 'ref_img_trans', 'ycorr', 'zcorr', 'ref_int',...
    'ycorr_blk', 'zcorr_blk', 'ref_int_blk', 'act_msk_y', 'act_msk_z',...
    'ref_sup', 'ref_blk','ref_avg_dat', 'yo', 'zo');

% figures
fig_id_str='ref';

if strcmp(output_output_eps_figures,'Yes')==1
    
    figure
    plot(flat_coords(:,1),z_pix-flat_coords(:,2)+1,'r+')
    hold on
    plot(coords(:,1),z_pix-coords(:,2)+1,'b+')
    plot(new_coords(:,1),z_pix-new_coords(:,2)+1,'k+')
    legend('Flat Image','Data Image','Corrected','Location','BestOutside')
    title('Location of registration points')
    xlabel('pixels')
    ylabel('pixels')
    axis equal
    axis manual
    axis([0 y_pix 0 z_pix])
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_registration_points.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    warning('off','Images:initSize:adjustingMag');
    figure
    imshow(ref_dat_img, 'Border', 'tight')
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_original_image.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    figure
    imshow(ref_img_trans_orig, 'Border', 'tight')
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_original_flat_image.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    figure
    imshow(ref_img_trans, 'Border', 'tight')
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_corrected_flat_image.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    if strcmp(output_limit_full_res_figures,'Yes')==0
        figure
        contourf(ycorr-yo,zcorr-zo, ref_img_trans,'LineStyle','none')
        colorbar
        
        title_str=['Flat image, intensities from color channel '...
            num2str(output_color_channel)];
        title(title_str)
        xlabel('inches')
        ylabel('inches')
        axis equal
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_flat_color_channel_image.eps');
        saveas(gcf,save_fig_name,'epsc')
    end
    
    figure
    contourf(ycorr_blk-yo,zcorr_blk-zo, ref_blk,'LineStyle','none')
    colorbar
    
    title_str=['Flat image, intensities from color channel '...
        num2str(output_color_channel)];
    title(title_str)
    xlabel('inches')
    ylabel('inches')
    axis equal
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_flat_color_channel_image_sub_block.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    if strcmp(output_limit_full_res_figures,'Yes')==0
        figure
        v=0:max(max(ref_int))/30:max(max(ref_int));
        contourf(ycorr-yo,zcorr-zo, ref_int,v,'LineStyle','none')
        colorbar
        switch ischar(reference_caxis_limits)
            case 1
            otherwise
                caxis(reference_caxis_limits)
        end
        title('Normalized intensity')
        xlabel('inches')
        ylabel('inches')
        axis equal
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_normalized_intensities.eps');
        saveas(gcf,save_fig_name,'epsc')
    end
    
    figure
    v=0:max(max(ref_int_blk))/30:max(max(ref_int_blk));
    contourf(ycorr_blk-yo,zcorr_blk-zo, ref_int_blk,v,'LineStyle','none')
    colorbar
    switch ischar(reference_caxis_limits)
        case 1
        otherwise
            caxis(reference_caxis_limits)
    end
    title('Normalized intensity')
    xlabel('inches')
    ylabel('inches')
    axis equal
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_normalized_intensities_sub_block.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    figure
    v=0:max(max(ref_sup))/30:max(max(ref_sup));
    contourf(act_msk_y-yo, act_msk_z-zo,ref_sup,v,'LineStyle','none')
    colorbar
    switch ischar(reference_caxis_limits)
        case 1
        otherwise
            caxis(reference_caxis_limits)
    end
    title('Super-pixel intensities')
    xlabel('inches')
    ylabel('inches')
    axis equal
    if output_super_pixels_yz(1)<30 && output_super_pixels_yz(2)<30
        hold on
        for ii=1:output_super_pixels_yz(1)
            for jj=1:output_super_pixels_yz(2)
                plot(act_msk_y(ii)-yo,act_msk_z(jj)-zo,'k+')
            end
        end
    end
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_super_pixel_intensities_contour.eps');
    saveas(gcf,save_fig_name,'epsc')
    
    if output_super_pixels_yz(1)<30 && output_super_pixels_yz(2)<30
        figure
        v=0:max(max(ref_avg_dat))/30:max(max(ref_avg_dat));
        contourf(ycorr_blk-yo,zcorr_blk-zo, ref_avg_dat,v)
        colorbar
        switch ischar(reference_caxis_limits)
            case 1
            otherwise
                caxis(reference_caxis_limits)
        end
        title('Super-pixel intensity, area, and location')
        xlabel('inches')
        ylabel('inches')
        axis equal
        if output_super_pixels_yz(1)<20 && output_super_pixels_yz(2)<20
            hold on
            for ii=1:output_super_pixels_yz(1)
                for jj=1:output_super_pixels_yz(2)
                    plot(act_msk_y(ii)-yo,act_msk_z(jj)-zo,'k+')
                end
            end
        end
    end
    save_fig_name=strcat(output_save_dir,...
        output_processing_info_identifier,'_',fig_id_str,...
        '_super_pixel_intensities.eps');
    saveas(gcf,save_fig_name,'epsc')
    
end

%% Process Angle data
% all coordinates are assumed to be the same for all the angles, as the
% "flat" image was defined at the beginning and kept constant; all other
% angles transformed to match that image

i_ratio=zeros(output_super_pixels_yz(2), output_super_pixels_yz(1),...
    num_angles);

for kk=1:num_angles
    if skip_dialog==0
        % Continue or not?
        button = MFquestdlg([ 0.6 , 0.1 ],...
            'Continue processing next angle?'); % note, non-matlab function
        
        switch button
            case 'Yes'
                
            otherwise
                return;
        end
    end
    
    close all
    
    info_string=['**** Angle' num2str(kk) ' Image **** \n'];
    mfprintf(fids,info_string);
    % read in the registration and data images
    registration_name_string=strcat('registration_name=angle',...
        num2str(kk),'_registration_image_name;');
    eval(registration_name_string)
    data_name_string=strcat('data_name=angle',num2str(kk),...
        '_data_image_name;');
    eval(data_name_string)
    
    angle_reg_img=imread(registration_name);
    angle_dat_img=imread(data_name);
    
    % read in the initialization information for these images
    angle_string=strcat('cam_angle(kk)=angle',num2str(kk),'_angle;');
    eval(angle_string)
    distance_string=strcat('cam_distance=angle',num2str(kk),'_distance;');
    eval(distance_string)
    distance_string=strcat('exp_time=angle',num2str(kk),...
        '_data_image_exposure_time;');
    eval(distance_string)
    cam_height_string=strcat('cam_height=angle',num2str(kk),...
        '_cam_height;');
    eval(cam_height_string)
    caxis_string=strcat('caxis_limits=angle',num2str(kk),'_caxis_limits;');
    eval(caxis_string)
    
    % If there were image distortions (other than perspective!), correct
    % angle_reg_img and angle_dat_img HERE:
    % ---> no current distortion correctons
    % -------------
    
    % dots are located in the reference registration image
    coords=redotlocator_v4(angle_reg_img,0);
    
    % camera position, as given
    camera_position(kk,:)=[cam_distance*sin(cam_angle(kk)*pi/180) ...
        cam_distance*cos(cam_angle(kk)*pi/180) cam_height];
    
    % correction for the offset between registration dots and image plane
    new_coords=find_reg_offset(coords, scaled_coords, ...
        camera_position(kk,:));
    
    % Find the red registration dots in image, and correct perspective
    % transform applied against the data image
    % The registration image and data image MUST be identical in every
    % way; location, resolution etc.
    [angle_img_trans_orig]=correct_perspective(angle_dat_img,...
        flat_coords, new_coords, y_pix, z_pix);
    
    %% Set everyhting below the threshold to zero
    angle_img_trans_orig(angle_img_trans_orig<output_threshold)=0;
    
    %% Apply camera correction
    angle_img_trans=make_image_correction(...
        double(angle_img_trans_orig(:,:,output_color_channel))/255,...
        output_correction_file_name);
    
    %% Histogram of data.  Check against too many saturated pixels.
    [counts bins]=imhist(angle_img_trans_orig(:,:,output_color_channel));
    
    if strcmp(output_output_eps_figures,'Yes')==1
        figure
        title('histogram')
        stem(bins(2:size(bins)),counts(2:size(bins)))
    end
    
    sat_percent=100*counts(size(counts,1))/sum(counts(2:size(counts,1)));
    mfprintf(fids,'Percentage saturated pixels: %f \n', sat_percent);
    if sat_percent>output_max_percent_sat_pixels
        mfprintf(fids,['    !!! WARNING !!! Percentage saturated pixels'...
            'above limit of %f \n'], output_max_percent_sat_pixels);
    end
    
    %% Interpret intensities from image into XYZ coordinates.
    
    [i_exact i_blk i_sup sup_cor_y sup_cor_z img_blk sup_avg]=...
        intensities2xyz_v2(angle_img_trans, cam_distance,...
        exp_time, pix_res, z_pix, y_inds, z_inds, ycorr, zcorr,...
        output_super_pixels_yz, output_super_pixels_size, ...
        output_super_pixel_shape, fids);
    
    %% Output data
    % raw data file
    save_ang_data_name=strcat(output_save_dir,...
       output_processing_info_identifier,'_angle',num2str(kk),'_data.mat');
    save(save_ang_data_name, 'angle_img_trans', 'ycorr', 'zcorr',...
        'i_exact','ycorr_blk', 'zcorr_blk', 'i_blk', 'act_msk_y',...
        'act_msk_z', 'i_sup', 'img_blk','sup_avg', 'yo', 'zo');
    
    % figures
    fig_id_str=strcat('angle',num2str(kk));
    
    if strcmp(output_output_eps_figures,'Yes')==1
        
        figure
        plot(flat_coords(:,1),z_pix-flat_coords(:,2)+1,'r+')
        hold on
        plot(coords(:,1),z_pix-coords(:,2)+1,'b+')
        plot(new_coords(:,1),z_pix-new_coords(:,2)+1,'k+')
        legend('Flat Image','Data Image','Corrected','Location',...
            'BestOutside')
        title('Location of registration points')
        xlabel('pixels')
        ylabel('pixels')
        axis equal
        axis manual
        axis([0 y_pix 0 z_pix])
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_registration_points.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        warning('off','Images:initSize:adjustingMag');
        figure
        imshow(angle_dat_img, 'Border', 'tight')
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_original_image.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        figure
        imshow(angle_img_trans_orig, 'Border', 'tight')
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_original_flat_image.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        figure
        imshow(angle_img_trans, 'Border', 'tight')
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_corrected_flat_image.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        if strcmp(output_limit_full_res_figures,'Yes')==0
            figure
            contourf(ycorr-yo,zcorr-zo, angle_img_trans,'LineStyle','none')
            colorbar
            
            title_str=['Flat image, intensities from color channel '...
                num2str(output_color_channel)];
            title(title_str)
            xlabel('inches')
            ylabel('inches')
            axis equal
            save_fig_name=strcat(output_save_dir,...
                output_processing_info_identifier,'_',fig_id_str,...
                '_flat_color_channel_image.eps');
            saveas(gcf,save_fig_name,'epsc')
        end
        
        figure
        contourf(ycorr_blk-yo,zcorr_blk-zo, img_blk,'LineStyle','none')
        colorbar
        
        title_str=['Flat image, intensities from color channel '...
            num2str(output_color_channel)];
        title(title_str)
        xlabel('inches')
        ylabel('inches')
        axis equal
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_flat_color_channel_image_sub_block.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        if strcmp(output_limit_full_res_figures,'Yes')==0
            figure
            v=0:max(max(i_exact))/30:max(max(i_exact));
            contourf(ycorr-yo,zcorr-zo, i_exact,v,'LineStyle','none')
            colorbar
            switch ischar(caxis_limits)
                case 1
                otherwise
                    caxis(caxis_limits)
            end
            title('Normalized intensity')
            xlabel('inches')
            ylabel('inches')
            axis equal
            save_fig_name=strcat(output_save_dir,...
                output_processing_info_identifier,'_',fig_id_str,...
                '_normalized_intensities.eps');
            saveas(gcf,save_fig_name,'epsc')
        end
        
        figure
        v=0:max(max(i_blk))/30:max(max(i_blk));
        contourf(ycorr_blk-yo,zcorr_blk-zo, i_blk,v,'LineStyle','none')
        colorbar
        switch ischar(caxis_limits)
            case 1
            otherwise
                caxis(caxis_limits)
        end
        title('Normalized intensity')
        xlabel('inches')
        ylabel('inches')
        axis equal
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_normalized_intensities_sub_block.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        figure
        v=0:max(max(i_sup))/30:max(max(i_sup));
        contourf(act_msk_y-yo, act_msk_z-zo,i_sup,v,'LineStyle','none')
        colorbar
        switch ischar(caxis_limits)
            case 1
            otherwise
                caxis(caxis_limits)
        end
        title('Super-pixel intensities')
        xlabel('inches')
        ylabel('inches')
        axis equal
        if output_super_pixels_yz(1)<30 && output_super_pixels_yz(2)<30
            hold on
            for ii=1:output_super_pixels_yz(1)
                for jj=1:output_super_pixels_yz(2)
                    plot(act_msk_y(ii)-yo,act_msk_z(jj)-zo,'k+')
                end
            end
        end
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_super_pixel_intensities_contour.eps');
        saveas(gcf,save_fig_name,'epsc')
        
        if output_super_pixels_yz(1)<30 && output_super_pixels_yz(2)<30
            figure
            v=0:max(max(sup_avg))/30:max(max(sup_avg));
            contourf(ycorr_blk-yo,zcorr_blk-zo, sup_avg,v)
            colorbar
            switch ischar(caxis_limits)
                case 1
                otherwise
                    caxis(caxis_limits)
            end
            title('Super-pixel intensity, area, and location')
            xlabel('inches')
            ylabel('inches')
            axis equal
            if output_super_pixels_yz(1)<20 && output_super_pixels_yz(2)<20
                hold on
                for ii=1:output_super_pixels_yz(1)
                    for jj=1:output_super_pixels_yz(2)
                        plot(act_msk_y(ii)-yo,act_msk_z(jj)-zo,'k+')
                    end
                end
            end
        end
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_super_pixel_intensities.eps');
        saveas(gcf,save_fig_name,'epsc')
        
    end
    
    % Fix "zero" superpixels to actual zero in i_sup.  leave ref alone to
    % avoid annoying divide-by-zero errors.
    i_sup(i_sup<0)=0;
    i_sup(ref_sup<0)=0;
    
    %% Ratio to reference values
    i_ratio(:,:,kk)=i_sup./ref_sup;
    
    if strcmp(output_output_eps_figures,'Yes')==1
        
        
        figure
        v=0:max(max(i_ratio(:,:,kk)))/30:max(max(i_ratio(:,:,kk)));
        contourf(act_msk_y-yo, act_msk_z-zo,i_ratio(:,:,kk),v,...
            'LineStyle','none')
        colorbar
        switch ischar(output_caxis_limits)
            case 1
            otherwise
                caxis(output_caxis_limits)
        end
        title('Super-pixel ratio with reference intensities')
        xlabel('inches')
        ylabel('inches')
        axis equal
        if output_super_pixels_yz(1)<30 && output_super_pixels_yz(2)<30
            hold on
            for ii=1:output_super_pixels_yz(1)
                for jj=1:output_super_pixels_yz(2)
                    plot(act_msk_y(ii)-yo,act_msk_z(jj)-zo,'k+')
                end
            end
        end
        save_fig_name=strcat(output_save_dir,...
            output_processing_info_identifier,'_',fig_id_str,...
            '_super_pixel_int_ratio_contour.eps');
        saveas(gcf,save_fig_name,'epsc')
        
    end
    
end

save_ratio_name=strcat(output_save_dir,...
    output_processing_info_identifier,'_ratio.mat');
save(save_ratio_name, 'act_msk_y', 'act_msk_z', 'i_ratio','yo', 'zo',...
    'ref_camera_position','camera_position');

warning('on','all');

if skip_dialog>0
    close all
end

%% Generate Sizing Initialization File

ini_filename=strcat(output_save_dir, output_processing_info_identifier,...
    '_sizing_ini.ini');
% inifile is a function written by by Primoz Cermelj, NOT part of
% matlab.

if exist(ini_filename)>0
    delete(ini_filename)
end

%Make comments section at top of file
inifile(ini_filename,'write',{'','',...
    ';- Lines starting with a semi-colon are treated as a comment.',''})
inifile(ini_filename,'write',{'','',...
    ';- Comment lines do NOT need the = sign at the end :)',''})
inifile(ini_filename,'write',{'','',...
    ';- Filenames should be given with a full pathname',''})
inifile(ini_filename,'write',{'','',...
    ';  if they are not in the current directory.',''})
inifile(ini_filename,'write',{'','',...
    ';- This .ini file may be hand edited; this file was created by',''})
inifile(ini_filename,'write',{'','',...
    ';  the inifile function written by by Primoz Cermelj.  Data is',''})
inifile(ini_filename,'write',{'','',...
    ';  easily read using the same function.',''})

% Input section
input_filename_str=strcat(' '' ',save_ratio_name,''' ');
y_ind_str=strcat('[',num2str(round(size(act_msk_y,2)/2)),']');
z_ind_str=strcat('[1:',num2str(size(act_msk_z,2)),']');
cam_ind_str=strcat('[1:',num2str(size(camera_position,1)),']');
% finds closest angle to 139 degrees for seed
[c sd_ind]=min(abs(139-cam_angle));
seed_ind_str=strcat('[',num2str(sd_ind),']');
scatter_calc_str=strcat(' '' ',test_section_scatter_calc_meth,''' ');

INPUT_cell_mat={'','Input', 'ratio_file',input_filename_str;...
    '','Input', 'y_indexes',y_ind_str;...
    '','Input', 'z_indexes',z_ind_str;...
    '','Input', 'camera_indexes',cam_ind_str;...
    '','Input', 'seed_index',seed_ind_str;...
    '','Input', 'wavelength',num2str(test_section_wavelength);...
    '','Input', 'xi',num2str(test_section_xi);...
    '','Input', 'gamma_ref',num2str(test_section_gamma_ref);...
    '','Input', 'scatter_calc_meth',scatter_calc_str;...
    '','Input', 'half_angle',num2str(test_section_half_angle)};
inifile(ini_filename,'write',INPUT_cell_mat)

% Output section
output_name_str=strcat(' '' ',output_save_dir,'sizing/output.mat',''' ');
OUT_cell_mat={'','Output', 'output_filename',output_name_str;...
    '','Output', 'make_a_plot','0';...
    '','Output', 'plot_angle_range','[135:0.5:150]'};
inifile(ini_filename,'write',OUT_cell_mat)

% Config section
scat_db_filename=strcat(' '' ',...
'../../database_files/scattering_coefficients_water_sub_angs.mat',...
   ''' ');
int_vs_size_db_str=strcat(' '' ',...
    '../../database_files/sz3_data_file_water.mat'' ');
Config_cell_mat={'','Config', 'int_vs_size_filename',...
                 int_vs_size_db_str;...
    '','Config', 'cross_section_database_filename',scat_db_filename;...
    '','Config', 'try_sigma',...
                 '[6.5 10 20 40 2.25 2.5 2.75 3.1 3.5 4 4.8 6]';...
    '','Config', 'use_poly',' ''No'' ';...
    '','Config', 'poly_order','5';...
    '','Config', 'try_dist_type',['{''logn'',''logn'',''logn'','...
                 '''logn'',''norm'',''norm'',''norm'',''norm'','...
                 '''norm'',''norm'',''norm'',''norm''}']};
inifile(ini_filename,'write',Config_cell_mat)


%% Close file(s)
fclose('all');

% creates a compressed archive of the output directory.  Over 20x savings
% typical of this data
zip_string=['nice -n 10 tar -jcf ' ...
    output_save_dir(1:(length(output_save_dir)-1)) '.tar.bz2 ' ...
    output_save_dir ' &'];
[status, result]=system(zip_string);

