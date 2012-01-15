% select a set of images which have been created
% by placing increasingly absorbing neutral density filters in front of the
% camera.  The set of ND filters used is input into the processing code,
% and images are selected in correlation to the input filters.  Or use
% shutter times.  Pick which by setting FOS variable.  The size
% and shape of the area to be processed is input, and any color channel and
% thresholding desired is configured.  A data file is created during the
% run which contains the ND filter information and intensity information
% from the calibration images
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

clear all
close all

current_path=pwd; 
startpath=current_path;  % sets the startpath in the dialog to current path
% or manually define the start path if you want

% Filters=0 or shutter times=1
FOS=1;

ND_filters=[0.0 0.1 0.2 0.3 0.5 1.0];  % value of the ND filters, in the 
%order you will select the associated images

SHT_times=[0.4 0.5 0.6 1.6 2 2.5];  % value of shutter times, seconds, 
% in the order you will select the associated images

crop_w=200;  % sets the size of the selection area, width by height
crop_h=200;
% can set thresholds prior to processing.  seemed like a good idea at first
% but I wound up never using it.
threshold_v=1;
max_v=255;
use_color=2;

% ---------------------------------------------------------------------- 
% No "user" configuration below.  Do not adjust things below here without 
% thinking it through.

if FOS==0
    transmit=10.^-ND_filters;
    
    % outputs the filters and transmittance before calling the dialog, just
    % to help in selection of the correct images.
    [ND_filters' transmit']
    
    [data_files, pathname, filterindex] = uigetfile('*.*', ['Select' ...
        ' image files corresponding to displayed ND filter strengths:'],...
        startpath ,'MultiSelect', 'On');
else
    transmit=SHT_times/SHT_times(1);
    
    % outputs the filters and transmittance before calling the dialog, just
    % to help in selection of the correct images.
    [SHT_times' transmit']
    
    [data_files, pathname, filterindex] = uigetfile('*.*', ['Select'...
        ' image files corresponding to displayed shutter times:'],...
        startpath ,'MultiSelect', 'On');
end

addpath(pathname);

if ischar(data_files)==1
    max_loop=1;
else
    max_loop=size(data_files,2);
end
for ii=1:max_loop
    figure(ii)
    if ischar(data_files)==1
        I=imread(data_files);
    else
        I=imread(cell2mat(data_files(ii)));
    end
    
    %I=rgb2gray(I);   % gray=bad, don't use this.
    
    warning('off','Images:initSize:adjustingMag');
    imshow(I);
    warning('on','all');
    % pick center of desired region, image of size defined above will be
    % stored in I2
    h=impoint;
    rect_pos = getPosition(h);
     cent_x=rect_pos(1)-crop_w/2;
     cent_y=rect_pos(2)-crop_h/2;
%    cent_x=rect_pos(1); cent_y=rect_pos(2);
    rect=[cent_x cent_y crop_w crop_h];
    I2=imcrop(I,rect);
    imshow(I2);
    % use ROI to create a mask of only G values between threshold and max
    BWg = roicolor(I2(:,:,use_color),threshold_v,max_v);
    I3=immultiply(I2(:,:,use_color), BWg);
    max_i(ii)=max(max(I3));
    imshow(I3);
    [output,xx]=imhist(I3);
    avg_i(ii)=sum(sum(I3))/(sum(output)-output(1));
 
end

if FOS==0
    CCD_intensity_ratio=avg_i/avg_i(1);
    CCD_percent=avg_i/255;
    
    GG=[ND_filters' transmit' CCD_intensity_ratio' CCD_percent' avg_i']
else
    CCD_intensity_ratio=avg_i/avg_i(1);
    CCD_percent=avg_i/255;
    
    GG=[SHT_times' transmit' CCD_intensity_ratio' CCD_percent' avg_i']
end
    
    
[save_file, save_pathname, filterindex] = uiputfile('*.*', ['Filename'...
    ' to save calibration data:'], strcat(pathname,'CCD_calibration.txt'));

save_string=strcat(save_pathname,save_file);

save(save_string, 'GG', '-ASCII', '-TABS')

close all
