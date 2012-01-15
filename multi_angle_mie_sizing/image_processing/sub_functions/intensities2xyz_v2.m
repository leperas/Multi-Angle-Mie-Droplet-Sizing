function [i_exact i_blk i_sup sup_cor_y sup_cor_z img_blk sup_avg]=...
    intensities2xyz_v2(flat_image, img_distance, img_exposure_time,...
    pix_res, z_pix, y_inds, z_inds, ycorr, zcorr,...
    output_super_pixels_yz, output_super_pixels_size,...
    output_super_pixel_shape, fids)
% Summary of output: 
% EXACT data: In the full-resolution immages, ycorr and zcorr
% contain the exact coordinates for the center of each pixel.  Intensities
% are in i_exact which is of the size y_pix and z_pix.
% The sub-block of data we want to anylize is selected from within the 
% above full resolution data.  ycorr_blk and zcorr_blk contain exact 
% coordinates for the intensities in i_blk, which is of size z_pix_blk and
% y_pix_blk.
% For example, the intensity at zcorr_blk(7), ycorr_blk(9) is
% i_blk(7,9)
% SUPERPIXEL data: sup_cor_y and sup_cor_z contain exact coordinates for
% the intensities in i_sup, which is of size output_super_pixels_yz.   
% For example, the intensity of the superpixel located at
% sup_cor_z(3) and sup_cor_y(5) is i_sup(3,5).
% IMAGES: img_blk is simply the sub-block image itself.
% sup_avg is a matrix of the same resolution as img_blk and i_blk, but 
% contains the values of/within the superpixels (values outside a 
% superpixel are zero), thus showing both the value and location of the 
% superpixels within the selected block.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% any superpixel that comes out to be "zero" will be assigned this value
% helps plotting zeros appear more clearly.
superpixel_zero_value=-0.05;  

% allocate the flat image
img_blk=flat_image(z_inds,y_inds);

% Reference intensity matrix, corrected for exposure
i_exact=flat_image/(img_exposure_time);
i_blk=img_blk/(img_exposure_time);

% interp block onto coarse grid/super-pixel averaging.

size_info_blk=size(i_blk);
z_pix_blk=size_info_blk(1);
y_pix_blk=size_info_blk(2);

if output_super_pixels_yz(1)>y_pix_blk
    mfprintf(fids,['Number of Y-axis super pixels greater than image'...
        ' resolution.\nAdjusted from %i to %i \n'],...
        output_super_pixels_yz(1),y_pix_blk);
    output_super_pixels_yz(1)=y_pix_blk;
end
if output_super_pixels_yz(2)>z_pix_blk
    mfprintf(fids,['Number of Z-axis super pixels greater than image'...
        ' resolution.\nAdjusted from %i to %i \n'],...
        output_super_pixels_yz(2),z_pix_blk);
    output_super_pixels_yz(2)=z_pix_blk;
end

% geometric center of super-pixel
%actual edges of selected area (vs requested left edge from .ini file)
lft=ycorr(y_inds(1))-0.5/pix_res;
rht=ycorr(y_inds(y_pix_blk))+0.5/pix_res;
tp=zcorr(z_inds(1))+0.5/pix_res;
bt=zcorr(z_inds(z_pix_blk))-0.5/pix_res;

mfprintf(fids,['Actual data block limits, Left, Right, Top, '...
    'Bottom: [%f %f %f %f] \n'],[lft rht tp bt]);
mfprintf(fids,'Block height: %f  width: %f \n',tp-bt,rht-lft);

% Find center of super pixels
for ii=1:output_super_pixels_yz(1)
    yc_s(ii)=((rht-lft)/output_super_pixels_yz(1))*(ii-1)+...
        0.5*((rht-lft)/output_super_pixels_yz(1))+lft;
end
for ii=1:output_super_pixels_yz(2)
    zc_s(ii)=((tp-bt)/output_super_pixels_yz(2))*...
        (output_super_pixels_yz(2)-ii)+0.5*((tp-bt)/...
        output_super_pixels_yz(2))+bt;
end

i_sup=zeros(output_super_pixels_yz(2),output_super_pixels_yz(1));

% protect against over-lapping superpixels - this is not allowed.
max_y_sup=floor(pix_res*(rht-lft)/output_super_pixels_yz(1))-1;
if max_y_sup<1
    max_y_sup=1;
end
if floor(output_super_pixels_size*pix_res)>max_y_sup
    mfprintf(fids,'Super pixels size adjusted from %6.4f to %6.4f \n',...
        output_super_pixels_size,max_y_sup/pix_res);
    output_super_pixels_size=max_y_sup/pix_res;
end
max_z_sup=floor(pix_res*(tp-bt)/output_super_pixels_yz(2))-1;
if max_z_sup<1
    max_z_sup=1;
end
if floor(output_super_pixels_size*pix_res)>max_z_sup
    mfprintf(fids,'Super pixels size adjusted from %6.4f to %6.4f \n',...
        output_super_pixels_size,max_z_sup/pix_res);
    output_super_pixels_size=max_z_sup/pix_res;
end

yc_s_pix=round((yc_s-lft)*pix_res+0.5);  % The 0.5 must be added due to 
                                         % using pixel CENTERS...
zc_s_pix=round((zc_s-bt)*pix_res+0.5);

% Make super-pixel mask
switch output_super_pixel_shape
    case 'square'
        square_size=floor(output_super_pixels_size*pix_res);
        mask=strel('square',square_size);
        mfprintf(fids,['Actual super pixel square side length: %e '...
            'inches.\n'],square_size/pix_res);
    otherwise
        disc_radius=floor(output_super_pixels_size*pix_res/2);
        mask=strel('disk',disc_radius,0);
        mfprintf(fids,['Actual super pixel circle diameter: %e '...
            'inches.\n'],2*(disc_radius+0.5)/pix_res);
end

msize=size(getnhood(mask));

num_pix_averaged=sum(sum(getnhood(mask)));
super_pixel_area=num_pix_averaged*(1/pix_res)^2;
mfprintf(fids,'Single pixel area: %e square inches.\n',(1/pix_res)^2);
mfprintf(fids,'Number of actual pixels averaged in a superpixel: %i\n',...
    num_pix_averaged);
mfprintf(fids,'Total super pixel area: %e square inches.\n',...
    super_pixel_area);

% Will be able to show the averaged superpixels
sup_avg=zeros(z_pix_blk,y_pix_blk);

for ii=1:output_super_pixels_yz(1)
    for jj=1:output_super_pixels_yz(2)
        mat_mask=zeros(z_pix_blk,y_pix_blk);
        
        d_col=yc_s_pix(ii)-floor(msize(2)/2);
        d_row=zc_s_pix(jj)-floor(msize(1)/2);
        
        mat_mask(d_row:d_row+msize(1)-1,d_col:d_col+msize(2)-1)=...
            mat_mask(d_row:d_row+msize(1)-1,d_col:d_col+msize(2)-1)+...
            getnhood(mask);
        
        sup_cor_y(ii)=(d_col-1)/pix_res+0.5*msize(2)/pix_res+lft;
        sup_cor_z(jj)=(z_pix_blk-d_row+1)/pix_res-0.5*msize(1)/pix_res+bt;
      
        sup_region=i_blk.*mat_mask;
        if max(max(sup_region))==0
            i_sup(jj,ii)=superpixel_zero_value; % gives enough contrast 
                                    % to zero point for contours to appear
        else
            i_sup(jj,ii)=mean(nonzeros(sup_region));
        end
        sup_avg=sup_avg+i_sup(jj,ii).*mat_mask;
    end
end
