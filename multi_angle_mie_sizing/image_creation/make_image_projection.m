function [corrected_img] = make_image_projection(target_coords,...
    camera_coords,image_width,y_pix, z_pix, orig_img)
% returns a perspective image given camera and target coordinates
%
% target_coords: [x y z] point that the camera is pointed toward.
% camera_coords: [x y z] point where camera is located.
% image_width: dimensional width of the flat image
% y_pix, z_pix: desired size of the corrected image
% orig_img: the "flat" image
% corrected_image: returned image with perspective transform applied
%
% Copyright (C) 2013, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% Need path to plane_line_intersect.m, a third-party function
addpath ../../multi_angle_mie_sizing_3rd_party/helper_functions

%camera_coords=[5 -20 0];
%target_coords=[0 0 0];

% normal vector to the camera image plane
dd=camera_coords-target_coords;
dd=dd./norm(dd);

%image_width=0.15;
%y_pix=225;
%z_pix=900;

% resolution, pixels/inch
pix_res=y_pix/image_width;
image_height=z_pix/pix_res;

% last "corner" is a point on the z-axis - effectively assuming that the
% camera orientation will always be vertical
% image_corners=[0  -image_width/2 -image_height/2
%                0   image_width/2 -image_height/2
%                0   image_width/2  image_height/2
%                0  -image_width/2  image_height/2
%                0   0              image_height/2];
         
image_corners=[0   0              0
               0   image_width    0
               0   image_width    image_height
               0   0              image_height
               0   image_width/2  image_height];
%

% Center the image about zero
image_corners(:,2)=image_corners(:,2)-target_coords(2);
image_corners(:,3)=image_corners(:,3)-target_coords(3);
% So now:
target_coords=[0 0 0];

for ii=1:size(image_corners,1)
    proj_pts(ii,:)=...
        plane_line_intersect(dd,target_coords,camera_coords,...
        image_corners(ii,:));
end

qd=[-1 -1
     1 -1
     1  1
    -1  1];

dist5=sqrt(sum((proj_pts(5,:)-target_coords).^2));
for ii=1:size(image_corners,1)-1
    dist_ii=sqrt(sum((proj_pts(ii,:)-target_coords).^2));
    dist_5ii=sqrt(sum((proj_pts(5,:)-proj_pts(ii,:)).^2));
    ang5=acos((dist5^2+dist_5ii^2-dist_ii^2)/(2*dist5*dist_5ii));
    dx=qd(ii,1)*dist_5ii*sin(ang5);
    dy=-dist_5ii*cos(ang5)+dist5;
    proj_corners(ii,:)=[dx dy];
end

% Shows the flat and the transformed image corners that will be used
% in the transform
% figure
% plot(image_corners(1:4,2),image_corners(1:4,3),'+')
% hold on
% plot(proj_corners(:,1),proj_corners(:,2),'o')

proj_corners(:,1)=proj_corners(:,1)+abs(min(image_corners(1:4,2)));
proj_corners(:,2)=proj_corners(:,2)+abs(min(image_corners(1:4,3)));

image_corners(1:4,2)=image_corners(1:4,2)+abs(min(image_corners(1:4,2)));
image_corners(1:4,3)=image_corners(1:4,3)+abs(min(image_corners(1:4,3)));

% AGain,shows the flat and the transformed image corners that will be used
% in the transform
% figure
% plot(image_corners(1:4,2),image_corners(1:4,3),'+')
% hold on
% plot(proj_corners(:,1),proj_corners(:,2),'o')

mytform = cp2tform(pix_res*image_corners(1:4,2:3),...
    pix_res*proj_corners, 'projective');

[corrected_img xdata ydata]=imtransform(orig_img, mytform,...
    'Xdata',[0 y_pix],'Ydata',[0 z_pix],'Size', [z_pix y_pix]);

% final projection correction
% figure
% imshow(corrected_img)

