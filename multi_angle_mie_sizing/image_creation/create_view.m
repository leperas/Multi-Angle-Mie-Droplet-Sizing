function [view_img H]=create_view(target_coords,camera_coords,...
    view_angle,image_width,orig_img)
% returns a perspective image given camera and target coordinates
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>
show_figs=0;

rot_img=imrotate(orig_img, -90);

d=size(rot_img);

pix_res=d(1)/image_width;
image_height=d(2)/pix_res;

% image translated so center of data distribution is at 0,0,0 origin.
zdata=[0 image_height; 0 image_height]-target_coords(3);
ydata=[0 0; image_width image_width]-target_coords(2);

% this way origin is at bottom left, not used
%zdata=[0 image_height; 0 image_height];
%ydata=[0 0; image_width image_width];

xdata=[0 0; 0 0];
cdata=rot_img;

H=figure;

surface('XData',xdata,'YData',ydata,...
        'ZData',zdata,'CData',cdata,...
        'FaceColor','texturemap','EdgeColor','none');

camproj('perspective')

axis equal
axis image
axis off

camva(view_angle);

camtarget([0 0 0]);

campos(camera_coords);
    
jj=getframe(gcf);

view_img=jj.cdata;

if show_figs==1
    figure
    title('Original image')
    imshow(orig_img)
    
    figure
    title('Viewed image')
    imshow(view_img)
end
