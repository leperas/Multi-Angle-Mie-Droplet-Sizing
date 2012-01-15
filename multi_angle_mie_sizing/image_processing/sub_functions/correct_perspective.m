function [corrected_img xdata ydata]=correct_perspective(data_img, ...
    flat_coords, registration_coords, y_pix, z_pix)
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% transform created back to the 'flat' image
mytform = cp2tform(registration_coords, flat_coords, 'projective');

[corrected_img xdata ydata]=imtransform(data_img, mytform,...
    'Xdata',[0 y_pix],'Ydata',[0, z_pix],'Size', [z_pix y_pix]);
