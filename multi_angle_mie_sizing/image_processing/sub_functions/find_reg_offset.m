function correct_coords=find_reg_offset(coords, sc_crd, camera_position)
% coordinates are in pixels, camera position inches [x y z]
% coords will be corrected for planar offset between data plane and 
% registration plane, and then returned again in pixels.
% planar offset is in the x axis, negative is for registration points
% "behind" the data plane, as is most likely.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

x1=camera_position(1);
y1=camera_position(2);
z1=camera_position(3);

delta_pix=coords(2,:)-coords(1,:);
delta_y=sc_crd(2,:)-sc_crd(1,:);
pix_res=norm(delta_pix)/norm(delta_y);

correct_coords=zeros(size(coords,1),2);

pln_crd=sc_crd;

for ii=1:size(coords,1)
    x2=pln_crd(ii,1);
    y2=pln_crd(ii,2);
    z2=pln_crd(ii,3);
    y=y1-x1*(y2-y1)/(x2-x1);
    z=z1-x1*(z2-z1)/(x2-x1);
    proj=[0 y z];
    %offset, in pixels
    offset=(sc_crd(ii,:)-proj)*pix_res;
    offset=offset.*[1  1 -1];
    correct_coords(ii,:)=coords(ii,:)+offset(2:3);
end


