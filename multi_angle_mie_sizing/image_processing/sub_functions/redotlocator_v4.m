function XY_coords = redotlocator_v4(img, show_figs) 
% if show_figs=1 figs will show, else they won't
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

warning('off','Images:initSize:adjustingMag');

max_dots=4;
cut_off=0.4;
d_size=10;

l=size(img);
img(img==0)=1; % makes math happier

RGB=double(img)/255;
R=RGB(:,:,1);
G=RGB(:,:,2);
B=RGB(:,:,3);

R_rel=R./(R+G+B);

R_rel_im=uint8(round(255*R_rel));

R_rel(R_rel<=cut_off)=0;
R_rel(R_rel>cut_off)=1;


img_both=uint8(255*R_rel);


warning('off', 'MATLAB:intConvertOverflow')
warning('off', 'MATLAB:intMathOverflow')

% gets rid of specs
filtered=medfilt2(img_both,[10 10]);
% fills in holes
background = imclose(filtered,strel('disk',d_size));

II=imadd(background,rgb2gray(img));


% Find the boundaries Concentrate only on the exterior boundaries.
% Option 'noholes' will accelerate the processing by preventing
% bwboundaries from searching for inner contours. 
[B,L] = bwboundaries(background, 'noholes');

% Determine objects properties

STATS = regionprops(L, 'all'); % we need Centroid'

if show_figs==1
    close all
    figure
    imshow(img)
    figure
    imshow(R_rel_im)
    figure
    imshow(img_both)
    figure
    imshow(filtered)
    figure
    imshow(background)
    figure
    imshow(II)
    figure
    imshow(img)
    hold on
end



Y_count = 0;
Z_count = 0;

for i=1:length(STATS)
   %DD(i)=size(B{i},1);
   DD(i)=STATS(i).Area;
end

JJ=sortrows(DD',-1);

for i = 1 : length(STATS)
   
    %if size(B{i},1)>=JJ(max_dots)
    if STATS(i).Area>=JJ(max_dots)
        
        centroid = STATS(i).Centroid;
        
        Y_coor = centroid(1);
        Y_count = Y_count+1;
        Z_coor = centroid(2);
        Z_count = Z_count+1;
        Y_inputs(Y_count) = Y_coor;
        Z_inputs(Z_count) = Z_coor;

    end
end

XY_coords=[Y_inputs' Z_inputs'];

if size(XY_coords,1)>max_dots
    err_strg=strcat('More than: ',num2str(max_dots),' red dots found!!!');
    error(err_strg)
end

% Sort coords so that "1" is bottom lefthand corner, then go CCW

XY_temp=sortrows(XY_coords,-2);

XY_sorted=sortrows(XY_temp(1:2,:),1);

XY_sorted2=sortrows(XY_temp(3:4,:),-1);

XY_coords=[XY_sorted;XY_sorted2];

for i=1:4        
        if show_figs==1
            plot(XY_coords(i,1),XY_coords(i,2),'w.');
            H=text(XY_coords(i,1)-25,XY_coords(i,2)-25,num2str(i));
            set(H,'Color',[1 1 1]);
        end
 end



warning('on','all');
return
