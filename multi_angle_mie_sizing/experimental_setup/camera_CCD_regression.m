% This code applies the desired regression to data generated in the
% camera_CCD_calibration code.  Multiple calibration files may be selected.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

clear all
close all

f_order=5;   % order of polynomial
cutoff=0.28;  % cutoff between high brightness polynomial and linear

current_path=pwd;
startpath=current_path;  % sets the startpath in the dialog to current path
% or manually define the start path if you want

% Probably no need to edit below here very often, unless you need to add
% "extra points" to the fit

[data_files, pathname, filterindex] = uigetfile('*.*', ['Select CCD'...
    ' calibration file:'], startpath ,'MultiSelect', 'On');

addpath(pathname);

figure(1)
hold on

if ischar(data_files)==1
    max_loop=1;
else
    max_loop=size(data_files,2);
end
for ii=1:max_loop
    if ischar(data_files)==1
        S=load(data_files);
    else
        S=load(cell2mat(data_files(ii)));
    end
    
    [ssi sj]=size(S);
    % add number to size to make room for "zero" point and other extra
    % points.  i.e. for one extra point, num_x=1.
    num_x=0;
    si=ssi+num_x;
    im=ii*si-si;
    ND_filters(im+1:im+si-num_x)=S(:,1);
    transmit(im+1:im+si-num_x)=S(:,2);
    CCD_intensity_ratio(im+1:im+si-num_x)=S(:,3);
    CCD_percent(im+1:im+si-num_x)=S(:,4);
    avg_i(im+1:im+si-num_x)=S(:,5);
    % add extra points to fit
    % zero point
%     CCD_percent(im+si)=0;
%     transmit(im+si)=0;
% other point
%%    CCD_percent(im+si)=CCD_percent(im+si-num_x)/1.2;
%%    transmit(im+si)=transmit(im+si-num_x)/1.2;
    % add three more extra points
%     CCD_percent(im+si-1)=CCD_percent(im+si-num_x)/2;
%     transmit(im+si-1)=transmit(im+si-num_x)/2;
%     CCD_percent(im+si-2)=CCD_percent(im+si-num_x)/4;
%     transmit(im+si-2)=transmit(im+si-num_x)/4;
%     CCD_percent(im+si-3)=CCD_percent(im+si-num_x)/8;
%     transmit(im+si-3)=transmit(im+si-num_x)/8;
    
    R=polyfit(CCD_percent(im+1:im+si),transmit(im+1:im+si),f_order);
    xrange=min(CCD_percent(im+1:im+si)):0.01:1;
    plot(CCD_percent(im+1:im+si), transmit(im+1:im+si),'o',xrange,...
        polyval(R,xrange),'-')
    P=R/polyval(R,1.0);
end

R=polyfit(CCD_percent,transmit,f_order);

xrange=0:0.01:1;
plot(xrange,polyval(R,xrange),'-r')
xlabel('CCD Percent')
ylabel('Transmittance')

% This is the high-order fit
P=R/polyval(R,1.0);

% this is the linear fit up to the cutoff CCD_percent point
yval_cut=polyval(P,cutoff);
L=polyfit([0 cutoff],[0, yval_cut],1);
xrange1=0:0.01:cutoff;
xrange2=cutoff:0.01:1;

figure(2)
plot(xrange1, polyval(L,xrange1),'-r',xrange2,polyval(P,xrange2),'-r')
ylim([0,1]);
xlabel('CCD Percent')
ylabel('Corrected Transmittance')

[save_file, save_pathname, filterindex] = uiputfile('*.mat',...
    'Filename to save polynomial:', strcat(pathname,'CCD_poly.mat'));

save_string=strcat(save_pathname,save_file);

save(save_string, 'P', 'L', 'cutoff');
