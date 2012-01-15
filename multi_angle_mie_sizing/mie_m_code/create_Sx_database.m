clear all
% creates database of scattering coefficients
%
% This may take a long time - status of the calculation is periodically
% stored in status_scat_calc.mat as the iteration, jj and the time up to
% that point.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% location of the mie function mie_S12.m and required subroutines
addpath ../../multi_angle_mie_sizing_3rd_party/mie_functions

% Desired database name
scattering_coefficients_database_name='scattering_coeff_database.mat';

% specify dimensional parameters for the desired database
max_d=200;  % maximum geometric diameter to calculate, micrometers

wavelength=514.5; % nanonmeters

% size resolution in the database
d_inc=0.05;  % micrometers

% measurement angle, radiens
theta=(pi/180)*[38:0.1:42 130:0.1:160]; % range 

m=1.333+0i; % complex index refraction for particle

% No need to edit below here...
tic

points=max_d/d_inc;
particle_size=d_inc:d_inc:max_d;

% size parameter will be saved in the database
x=2*pi*(particle_size/2)*1e3/wavelength;

S1=zeros(size(particle_size,2),1)';
S2=zeros(size(particle_size,2),1)';
V1 = struct('S1',S1, 'S2', S2);
A = repmat(V1,1,size(theta,2));

B.m=m;
B.x=x;
B.theta=theta;

u=cos(theta);

result=zeros(points,2);
tm=zeros(size(theta,2),1)';

for jj=1:size(theta,2);
    theta(jj)*180/pi;
    for ii=1:points
        result(ii,:)=mie_S12(m,x(ii),u(jj));
    end
    A(jj).S1=result(:,1);
    A(jj).S2=result(:,2);
    calc_time(jj)=toc;
    save(scattering_coefficients_database_name, 'A', 'B');
    save status_scat_calc jj calc_time
end


