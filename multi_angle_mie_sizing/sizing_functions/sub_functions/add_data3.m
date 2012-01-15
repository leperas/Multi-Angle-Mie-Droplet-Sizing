function result=add_data3(theta_in,phi_in,sigma,dist_type,ini_filename,...
    sz3_data, d_mean)
% If called to create a "virgin" size database look-up-table, input the
% below conditions for the 1st entry in the database.
% theta_in: scattering angle, radians
% phi_in: angle between scattering plane and reference plane, radians
% sigma: should be in terms of the size parameter, x, as so:
%       sigma=2*pi*(sigma_dimensional)*1e3/wavelength
%       where wavelength is in nanometers
% dist_type: 'logn' or 'norm'
% ini_filname: the filename of the sizing initialization file containing
%       the relavant optical setup parameters
% sz3_data: set to zero (0) if virgin, else this is the database structure
% d_mean: vector of size parameters, definies size resolution of the 
%       database entry; size parameters will be interpolated between 
%       adjacent database values.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% where A and B are from the cross section database need to be global
global A B

[keys,sections,subsections] = inifile(ini_filename,'readall');

key_assignment_string=strcat(keys(:,2),'_',keys(:,3),'=',keys(:,4),';');

% Assign values to the initalization file variables
for ii=1:size(key_assignment_string,1)
    eval(char(key_assignment_string(ii)))
end

% path to the Irr_int.m function is required
addpath ../../mie_m_code

fids=1;  % 1 for output to screen, add additional fids as desired

[c ind]=min(abs(B.theta-theta_in));

% create an initial structure if sz3_data=0, i.e. if starting from scratch
 if isstruct(sz3_data)==0
     clear sz3_data
     mfprintf(fids,'Creating an empty record in a new structure.\n');
     intensity=zeros(size(d_mean));
     V1 = struct('theta', B.theta(ind), 'phi', phi_in, 'half_angle', ...
         input_half_angle, 'gamma_ref', num2str(input_gamma_ref), 'xi',...
         num2str(input_xi), 'method',input_scatter_calc_meth, 'sigma', ...
         sigma, 'd_mean', d_mean, 'intensity', intensity, ...
         'dist_type',dist_type);
     sz3_data = repmat(V1,1,1);
 end

% See if data record exists.
[chk_data rec_ind]=chk4data3(B.theta(ind),phi_in,...
    input_half_angle,num2str(input_gamma_ref),num2str(input_xi),...
    input_scatter_calc_meth, sigma,dist_type,sz3_data);

if isstruct(chk_data)==0
    %add a new data record
    mfprintf(fids,'Adding new data record.\n');
    st_ind=1+size(sz3_data,2);
else
    %overwrite the existing data record
    mfprintf(fids,'Over-writing data record.\n');
    st_ind=rec_ind;
end

% Calculate the size vs. intensity data
mfprintf(fids,'Working on: %5.1f,  %5.1f, Completion %%:\n0.0',...
    180*theta_in/pi(),180*phi_in/pi());

cnt=0;
for kk=1:size(d_mean,2)
    cnt=cnt+1;
    if cnt>(size(d_mean,2)/10)
        mfprintf(fids,'...%2.0f',100*kk/size(d_mean,2));
        cnt=0;
    end
    
    % Must convert to micrometers as Irr_int.m uses that as input.  Might
    % adjust this in the future to accept x.  Database is in x anyway.
    %rl_x=2*input_wavelength*d_mean/(2*pi*1e3);    
    [II Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
                Irr_int(d_mean(kk), sigma, dist_type, B.theta(ind),...
                phi_in, input_half_angle, input_gamma_ref,input_xi,...
                input_scatter_calc_meth);
    R_II(kk)=II;
    
end
mfprintf(fids,'...Done!\n');

% assign values actually used in calc, not the requested angles
sz3_data(st_ind).theta=B.theta(ind);
sz3_data(st_ind).phi=phi_in;
sz3_data(st_ind).half_angle=input_half_angle;
sz3_data(st_ind).gamma_ref=num2str(input_gamma_ref);
sz3_data(st_ind).xi=num2str(input_xi);
sz3_data(st_ind).method=input_scatter_calc_meth;
sz3_data(st_ind).dist_type=dist_type;
sz3_data(st_ind).sigma=sigma;
sz3_data(st_ind).d_mean=d_mean;
sz3_data(st_ind).intensity=R_II;


% return the sturcture with the new data added/written
result=sz3_data;