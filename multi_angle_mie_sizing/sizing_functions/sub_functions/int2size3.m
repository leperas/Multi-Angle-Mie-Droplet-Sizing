function [mean_size]=int2size3(varargin)
%  input parameters are:
%    (data_filename, int_vs_size_filename, cross_section_database_filename,
%                     half_angle, sigma, dist_type)
%
% if data_filename is a number, set
% skip_dialog=0, ask for files, print output
% skip_dialog=1, ask for files, don't print output
% else
% skip_dialog=2, don't ask, don't print output
%
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% A and B must already be global and loaded from the database
global A B

data_filename=varargin{1};

if isnumeric(data_filename)==1
    skip_dialog=data_filename;
    current_path=pwd;   % gets current path for use in dialogs
    [data_file, current_path_d, filterindex] = uigetfile...
        ('*.mat', 'Select intensity ratio data file:',...
        current_path ,'MultiSelect', 'Off');
    data_filename=strcat(current_path_d,data_file);
    pause(1);
    [size_file, current_path_s, filterindex] = uigetfile...
        ('*.mat', 'Select intensity vs. size data file:',...
        current_path ,'MultiSelect', 'Off');
    int_vs_size_filename=strcat(current_path_s,size_file);
    pause(1);
    [cross_file, current_path_c, filterindex] = uigetfile...
        ('*.mat', 'Select intensity cross section database file:', ...
        current_path ,'MultiSelect', 'Off');
    cross_section_database_filename=strcat(current_path_c,cross_file);
    
    prompt = {'Enter camera lens half-angle (degrees):',...
        'Enter size distribution standard deviation:'};
    dlg_title = 'Input for calcuation:';
    num_lines = 1;
    def = {'0.4','10'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    half_angle=str2num(answer{1});
    sigma=str2num(answer{2});
    
else
    skip_dialog=2;
    
    output_output_filename=varargin{2};
    int_vs_size_filename=varargin{3};
    half_angle=varargin{4};
    gamma_ref=varargin{5};
    xi=varargin{6};
    method=varargin{7};
    sigma=varargin{8};
    dist_type=varargin{9};
    ini_filename=varargin{10};
    
end

% create output filename and logfile name 
result=strfind(output_output_filename,'.');
if numel(result)==0
    output_filename=strcat(output_output_filename,'_size_sig',...
        num2str(sigma),'_',dist_type,'.mat');
    log_filename=strcat(output_output_filename,'_size_sig',...
        num2str(sigma),'_',dist_type,'.log');
else
    lst_dt=result(size(result,2));
    output_filename=strcat(output_output_filename(1:lst_dt-1),...
        '_size_sig',num2str(sigma),'_',dist_type,...
        output_output_filename(lst_dt:length(output_output_filename)));
    log_filename=strcat(output_output_filename(1:lst_dt-1),'_size_sig',...
        num2str(sigma),'_',dist_type,'.log');
end
[out_fid message]=fopen(log_filename,'w');
if out_fid<0
    error(strcat('Failure to create log file: ',message))
end

if skip_dialog>0
    fids=[out_fid];
else
    fids=[1 out_fid];
end

% load the intensity ratio data
load(data_filename,'-mat');

% A B must be global and available
% This determines the angular and size resolution available for the
% calcualtion.

% load the structure sz_data, which contains all previously made intensity
% ratio vs. size calcualtions.  It is the user's responsibilty to pick a
% file which contains correct data for the given situation. (i.e correct
% index of refraction, wavelength etc.
use_lockfile('sz_data_lock_file.txt','lock');
load(int_vs_size_filename);
use_lockfile('sz_data_lock_file.txt','unlock');

% Find starting x point, start 'norm' distributions 3 sigma above zero.
switch dist_type
    case 'logn'
        ind_Bmin=1;
    case 'norm'
        B_min=(3/2)*sigma;
        [c ind_Bmin]=min(abs(B.x-B_min));
    otherwise
        error('Unrecognized distribution')
end

% how course the database spacing will be in x
d_mean_skip=4;
d_mean=B.x(ind_Bmin:d_mean_skip:end);
% d_mean seach index; one less than the total size of d_mean, makes search
% function simple
d_si=1:(size(d_mean,2)-1);

% Courseness in phi of the database (radians)
phi_inc=1*pi()/180;

tolp=phi_inc/2;
n_phis=round((pi()/2)/phi_inc);
phi=zeros(1,2*n_phis+1,1);
phi(n_phis+1)=0;
phi(1:n_phis)=phi_inc*(n_phis+1-(1:n_phis));
phi(n_phis+2:2*n_phis+1)=phi_inc*(1:n_phis);

% for each value in i_ratio, there first must be an entry in the local
% intensity vs size database

ny=size(act_msk_y,2);
nz=size(act_msk_z,2);
n_angs=size(i_ratio,3);

r_dist=zeros(nz,ny);
d_dist=zeros(nz,ny);
r_ang=zeros(nz,ny);
d_ang=zeros(nz,ny);

% makes i_ratio_c a 3-D array, even if data in i_ratio was only 2-D
i_ratio_c=zeros(nz,ny,n_angs);
if size(i_ratio,3)==1
    i_ratio_c(:,:,1)=i_ratio;
else
    i_ratio_c=i_ratio;
end

% empty local data struc, will be populated from sz_data file database as
% records are needed
l_sz_data = struct('theta',0,'phi',0,'half_angle', 0, 'gamma_ref','0',...
    'xi', '0', 'method','0', 'sigma', 0,...
    'd_mean', 0, 'intensity', 0,'dist_type','logn');
l_sz_data = repmat(l_sz_data,1,1);
l_sz_ind=0;

for kk=1:n_angs
    for jj=1:nz
        for ii=1:ny
            point=[0, act_msk_y(ii) act_msk_z(jj)];
            
            mfprintf(fids,'Working on Angle %i, point [%f %f %f]\n',...
                kk,point);
            
            % Reference and data distances
            r_dist(jj,ii)=sqrt(sum((ref_camera_position-point).^2));
            d_dist(jj,ii)=sqrt(sum((camera_position(kk,:)-point).^2));
            
            ry_dist=ref_camera_position(2)-point(2);
            dy_dist=camera_position(kk,2)-point(2);
            
            rz_dist=ref_camera_position(3)-point(3);
            dz_dist=camera_position(kk,3)-point(3);
            
            rx_dist=ref_camera_position(1)-point(1);
            dx_dist=camera_position(kk,1)-point(1);
            
            % angle is in the r-y plane that lies the data point
            % (scattering plane)
            % need to establish direction so angle is right!
            r_ang(jj,ii)=acos(ry_dist/r_dist(jj,ii));
            d_ang(jj,ii)=acos(dy_dist/d_dist(jj,ii));
            
            r_phi=atan(rz_dist/rx_dist);
            d_phi=atan(dz_dist/dx_dist);
            
            % use B to find closest discrete angle in the data file
            [c r_ind]=min(abs(B.theta-r_ang(jj,ii)));
            r_d_ang=B.theta(r_ind);
            [c d_ind]=min(abs(B.theta-d_ang(jj,ii)));
            d_d_ang=B.theta(d_ind);
            
            [c r_phi_ind]=min(abs(phi-r_phi));
            r_d_phi=phi(r_phi_ind);
            [c d_phi_ind]=min(abs(phi-d_phi));
            d_d_phi=phi(d_phi_ind);
            
            mfprintf(fids,['...theta: %f, phi: %f, solid ang: %f,'...
                'sigma: %f, dist: %s\n'],r_d_ang,r_d_phi,half_angle,...
                sigma,dist_type);
            
            % corrected intensity ratio for distance
            i_ratio_c(jj,ii,kk);
            i_ratio_c(jj,ii,kk)=(d_dist(jj,ii)^2/...
                r_dist(jj,ii)^2)*i_ratio_c(jj,ii,kk);
            
            % Look in the local-to-this-run database.  If data is there
            % then a useable data record for this position is returned
            
            % Ref Angle     
            [sz_data_record_ref ind]=chk4data3(r_d_ang,r_d_phi,...
                half_angle,gamma_ref,xi,method,sigma,dist_type,l_sz_data);
            
            if isstruct(sz_data_record_ref)==0
                % if no useable data returned locally then the database
                % from the file is checked
                
                mfprintf(fids,'...No local data for: %f %f %f %f %s\n',...
                    r_d_ang,r_d_phi,half_angle,sigma,dist_type);
 
                [sz_data_record_ref ind]=chk4data3(r_d_ang,r_d_phi,...
                  half_angle,gamma_ref,xi,method,sigma,dist_type,sz3_data);
                
                if isstruct(sz_data_record_ref)==0
                    % if data is not in the file, calculate a new record
                    
                    mfprintf(fids,['...No file data for:'...
                        ' %f %f %f %f %s\n'],r_d_ang,r_d_phi,half_angle,...
                        sigma,dist_type);
                    
                    sz3_data=add_data3(r_d_ang,r_d_phi, sigma,dist_type,...
                        ini_filename,sz3_data,d_mean);
                    
                    % This will DEFINITELY return a record now
                    [sz_data_record_ref ind]=chk4data3(r_d_ang,r_d_phi,...
                        half_angle,gamma_ref,xi,method,sigma,...
                        dist_type,sz3_data);
                    
                    % update the disk copy of sz3_data
                    use_lockfile('sz_data_lock_file.txt','lock');
                    save(int_vs_size_filename,'sz3_data');
                    use_lockfile('sz_data_lock_file.txt','unlock');

                end
                
                % Add the "new" record to the local database
                l_sz_ind=l_sz_ind+1;
                l_sz_data(l_sz_ind)=sz_data_record_ref;
            end
            
            % Data Angle     
            [sz_data_record_dat ind]=chk4data3(d_d_ang,d_d_phi,...
                half_angle,gamma_ref,xi,method,sigma,dist_type,l_sz_data);
            
            if isstruct(sz_data_record_dat)==0
                % if no useable data returned locally then the database
                % from the file is checked
                
                mfprintf(fids,'...No local data for: %f %f %f %f %s\n',...
                    d_d_ang,d_d_phi,half_angle,sigma,dist_type);
               
                [sz_data_record_dat ind]=chk4data3(d_d_ang,d_d_phi,...
                  half_angle,gamma_ref,xi,method,sigma,dist_type,sz3_data);
                
                if isstruct(sz_data_record_dat)==0
                    % if data is not in the file, calculate a new record
                    
                    mfprintf(fids,['...No file data for:'...
                        ' %f %f %f %f %s\n'],d_d_ang,d_d_phi,half_angle,...
                        sigma,dist_type);
                     
                    sz3_data=add_data3(d_d_ang,d_d_phi,sigma,dist_type,...
                        ini_filename,sz3_data,d_mean);
                    % This will definitely return a record now
                    [sz_data_record_dat ind]=chk4data3(d_d_ang,d_d_phi,...
                        half_angle,gamma_ref,xi,method, sigma,dist_type,...
                        sz3_data);
                    
                    % update the disk copy of sz3_data
                    use_lockfile('sz_data_lock_file.txt','lock');
                    save(int_vs_size_filename,'sz3_data');
                    use_lockfile('sz_data_lock_file.txt','unlock');
                end
                % Add the "new" record to the local database
                l_sz_ind=l_sz_ind+1;
                l_sz_data(l_sz_ind)=sz_data_record_dat;
            end
            
            % No matter the initital state of the local database or file
            % database, sz_data_record now contains the required record.
            
            % Find index (or indexes) which are closest to given itensity
            % ratio, rat (and who doesn't like RAT as a variable)
            rat=i_ratio_c(jj,ii,kk);
            
            data_rec_i_ratio=sz_data_record_dat.intensity./...
                sz_data_record_ref.intensity;
            int_inds=find((data_rec_i_ratio(d_si)<rat & ...
                data_rec_i_ratio(d_si+1)>rat | ...
                data_rec_i_ratio(d_si)>rat & ...
                data_rec_i_ratio(d_si+1)<rat));
            
            if size(int_inds,2)==0
                mean_size{jj,ii,kk}(1)=0;
            else
                for ll=1:size(int_inds,2)
                    % cell array at each y-z location of the possible mean
                    % diameters
                    % linear interplation between nearest data points
                    
                    mean_size{jj,ii,kk}(ll)=interp1...
                        (data_rec_i_ratio(int_inds(ll):int_inds(ll)+1),...
                        sz_data_record_dat.d_mean(int_inds(ll):...
                        int_inds(ll)+1),rat);
                end
            end
        end
    end
end

fclose(out_fid);

% save output
%save(output_filename, 'act_msk_y', 'act_msk_z', 'mean_size','yo',...
%    'zo','sigma','half_angle','dist_type');
