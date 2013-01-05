function [best_mean best_std best_sigma best_dist match_fraction ...
    confidence_number]=...
    analyze_sizes_f(filename)
% Reads the ini_filname and estimates droplet sizes using input intensity
% ratios
% Process the given processing initialization file, or
% if ini_filename is a number, open a GUI
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>
tic;

% location of Irr_int.m
addpath ../../mie_m_code

if isnumeric(filename)==1
    current_path=pwd;   % gets current path for use in dialogs
    [data_file, current_path, filterindex] = uigetfile('*.ini',...
        'Select .ini file:', current_path ,'MultiSelect', 'Off');
    ini_filename=strcat(current_path,data_file);
else
    result=strfind(filename,'/');
    if numel(result)==0
        current_path=strcat(pwd,'/');
        ini_filename=filename;
    else
        lst_sl=result(size(result,2));
        current_path=filename(1:lst_sl);
        data_file=filename(lst_sl+1:length(filename));
        ini_filename=strcat(current_path,data_file);
    end
end

[keys,sections,subsections] = inifile(ini_filename,'readall');

key_assignment_string=strcat(keys(:,2),'_',keys(:,3),'=',keys(:,4),';');

% Assign values to the initalization file variables
for ii=1:size(key_assignment_string,1)
    eval(char(key_assignment_string(ii)))
end

% create output directory
[pathstr, name_ofile, file_ext, versn] = fileparts(output_output_filename);
[success message messageid]=mkdir(pathstr);
if success~=1
    error(strcat('Failure to create output directory: ',message))
end

% create temporary data file which contains limited data block
temp_datafile_name=strcat(input_ratio_file,'.tmp');

% load the data
load(input_ratio_file)

global A B
load(config_cross_section_database_filename)

ic=0;
for ii=input_y_indexes
    ic=ic+1;
    y_pos(ic)=act_msk_y(ii);
    jc=0;
    for jj=input_z_indexes
        jc=jc+1;
        z_pos(jc)=act_msk_z(jj);
        kc=0;
        for kk=input_camera_indexes
            kc=kc+1;
            int_ratio(jc,ic,kc)=i_ratio(jj,ii,kk);
            cam_loc(kc,:)=camera_position(kk,:);
        end
    end
end

% clear old variables and re-assign
clear act_msk_y act_msk_z i_ratio camera_position
act_msk_y=y_pos;
act_msk_z=z_pos;
i_ratio=int_ratio;
camera_position=cam_loc;

switch config_use_poly
    case 'Yes'
        % curve fitted i_ratios
        if size(camera_position,1)>1
            for ii=1:size(act_msk_y,2)
                for jj=1:size(act_msk_z,2)
                    for kk=1:size(camera_position,1)
                        point=[0, act_msk_y(ii) act_msk_z(jj)];
                        d_dist=sqrt(sum((camera_position(kk,:)-point).^2));
                        dy_dist=camera_position(kk,2)-point(2);
                        d_ang(kk)=acos(dy_dist/d_dist);
                    end
                    pp=spline(180*d_ang/pi(),squeeze(i_ratio(jj,ii,:)));
                    data_range_pts=...
                        [180*d_ang(1)/pi():0.05:180*d_ang(end)/pi()];
                    curve_pts=ppval(pp,data_range_pts);
                    
                    [p S smu]=polyfit(data_range_pts,curve_pts,...
                        config_poly_order);
                    
                    for kk=1:size(camera_position,1)
                        i_ratio(jj,ii,kk)=polyval(p,...
                            ((180*d_ang(kk)/pi())-smu(1))/smu(2));
                    end
                end
            end
        end
    otherwise
end
% temporary file with limited data set created
save(temp_datafile_name, 'act_msk_y', 'act_msk_z', 'i_ratio','yo',...
    'zo','ref_camera_position','camera_position');


%% Begin determination of sizes

ny=size(act_msk_y,2);
nz=size(act_msk_z,2);
ns=size(config_try_sigma,2);
n_angs=size(camera_position,1);

all_angs_mean_size=zeros(nz,ny,ns);
all_angs_std_dev=zeros(nz,ny,ns);

for mm=1:size(config_try_sigma,2)
    
    [mean_size]=int2size3(temp_datafile_name,...
       output_output_filename, config_int_vs_size_filename,...
       input_half_angle, num2str(input_gamma_ref), num2str(input_xi),...
       input_scatter_calc_meth,...
       2*pi*config_try_sigma(mm)*1e3/input_wavelength,...
       config_try_dist_type{mm},ini_filename);
    
  
   
    for ii=1:ny
        for jj=1:nz
            
            for kk=1:n_angs
                mean_size_groups{kk}=mean_size{jj,ii,kk};
            end
            
            [groups angles group_mean_size group_std_dev]=...
                make_groups_f(mean_size_groups,input_seed_index);
            
            
            for ll=1:size(groups,2)
                grp_elements(ll)=size(groups{ll},2);
            end
            
            % create a sorting matrix to pick out most likely diameter.
            % most likely is group with the most elements matching, 
            % smallest normalized standard deviation, and largest diameter.
            sorting_matrix=zeros(size(groups,2),4);
            
            % 1st column mean diameter in each group
            sorting_matrix(:,1)=group_mean_size';
            % 2nd column std deviation in each group
            sorting_matrix(:,2)=group_std_dev';
            % 3rd column number elements in group
            sorting_matrix(:,3)=grp_elements';
            % 4th column normalized std dev in group
            for nn=1:size(groups,2)
                if group_mean_size(nn)==0
                    sorting_matrix(nn,4)=0;
                else
                    sorting_matrix(nn,4)=...
                        group_std_dev(nn)'./group_mean_size(nn)';
                end
            end
            
            [sorted_matrix II]=sortrows(sorting_matrix,[-3 4 2 -1]);
           
            ll=1;
            while (ll<size(groups,2) && sorted_matrix(ll,1)<0.1 )
                ll=ll+1;
            end
            
            all_angs_mean_size(jj,ii,mm)=sorted_matrix(ll,1);
            all_angs_std_dev(jj,ii,mm)=sorted_matrix(ll,2);
            all_angs_ang_inds{jj,ii,mm}=angles{II(ll)};
            all_angs_num_matches(jj,ii,mm)=sorted_matrix(ll,3);
            
            clear groups angles group_mean_size group_std_dev ...
                mean_size_groups grp_elements
        end
    end
    
end

best_mean=zeros(nz,ny);
best_std=zeros(nz,ny);
best_sigma=zeros(nz,ny);
confidence_number=zeros(nz,ny);

for ii=1:ny
    for jj=1:nz
        sorting_matrix2=...
            zeros(size(squeeze(all_angs_mean_size(jj,ii,:)),1),4);
        
        % 1st column mean diameter in each group
        sorting_matrix2(:,1)=squeeze(all_angs_mean_size(jj,ii,:))';
        % 2nd column std deviation in each group
        sorting_matrix2(:,2)=squeeze(all_angs_std_dev(jj,ii,:))';
        % 3rd column number elements in group
        sorting_matrix2(:,3)=squeeze(all_angs_num_matches(jj,ii,:))';
        % 4th column normalized std dev in group
        for nn=1:size(squeeze(all_angs_mean_size(jj,ii,:)),1)
            if squeeze(all_angs_mean_size(jj,ii,nn))==0
                sorting_matrix2(nn,4)=0;
            else
                sorting_matrix2(nn,4)=...
                    squeeze(all_angs_std_dev(jj,ii,nn))'./...
                    squeeze(all_angs_mean_size(jj,ii,nn))';
            end
        end
        
        [sorted_matrix2 II2]=sortrows(sorting_matrix2,[-3 4 2 -1]);
        
        nn=1;
        while(sorted_matrix2(nn,3)>=round(0.8*sorted_matrix2(1,3)) && ...
                nn<size(squeeze(all_angs_mean_size(jj,ii,:)),1))
             nn=nn+1;
        end
                
        resort_size=nn-1;
        if resort_size<2
            resort_size=1;
        end
        
        [sorted_matrix3 II3]=...
            sortrows(sorted_matrix2(1:resort_size,:),[4 -3 2 -1]);
        
        ll=1;
        while (ll<resort_size && sorted_matrix2(ll,1)<0.1 )
            ll=ll+1;
        end
        
        best_mean(jj,ii)=sorted_matrix3(ll,1);
        best_std(jj,ii)=sorted_matrix3(ll,2);
        best_sigma(jj,ii)=config_try_sigma(II2(II3(ll)));
        best_dist{jj,ii}=config_try_dist_type{II2(II3(ll))};
        match_fraction(jj,ii)=sorted_matrix3(ll,3)/n_angs;
        
        
        point=[0, act_msk_y(ii) act_msk_z(jj)];
        r_dist=sqrt(sum((ref_camera_position-point).^2));
        ry_dist=ref_camera_position(2)-point(2);
        r_ang=acos(ry_dist/r_dist);
        
        rz_dist=ref_camera_position(3)-point(3);
        rx_dist=ref_camera_position(1)-point(1);
        
        r_phi=atan(rz_dist/rx_dist);
        
        % use B to find closest discrete angle in the data file
        [c r_ind]=min(abs(B.theta-r_ang));
        r_d_ang=B.theta(r_ind);
        
        bst_num=0;
        ths_num=0;
        ex_cnt=0;
        for kk=1:n_angs
            d_dist=sqrt(sum((camera_position(kk,:)-point).^2));
            dy_dist=camera_position(kk,2)-point(2);
            d_ang(kk)=acos(dy_dist/d_dist);
            
            % difference from current angle to 139
            delt_t=180*abs(d_ang(kk)-2.426)/pi();
            bst_num=bst_num+1/(1+delt_t);
            
            % find indexes of the "excluded points"
            if max(ismember(...
                    all_angs_ang_inds{jj,ii,II2(II3(ll))},kk))==0
                ex_cnt=ex_cnt+1;
                ex_inds(ex_cnt)=kk;
            else
                ths_num=ths_num+1/(1+delt_t);
            end
        end
        
        % Weighted fraction; matches closer to 139 degrees are worth more
        confidence_number(jj,ii)=ths_num/bst_num;
        
        % Make a plot (won't plot failures)
        if (output_make_a_plot==1 && best_mean(jj,ii)>1e-7)
            
            figure
            hold on
            plot(180*d_ang/pi(),squeeze(int_ratio(jj,ii,:)),'ko')
            legend_str{1}='Data';
            lgnd_cnt=1;
            
            if ex_cnt>0
                % plot excluded points
                plot(180*d_ang(ex_inds)/pi(),...
                    squeeze(int_ratio(jj,ii,ex_inds)),'rx')
                lgnd_cnt=lgnd_cnt+1;
                legend_str{lgnd_cnt}='No match (by Grouping)';
            end
            
            if n_angs>1
                pp=spline(180*d_ang/pi(),squeeze(i_ratio(jj,ii,:)));
                data_range_pts=...
                    [180*d_ang(1)/pi():0.05:180*d_ang(end)/pi()];
                curve_pts=ppval(pp,data_range_pts);
                
                [p S smu]=...
                    polyfit(data_range_pts,curve_pts,config_poly_order);
                poly_pts=polyval(p,(data_range_pts-smu(1))/smu(2));
            end
            
            
            sig_gx=2*pi*(best_sigma(jj,ii))*1e3/input_wavelength;
            for kk=1:size(output_plot_angle_range,2)
                [IntI Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
                    Irr_int(best_mean(jj,ii), sig_gx, best_dist{jj,ii},...
                    r_d_ang, r_phi, input_half_angle, input_gamma_ref,...
                    input_xi,input_scatter_calc_meth);
                
                [AIntI Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
                    Irr_int(best_mean(jj,ii), sig_gx, best_dist{jj,ii},...
                    output_plot_angle_range(kk)*pi/180, r_phi,...
                    input_half_angle, input_gamma_ref,input_xi,...
                    input_scatter_calc_meth);
                
                R_II(kk)=AIntI/IntI;
                
            end
            plot(output_plot_angle_range,R_II,'b-')
            lgnd_cnt=lgnd_cnt+1;
            legend_str{lgnd_cnt}=strvcat(['Mie theory: d=' ...
               num2str(2*input_wavelength*best_mean(jj,ii)/(2*pi*1e3))],...
               ['\sigma=' num2str(best_sigma(jj,ii)) ', ' ...
               best_dist{jj,ii} ' distribution']);
            
%            % Uncomment to also plot a "known" Mie theory curve
%            known_d=7.125;
%            known_dist='logn';
%            known_sig=10;
%            kdx=2*pi*(known_d/2)*1e3/input_wavelength;
%            sig_gx=2*pi*(known_sig)*1e3/input_wavelength;
%             for kk=1:size(output_plot_angle_range,2)
%               [IntIk Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
%                     Irr_int(kdx, sig_gx, known_dist,...
%                     r_d_ang, r_phi, input_half_angle, input_gamma_ref,...
%                     input_xi,input_scatter_calc_meth);
%                 
%              [AIntIk Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
%                     Irr_int(kdx, sig_gx, known_dist,...
%                     output_plot_angle_range(kk)*pi/180, r_phi,...
%                     input_half_angle, input_gamma_ref,input_xi,...
%                     input_scatter_calc_meth);
%                 
%                R_IIk(kk)=AIntIk/IntIk;
%                 
%             end
%             plot(output_plot_angle_range,R_IIk,'g-')
%             lgnd_cnt=lgnd_cnt+1;
%             legend_str{lgnd_cnt}=strvcat(['Known Mie theory: d=' ...
%                num2str(known_d)],...
%                ['\sigma=' num2str(known_sig) ', ' ...
%                known_dist ' distribution']);
%            % End "Known" section
           
            switch config_use_poly
                case 'Yes'
                    if n_angs>1
                        plot(data_range_pts,poly_pts,'r-')
                        lgnd_cnt=lgnd_cnt+1;
                        legend_str{lgnd_cnt}='Poly fit to data';
                    end
                otherwise
            end
            
            % find error at each data point
            for kk=1:n_angs
                [c d_ind]=min(abs(B.theta-d_ang(kk)));
                dd_ang(kk)=180*B.theta(d_ind)/pi();
                [c ind]=min(abs(output_plot_angle_range-dd_ang(kk)));
                err_int(kk)=(R_II(ind)-i_ratio(jj,ii,kk));
            end
            
            num_outliers=0;
            lst_outliers=-1;
            err_cnt=0;
            tmp_err_int=err_int;
            while num_outliers>lst_outliers
                lst_outliers=num_outliers;
                
                err_mean=mean(tmp_err_int);
                err_std=std(tmp_err_int);
                
                % An outlier will arbitrarily be considered 
                % 3 standard deviations
                for kk=1:n_angs
                   if abs(err_int(kk))>abs(3*err_std)
                        err_cnt=err_cnt+1;
                        err_inds(err_cnt)=kk;
                   end
                end
                
                % remove outliers from calculation
                tmp_err_int=err_int;
                if err_cnt>0
                    tmp_err_int(err_inds)=[];
                    num_outliers=size(unique(err_inds),2);
                end
            end
            
            if err_cnt>0
                lgnd_cnt=lgnd_cnt+1;
                legend_str{lgnd_cnt}='Outlier, Mie Theory';
                plot(180*d_ang(unique(err_inds))/pi(),...
                    squeeze(int_ratio(jj,ii,unique(err_inds))),'ms')
            end
            
            
            % location given as relative to the "figure origin"
            title_str=['Angle vs. Intensity Ratio, ' 'y=' ...
               num2str(act_msk_y(ii)-yo) '  z=' num2str(act_msk_z(jj)-zo)];
            title(title_str)
            xlabel('Angle, degress')
            ylabel('Intensity Ratio')
            legend(legend_str,'Location','Best')
            
            clear legend_str err_inds ex_inds
         end
        
    end
end

processing_time=toc;

% save results
save(output_output_filename,'act_msk_y','act_msk_z','yo','zo',...
    'best_mean','best_std','best_sigma','best_dist','all_angs_std_dev',...
    'all_angs_mean_size','match_fraction','confidence_number',...
    'processing_time');
%cleanup temp file
delete(temp_datafile_name);
