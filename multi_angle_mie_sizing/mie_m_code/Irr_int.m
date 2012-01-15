function [intensity Isr2_dist Qsr2_dist Usr2_dist Vsr2_dist matches]=...
             Irr_int(varargin)
% Irr_int(dgx, sig_gx, PDF_type, theta, phi, half_cone_ang, gamma_ref,...
%                        xi, method)
%
% Calculation of the itegrated irradiance over a designated solid angle 
% from log-normal or normal distribution of homogeneous spherical droplets
%
% Also returns the Stokes vector, (I, Q, U, V) and closest matching theta
%    (and size parameter in 'single' particle cases)
%
% All inputs are required, even in cases where they may be ignored due to
% a particular method.
%
% dgx is mean particle size parameter, see size parameter defined in code
% sig_gx is standard deviation size parameter, ditto
% PDF_type a string "single, "logn" or "norm", determines which PDF to use.
% theta is the scattering angle, all angles in radians
% phi is the angle between the scattering plane and referecne plane
% half_cone_ang is half the solid angle of the detector
% gamma_ref is the laser linear polarization angle to reference plane.  Set
%    to 'unpolarized' if the incident light is unpolarized
% xi is the polarizing filter transmission axis angle to reference plane
%   NOTE: assign string 'no filter' if no filter was used
% method is string, can be 'full', '1D', or 'fast'
% 'full': valid at any phi, 2D integration of detector area, allows filter
% '1D': valid at any phi, 1D integration (over theta), allows filter
% 'fast': vertical polarized incident light (scatter plane at 90 degrees)
%    1D integration, and xi set to 'no filter' (not needed at 90 degree
%    scatter plane)
% The scattering crosssection database variables A and B must be available
% A and B are data structures created by the program, create_Sx_database.m
% and contain the index of refraction of the particles, the size 
% parameter of particles, and S1, S2 scattering coeffecients
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

global A B

dgx=varargin{1}; % mean diameter size parameter
sig_gx=varargin{2}; % standard deviation size parameter
PDF_type=varargin{3}; % string, 'logn', 'norm', or 'single'.  
%                      Single is calculation for a individual particle of
%                      size dgx, sigma is ignored
theta=varargin{4}; % scattering angle, all angles in radians 
phi=varargin{5}; % angle of scattering plane to reference plane
half_cone_ang=varargin{6}; % half the solid cone angle of detector
gamma_ref=varargin{7}; % laser polarization angle to reference plane
                % assign string 'unpolarized' for unpolarized light
xi=varargin{8}; % polarizing filter transmission axis angle to reference
                % plane, assign string 'no filter' if no filter was used
method=varargin{9}; % can be 'full', '1D', or 'fast'
 % 'full': valid at any phi, 2D integration of detector area, allows filter
 % '1D': valid at any phi, 1D integration (over theta), allows filter
 % 'fast': vertical polarized incident light (scatter plane at 90 degrees)
 %    1D integration, and xi set to 'no filter' (not needed at 90 degree 
 %    scatter plane)

% DEFINITION OF SIZE PARAMETER!
% Before calling Irr_int.m, use this to convert the distribution params 
% into size params, where:
% dg is DIAMETER of a droplet, as given below in micrometers, and sig_g
% is the standard deviation of the normal distribution with mean dg, also
% in micrometers.
%dgx=2*pi*(dg/2)*1e3/wavelength;
%sig_gx=2*pi*(sig_g)*1e3/wavelength; 

 
if isnumeric(xi)==1
    filter=1;
else
    filter=0;
end

% Find characteristics of the database
dtheta=(B.theta(2)-B.theta(1));

tol=(B.theta(2)-B.theta(1))/2;

% Find closest match to requested angle within the database
[c ind] = min(abs(B.theta-theta));

% This is the value of theta at the closest point in the database
exact_theta=B.theta(ind);
matches=[exact_theta   0];  % [ closest theta      closest x]

% Find range of all inds within requested solid cone and center on closest
% point in the database
inds_cone=find(exact_theta-half_cone_ang-tol<=B.theta & ...
    exact_theta+half_cone_ang+tol>=B.theta);

% determines size of "mask" to apply, defining area of detector and making
% the database point exact_theta the exact center of the observation area
mask_size=min([ind-min(inds_cone) max(inds_cone)-ind]);

theta_inds=ind-mask_size:ind+mask_size;

% Generate grid of evenly spaced coordinates in theta-phi space
theta_matrix=B.theta(theta_inds);

switch method
    case '1D'
        phi_matrix=phi;
        irr_mat=ones(2*mask_size+1,1);
    case 'fast'
        phi_matrix=phi;
        irr_mat=ones(2*mask_size+1,1);
    case 'full'
        phi_matrix=zeros(1,2*mask_size+1,1);
        phi_matrix(mask_size+1)=phi;
        phi_matrix(1:mask_size)=phi-dtheta*(mask_size+1-(1:mask_size));
        phi_matrix(mask_size+2:2*mask_size+1)=phi+dtheta*(1:mask_size);
        
        % Creates circular mask for the receiving optics area.
        mask=getnhood(strel('disk',mask_size,0));
        
        % Create matrix to store irradiance results, same size as mask,
        %  and matrices for the Stokes vector at each location;
        irr_mat=ones(size(mask)).*mask;
    otherwise
        error('Method string is not valid or recognized.')
end

Isr2_dist=irr_mat;
Qsr2_dist=irr_mat;
Usr2_dist=irr_mat;
Vsr2_dist=irr_mat;

% Properites of the size distribution used

M = log((dgx^2)/sqrt(sig_gx+dgx^2));
V = sqrt(log(sig_gx/(dgx^2)+1));

switch PDF_type
    case 'norm'
        p=normpdf(B.x,dgx,sig_gx);
    case 'logn'
        p=lognpdf(B.x,M,V);
    case 'single'
        % Find closest match to requested angle within the database
        [c indx] = min(abs(B.x-dgx));
        matches(2)=B.x(indx);
    otherwise
        error('Unknown or unrecognized PDF type.')
end

% Find the irradiance at each location within the mask

for ii=1:size(theta_matrix,2)
    for jj=1:size(phi_matrix,2)
        % find input stokes vector
        if isnumeric(gamma_ref)==1
            gamma=gamma_ref-phi_matrix(jj);
            stokes=[1 cos(2*gamma) sin(2*gamma) 0]';
        else
            stokes=[1  0  0  0]';
        end
        % skip values where mask==0
        if irr_mat(ii,jj)~=0
            switch PDF_type
                case 'single'
                    S1=A(theta_inds(ii)).S1(indx);
                    S2=A(theta_inds(ii)).S2(indx);
                otherwise
                    S1=A(theta_inds(ii)).S1;
                    S2=A(theta_inds(ii)).S2;
            end
            S1s=abs(S1).^2;
           
            switch method
                case 'fast'
                    % assumes vertically polarized incident light and
                    % detector in 90 degree scattering plane
                    switch PDF_type
                        case 'single'
                            Isr2_dist(ii,jj)=S1s;
                        otherwise
                            Isr2=S1s;
                            
                            Isr2_dist(ii,jj)=trapz(B.x,Isr2.*p');
                    end
                    Qsr2_dist(ii,jj)=-Isr2_dist(ii,jj);
                    Usr2_dist(ii,jj)=0;
                    Vsr2_dist(ii,jj)=0;
                otherwise
                    S2s=abs(S2).^2;
                    
                    S11=0.5*(S2s+S1s);
                    S12=0.5*(S2s-S1s);
                    
                    % for unpolarized incident light
                    if isnumeric(gamma_ref)==0
                        % scattering matrix=[S11 S12
                        %                    S12 S11];
                        Isr2_d=S11*stokes(1)+S12*stokes(2);
                        Qsr2_d=S12*stokes(1)+S11*stokes(2);
                        switch PDF_type
                            case 'single'
                            otherwise
                                Isr2_d=trapz(B.x,Isr2_d.*p');
                                Qsr2_d=trapz(B.x,Qsr2_d.*p');
                        end
                        Usr2_d=0;
                        Vsr2_d=0;
                    else
                        % for linearly polarized incident light
                        S33=0.5*(conj(S2).*S1+S2.*conj(S1));
                        S34=1i*0.5*(S1.*conj(S2)-S2.*conj(S1));
                        
                        % scattering matrix=[S11 S12 0    0
                        %                    S12 S11 0    0
                        %                    0   0   S33  S34
                        %                    0   0   -S34 S33];
                        Isr2_d=S11*stokes(1)+S12*stokes(2);
                        Qsr2_d=S12*stokes(1)+S11*stokes(2);
                        Usr2_d=S33*stokes(3)+S34*stokes(4);
                        Vsr2_d=-S34*stokes(3)+S33*stokes(4);
                        
                        switch PDF_type
                            case 'single'
                            otherwise
%                                 % For plot later, comment out otherwise
%                                 Isr2=Isr2_d;
                                
                                % Integrate over the whole distribution
                                Isr2_d=trapz(B.x,Isr2_d.*p');
                                Qsr2_d=trapz(B.x,Qsr2_d.*p');
                                Usr2_d=trapz(B.x,Usr2_d.*p');
                                Vsr2_d=trapz(B.x,Vsr2_d.*p');
                        end
                    end
                    % If polarization filter used at detector
                    if filter==1
                        omega=asin(cos(theta_matrix(ii))*...
                            sin(phi_matrix(jj)));
                        xi_c=xi+omega;
                        
                        cos2x=cos(2*xi_c);
                        sin2x=sin(2*xi_c);
                        
                        Isr2_dist(ii,jj)=0.5*...
                          (Isr2_d+Qsr2_d*cos2x+Usr2_d*sin2x);
                        Qsr2_dist(ii,jj)=0.5*...
                          (Isr2_d*cos2x+Qsr2_d*cos2x^2+Usr2_d*cos2x*sin2x);
                        Usr2_dist(ii,jj)=0.5*...
                          (Isr2_d*sin2x+Qsr2_d*sin2x*cos2x+Usr2_d*sin2x^2);
                        Vsr2_dist(ii,jj)=0; 
                    else
                        Isr2_dist(ii,jj)=Isr2_d;
                        Qsr2_dist(ii,jj)=Qsr2_d;
                        Usr2_dist(ii,jj)=Usr2_d;
                        Vsr2_dist(ii,jj)=Vsr2_d;
                    end
                    
            end
        end
    end
end

% Use trapz to perform descrete integration of the irradiance matrix, 
% returns intensity at detector

if (size(Isr2_dist,1)>1 && size(Isr2_dist,2)>1)
    % note='2d' % intensity, Watts/rad^2
    intensity = trapz(phi_matrix,trapz(theta_matrix,Isr2_dist));
elseif size(Isr2_dist,1)>1 && size(Isr2_dist,2)==1
    % note='1d' % Watts/rad
    intensity = trapz(theta_matrix,Isr2_dist);
else
    % note='single point' % Watts
    intensity=Isr2_dist;
end

% %% Uncomment to generate some informative plots during each run.
% % Works best with 0 cone angle and 'norm' or 'logn'
% 
% wavelength=514.5; % nanometers
% 
% % Conversion back to actual diameters
% rl_x=2*wavelength*B.x/(2*pi*1e3);
% dg=2*wavelength*dgx/(2*pi*1e3);
% rl_sig=wavelength*sig_gx/(2*pi*1e3);
% 
% figure
% range=1:size(Isr2,1);
% 
% plot(rl_x(range),Isr2(range)/max(Isr2(range)),'--',...
%       rl_x(range),p(range)/max(p(range)),'-')
% title_str=strcat('Particle Distribution, compared to Is, theta=',...
%       num2str(180*theta/pi));
% title(title_str)
% xlabel('Particle diameter, micrometers')
% ylabel('Normalized irradiance, normalized probability')
% legend_str=strcat('p, mn=',num2str(dg),' dv=',num2str(rl_sig));
% legend('Is', legend_str,'Location','NorthWest')
% ylim([0 1.1])
% %saveas(gcf,'scattering_pdf_1.eps','epsc')
%    
% figure
% plot(rl_x(range),(Isr2(range).*p(range)')/max(Isr2(range).*p(range)'),...
%         '--',rl_x(range),p(range)/max(p(range)),'-')
% title_str=strcat('Particle Dist, compared to p*Is.  Intensity=',...
%         num2str(floor(intensity)));
% title(title_str)
% xlabel('Particle diameter, micrometers')
% ylabel('Normalized irradiance, normalized probability')
% legend('p * Is', legend_str,'Location','NorthWest')
% ylim([0 1.1])
% %saveas(gcf,'scattering_pdf_2.eps','epsc')
% %% End informative plots

return

