function tex_out(dg, best_mean,best_sigma,best_dist,best_std,conf_num)
% Takes the output from analyze_sizes_f.m and returns table-formated 
% data to the screen.  Also calculates an 'Error percent' based on
% known dg, in micrometers.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

wavelength=514.5; % nanometers

% dimensions back to micrometers for mean diameter and standard dev
rl_x=wavelength*best_mean/(pi*1e3);
rl_std=2*wavelength*(best_std)/(pi*1000);

bd=flipud(best_dist);

mfprintf(1,' & %3.1f',flipud(rl_x));
mfprintf(1,' \\\\ \n');
mfprintf(1,'std dev ');
mfprintf(1,' & %4.2f',flipud(rl_std));
mfprintf(1,' \\\\ \n');
mfprintf(1,' $\\sigma$ ');
mfprintf(1,' & %3.1f',flipud(best_sigma));
mfprintf(1,' \\\\ \n');
mfprintf(1,'Distribution');
mfprintf(1,' & %s',bd{:});
mfprintf(1,' \\\\ \n');
mfprintf(1,' CN ');
mfprintf(1,' & %3.2f',flipud(conf_num));
mfprintf(1,' \\\\ \n');

err_d=100*abs(dg-flipud(rl_x))./dg;

%err_d(err_d==100)=0;
mfprintf(1,'Error \\%% ');
mfprintf(1,' & %3.2f',err_d);
mfprintf(1,' \\\\ \\hline \n \n');