% Simple GUI interface tool that calls redotlocator_v4.m, a good way
% to verify registration images are of correct exposure to yield correct
% spatial results.
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

clear all
close all

addpath ./sub_functions

current_path=pwd;
pathname=current_path;  % sets the startpath in the dialog to current path

keep_going=1;
while keep_going<2
    keep_going=keep_going+1;
    
    
    [data_file, pathname, filterindex] = uigetfile('*.*',...
        'Select files:', pathname ,'MultiSelect', 'Off');
    
    filename=strcat(pathname,data_file);
    img=imread(filename);
    
    data_file
    redotlocator_v4(img,1)
    
    button = MFquestdlg([ 0.6 , 0.1 ],'Open another file?'); 
    
    switch button
        case 'Yes'
            keep_going=1;
            close all
        otherwise
            keep_going=keep_going+1;
    end

    pause(1)
end



