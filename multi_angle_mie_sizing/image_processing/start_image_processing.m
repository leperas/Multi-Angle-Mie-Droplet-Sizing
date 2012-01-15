clear all
close all

% Set the location of image procesing sub-functions: processor2_v82_f.m,
% redotlocator_v4.m, find_reg_offset.m, correct_perspective.m,
% make_image_correction.m, intensities2xyz_v2.m
current_path=pwd;   % gets current path for use in dialogs too
dir_str=strcat(current_path,'/sub_functions');

% Control the number of simultanious processes running
% This program is one process, so usually set this value to the
% number of processors +1 
num_processes_allowed=3;

%% Read in the initialization, no more user-settings below

current_path=pwd;   % gets current path for use in dialogs

warning('off','Images:initSize:adjustingMag');


[data_file, current_path, filterindex] = uigetfile('*.txt',...
    'Select .txt file with list of .ini files:', current_path ,...
    'MultiSelect', 'Off');

filename=strcat(current_path,data_file);

A=importdata(filename);

[status, result]=system('ps -ef |grep -c [M]ATLAB');
num_procs=str2num(result);

tic
for ii=1:size(A,1)
    if num_procs<num_processes_allowed
        
        system_string=strcat(...
            'nice -n 5 /opt/matlab/bin/matlab -r "cd ''',dir_str,...
            ''';processor2_v82_f(''',cell2mat(A(ii)),''');exit" &');
        
        [status0, result0]=system(system_string);
       
        pause(2);
        
        [status1, result]=system('ps -ef |grep -c [M]ATLAB');
        num_procs=str2num(result);
                
        mfprintf(1, '%s \n',cell2mat(A(ii)));
        mfprintf(1, ['Matlab procs: %i, Jobs started: %3.0f%%, '...
            'Status: %i,  %s\n\n'], num_procs, 100*ii/size(A,1),...
            status0, result0);
        
    end
    
    while num_procs>=num_processes_allowed
        
        pause(15);
        
        [status2, result]=system('ps -ef |grep -c [M]ATLAB');
        num_procs=str2num(result);
    end
    toc
end
