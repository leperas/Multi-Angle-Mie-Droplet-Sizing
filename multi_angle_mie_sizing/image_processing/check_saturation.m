% Simple GUI interface tool that checks the histogram of an image or set of
% images for saturation; checks all three color channels.
% A good way to verify that data images are of correct exposure
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>
clear all
close all

current_path=pwd;
pathname=current_path;  % sets the startpath in the dialog to current path

keep_going=1;
while keep_going<2
    keep_going=keep_going+1;
    
    [data_file, pathname, filterindex] = uigetfile('*.*',...
        'Select files:', pathname ,'MultiSelect', 'On');
    
    if ischar(data_file)==1
        max_loop=1;
    else
        max_loop=size(data_file,2);
    end
    
    for jj=1:max_loop
        if ischar(data_file)==1
            mfprintf(1,'%s \n',data_file);
            filename=strcat(pathname,data_file);
            img=imread(filename);
        else
            mfprintf(1,'%s \n',cell2mat(data_file(jj)));
            filename=strcat(pathname,cell2mat(data_file(jj)));
            img=imread(filename);
        end
        
        
        figure
        subplot(2,2,4)
        imshow(img)
        %% Histogram of data.  Check against too many saturated pixels.
        nonzero=0;
        for ii=1:3
            [counts bins]=imhist(img(:,:,ii));
            subplot(2,2,ii)
            title('histogram')
            hh=stem(bins,counts);
            % only count "signal", threshold of 25 (10%)
            sat_percent=100*counts(size(counts,1))/...
                sum(counts(26:size(counts,1)));
            if sat_percent>0.0
                mfprintf(1,['Channel %i, Percentage saturated pixels:' ...
                    '%f.  Pixels: %i \n'], ii, sat_percent,...
                    counts(size(counts,1)));
                nonzero=1;
            end
            set(gca, 'yscale', 'log')
            set(hh,'BaseValue',1)
        end
        if nonzero==0
            mfprintf(1,'  No saturated pixels \n');
        end
     
    end
    
    button = MFquestdlg([ 0.8 , 0.1 ],'Open another file?');  
    
    switch button
        case 'Yes'
            keep_going=1;
            close all
        otherwise
            keep_going=keep_going+1;
    end
    mfprintf(1,'\n');
    pause(1)
end

close all

