function [myID]=use_lockfile(shared_lock_file,action_item)
% FUNCTION [myID]=use_lockfile(shared_lock_file,action_item)
% Written based on code suggested by James O'Connell via MatlabCentral
% 
% action_item: either 'lock' or 'unlock'
%
% Use: in all programs sharing any common file, each program should
% call 'lock', to the same shared_lock_file, before accessing the 
% common file.  After accessing the common file, the program should call
% 'unlock'
%
% Function logic is as follows:
% To access a whatever common file, run the 'lock' action of this function.
% The function first waits until the lock-file doesn't exist, then 
% attempts to create the lock-file and write the uniquely generated 'myID'.
% The function then waits for a given amount of time and checks to 
% see if the lock-file still contains 'myID'. If it does, the function
% returns 'myID' and it's now safe to do whatever needs doing with the 
% common file.  If it doesn't, then another session has obviously written 
% the lock-file a fraction of a second after.  The function concedes and
% goes back to waiting for the lock-file to disappear before trying
% again.
% After recieving the return value from 'lock', the current session has
% the lock and it's safe to use the common file.  When finished with the
% common file, run the 'unlock' action.  The lock-file is deleted, and any
% other simultanious functions waiting to access the common file will get
% their chance.
%
% The returned 'myID' is a unique string generated to ensure that a lock
% actually occured; the return value is unused; a new string is defined on
% every call.
%
% Steve LePera, copyright 2012.  Free for distribution and use to anyone.
%   Use at your own risk!  NO WARRENTY!
%

% Throws error if a wait for the lockfile exceeds this number of seconds.
max_wait=30; 
tic;
switch action_item
    case 'lock'
        % Create lock file
        lockSuccess = 0;
        myID=strcat(num2str(rand),num2str(rand),num2str(cputime));
        while ~lockSuccess
            % Wait for lock-file to disappear
            %fprintf('Checking lock-file availability...');
            while exist(shared_lock_file,'file')
                pause(0.5);
                if toc>max_wait
                    error_str=['Stuck on wait-for-lockfile: ''' ...
                        shared_lock_file '''. Call ''unlock'' or delete'...
                        ' this file manually if certain no other '...
                        'process has access.'];
                    error(error_str);
                end
            end
            %fprintf(' Available.\n');
            % Attempt to create a lock
            %fprintf('Attempting to secure lock file...');
            fid = fopen(shared_lock_file,'w');
            fprintf(fid,'%s',myID);
            fclose(fid);
            pause(1);
            % Was I successful?
            if strcmp(textread(shared_lock_file,'%s'),myID)
                lockSuccess = 1;
                %fprintf(' Successful.\n');
            else
                %fprintf(' Unsuccessful, trying again.\n');
            end
        end
    case 'unlock'
        % Remove the lock-file
        unix(['rm ' shared_lock_file]);
        %fprintf('Lock file deleted.\n\n');
    otherwise
        error('Unknown action item.')
end
