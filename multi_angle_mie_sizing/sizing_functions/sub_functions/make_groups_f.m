function [groups angles group_mean_size group_std_dev]=...
    make_groups_f(mean_size,seed_index)
% A custum cluster analysis program.
% Finds groups of nearby data points and sorts them into groups
% Input is a cell-array, each entry is an array of possible droplet
% diameters, and which set to use as a starting "seed"
% Output is a cell array of the groups, a cell array which keeps track of
% which angle each data point originated from, the mean diameter of each 
% group, and standard deviation of each group
%
% This "sort" is essentially a specialized bucket sort with:
% 1) assumptions that no angular data set has any points that should be 
%    in the same group
% 2) buckets are added as needed when a data set is larger than the
%    comparison data set
% 3) buckets are added if the difference between the comparison set and the
%    data set is larger than the difference between elements in the data.
% 
% Code is written for specific data processing, not general use.  Example
% data is given below.  Documentation on method is found throughout the 
% code
%
%
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% -------------code begins------------
angs=size(mean_size,2);

% Use the seed_index data set as a starting comparison.  Thus there will be
% at least as many groups as there are entries in the seed set.  More
% groups are added as needed automatically
% The running mean of each group is used as the comparison against new data
% sets
group_mean_size=mean_size{seed_index};
ngroups=0;  % will be incremented as groups are added

kk_index_order=1:angs;
kk_index_order(1)=seed_index;
kk_index_order(seed_index)=1;

%for kk=1:angs
for kk=kk_index_order    
    n_el1=size(group_mean_size,2);
    n_el2=size(mean_size{kk},2);
    
    match_set1=group_mean_size;
    match_set2=mean_size{kk};
    adj_diff=[match_set2(1) abs(diff(match_set2)) match_set2(end)];
 
    % Finds the differece between every element in the data set compared to
    % the comparison set
    cnt=0;
    for ll=1:n_el2
        for mm=1:n_el1
            cnt=cnt+1;
            ind(cnt,:)=[ll mm];
            a_diff(cnt)=match_set2(ll)-match_set1(mm);
        end
    end
    % The differences are sorted in acending order
    [Y I]=sort(abs(a_diff));
    
    cnt=0;
    
    % Each time a data point is added to a group, it's indexes are stored
    % so that it cannot be added to any other groups
    used_ind=zeros(size(ind));
    for ll=1:size(a_diff,2)
        % Every index is checked in acending order.  If it is unused it
        % will be added to the closest group
        if (ismember(ind(I(ll),1), used_ind(:,1))==0 && ...
                ismember(ind(I(ll),2), used_ind(:,2))==0)
            if exist('groups')>0
                if ind(I(ll),2)<=size(groups,2)
                    grp_sz=size(groups{ind(I(ll),2)},2);
                else
                    grp_sz=0;
                end
            else
                grp_sz=0;
            end
            
            % Prevent jumping group membership over adjacent data
            % instead just create new group
            match_diff=abs(match_set2(ind(I(ll),1))-...
                match_set1(ind(I(ll),2)));
            if (match_diff<adj_diff(ind(I(ll),1)) && ...
                    match_diff<adj_diff(ind(I(ll),1)+1) )
                % Prevent adding values to a group that are way off
                % scales arbitrarily from 100% at small diameters down to 
                % 25% at larger diameters.
                s_fact=(-0.7/1000)*match_set2(ind(I(ll),1))+1;
                if s_fact<0.20
                    s_fact=0.20;
                end
                %s_fact=0.15;
                if match_diff<s_fact*match_set2(ind(I(ll),1))
                   groups{ind(I(ll),2)}(grp_sz+1)=match_set2(ind(I(ll),1));
                   angles{ind(I(ll),2)}(grp_sz+1)=kk;
                   used_ind(I(ll),:)=ind(I(ll),:);
                   cnt=cnt+1;
                   ngroups=size(groups,2);
                end
            end
        end
        
    end
    
    % If there are unused indexes (for case when data jumped as above or
    % cases where the current data set is larger than the comparison set)
    % then new groups are added here in acending order
    if cnt<n_el2
        grp_sz=0;
        for ll=1:size(I,2)
            if (ismember(ind(I(ll),1), used_ind(:,1))==0)
                unused_ind(I(ll),:)=ind(I(ll),:);
            end
        end
        
        add_ind=unique(unused_ind(:,1));
        for ll=1:size(add_ind,1)
            if add_ind(ll)~=0
                ngroups=ngroups+1;
                groups{ngroups}(grp_sz+1)=match_set2(add_ind(ll));
                angles{ngroups}(grp_sz+1)=kk;
            end
        end
    end
     
    % find mean of each group.  This carries forward and is used to compare
    % against each new data set.
    for ll=1:ngroups
        if size(groups{ll},2)>0  % protects against empty groups
            group_mean_size(ll)=mean(groups{ll});
        else
            group_mean_size(ll)=0;
        end
    end
   
    % Sort averages and re-arrange groups in acending order.  Order of the
    % data within each group is left as-is.
    [Z J]=sort(group_mean_size);
    tempgroups=groups;
    tempangles=angles;
    group_mean_size=Z;
    
    for ll=1:ngroups
        groups{ll}=tempgroups{J(ll)};
        angles{ll}=tempangles{J(ll)};
    end
    
    % Must clear certain variables each loop
    clear a_diff Y I ind unused_ind add_ind
    
end

% Calculate the standard deviation of each group.
for ll=1:size(groups,2)
    group_std_dev(ll)=std(groups{ll});
end

% move outliers from group to group
grp_num=1;
num_passes=0;
while grp_num<=size(groups,2)
    mvd_element=0;
    num_passes=num_passes+1;
    ll=1;
    while ll<=size(groups{grp_num},2);
        diff_m=abs(group_mean_size(grp_num)-groups{grp_num}(ll));
        
        if size(groups,2)==1
            diff_h=1e6;
            diff_l=1e6;
        elseif grp_num==1
            diff_h=abs(group_mean_size(grp_num+1)-groups{grp_num}(ll));
            diff_l=1e6;
        elseif grp_num==size(groups,2)
            diff_h=1e6;
            diff_l=abs(group_mean_size(grp_num-1)-groups{grp_num}(ll));
        else
            diff_h=abs(group_mean_size(grp_num+1)-groups{grp_num}(ll));
            diff_l=abs(group_mean_size(grp_num-1)-groups{grp_num}(ll));
        end
        
        % move element one group higher
        %if (diff_l>diff_h && diff_h<diff_m && ...
        %           diff_h<3*group_std_dev(grp_num+1))  % old method crit
        if (diff_l>diff_h && diff_h<diff_m)
            mvd_element=1;
            groups{grp_num+1}(size(groups{grp_num+1},2)+1)=...
                groups{grp_num}(ll);
            angles{grp_num+1}(size(angles{grp_num+1},2)+1)=...
                angles{grp_num}(ll);
            groups{grp_num}(ll)=[];
            angles{grp_num}(ll)=[];
        % move element one group lower
        %elseif (diff_h>diff_l && diff_l<diff_m && ...
        %           diff_l<3*group_std_dev(grp_num-1))  % old method crit
        elseif (diff_h>diff_l && diff_l<diff_m)
            mvd_element=1;
            groups{grp_num-1}(size(groups{grp_num-1},2)+1)=...
                groups{grp_num}(ll);
            angles{grp_num-1}(size(angles{grp_num-1},2)+1)=...
                angles{grp_num}(ll);
            groups{grp_num}(ll)=[];
            angles{grp_num}(ll)=[];
        else
            ll=ll+1;
        end
    end
    
    if (mvd_element==0 || num_passes>20)
        grp_num=grp_num+1;
        num_passes=0;
    end
    
    % refind mean of each group.
    for mm=1:size(groups,2)
        if size(groups{mm},2)>0
            group_mean_size(mm)=mean(groups{mm});
        else
            group_mean_size(mm)=0;
        end
    end
end

% delete outliers, one-pass
for mm=1:size(groups,2)
    ll=1;
    while ll<=size(groups{mm},2)
        cnt=0;
        clear elements
        elements(1)=1;
        for nn=1:size(groups{mm},2)
            if nn~=ll
                cnt=cnt+1;
                elements(cnt)=nn;
            end
        end
        std_dev_g=std(groups{mm}(elements));
        diff_dev=abs(groups{mm}(ll)-group_mean_size(mm));
        if (diff_dev>3*std_dev_g && std_dev_g>0)
            groups{mm}(ll)=[];
            angles{mm}(ll)=[];
        else
            ll=ll+1;
        end
    end
end

% re-find mean of each group.
for mm=1:size(groups,2)
    if size(groups{mm},2)>0
        group_mean_size(mm)=mean(groups{mm});
    else
        group_mean_size(mm)=0;
    end
end

% re-calculate the standard deviation of each group.
for ll=1:size(groups,2)
    group_std_dev(ll)=std(groups{ll});
end

