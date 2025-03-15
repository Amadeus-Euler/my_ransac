function [arclength,br,p_proj_sel,long_ip] = findStartandEndPerSet(points,radii,center,angle_eps)  
%     Purpose: To accurately calculate the radian of each cluster
%     Input:  points,radii,center,angle_eps
%     points is the point data, radii is the radii of the fitted circle,
%     center is the coordinate of circle center, angle_eps is the threshold of angel 

%     Output: [arclength,br,p_proj_sel,long_ip]   
%       arclength is the arclength of the longest cluster

    br=0;p_proj_sel=[];long_ip=[];
    [thetas,~]=cart2pol(points(:,1)-center(1),points(:,2)-center(2));
    [x_on_pre,y_on_pre]=pol2cart(thetas,repmat(radii,1,length(thetas))');
    x_on_circle=x_on_pre+center(1);
    y_on_circle=y_on_pre+center(2);
    p_proj=[x_on_circle,y_on_circle];

    epsilon=2*radii*sin(angle_eps*pi/360); %15 degree
    minPts = 5;
    [idx, ~] = dbscan(p_proj, epsilon, minPts);
    if(max(idx)==-1)
        arclength=-1;
        return
    end
    p_proj(:,3)=idx;
    p_proj_sel=p_proj(idx>0,:);
    %calculate radian of each cluster
    [thetas_sets,~]=cart2pol(p_proj_sel(:,1)-center(1),p_proj_sel(:,2)-center(2));
    thetas_sets(thetas_sets < 0) = thetas_sets(thetas_sets < 0) + 2*pi;
    thetas_sets(:,2)=p_proj_sel(:,3);
    lengths=unique(thetas_sets(:, 2));
    range_per_set = zeros(length(lengths), 3);
    for j = 1:length(lengths)
        current_class = lengths(j);
        set_scores = thetas_sets(thetas_sets(:, 2) == current_class, 1);
        set_range = max(set_scores) - min(set_scores);
        set_sd=std(set_scores);
        range_per_set(j, :) = [current_class, set_range, set_sd];
    end
    %Just one cluster
        % When there is only one class, it is difficult to judge whether it has
        % passed 0 degrees, so the line segment sum is used as the arc length to 
        % calculate the corresponding radian, and the original radian is compared
        % with the carefully designed threshold value as the judgment condition.
    if (max(idx)==1)
        order_thetas_sets=sortrows(thetas_sets,1);
        order_thetas_sets(:,2)=radii;
        [xx,yy]=pol2cart(order_thetas_sets(:,1),order_thetas_sets(:,2));
        xxyy=[xx,yy];
        xxyy_cuo=xxyy([length(xxyy),1:(length(xxyy)-1)],:);
        diff_mat=xxyy-xxyy_cuo;
        arclengths=zeros(length(xxyy),1);
        for k=1:length(xxyy)
            arclengths(k)=norm(diff_mat(k,:));
        end
        arclength=sum(arclengths)-max(arclengths);
        range_a_set=arclength/radii;
        if abs(range_a_set-range_per_set(1,2))>0.02 % pass the 0 degree
            rows=length(order_thetas_sets);
            diff=order_thetas_sets(:,1)-order_thetas_sets([2:rows,1],1);
            [~, minIndex] = min(diff);
            fenge=(order_thetas_sets(minIndex,1)+order_thetas_sets(minIndex+1,1))/2;
            largeset=thetas_sets(thetas_sets(:,1)>fenge,:);
            smallset=thetas_sets(thetas_sets(:,1)<fenge,:);
            range_breakset=2*pi-min(largeset(:,1))+max(smallset(:,1));
            arclength=range_breakset*radii;
            range_per_set(1,2)=range_breakset; % revamp BR
        end
        long_ip=p_proj_sel;
    else% more than one clusters
        if (sum(range_per_set(:,2))>2*pi-angle_eps*pi/180) % one of culsters passes the 0 degree
            maxrange=max(range_per_set(:,2));
            maxrange_idx=range_per_set(range_per_set(:,2)==maxrange,1);
            remain_idx=range_per_set(range_per_set(:,2)<maxrange,1);
            sets2selpoint=thetas_sets(ismember(thetas_sets(:,2),remain_idx));
            fenge=sets2selpoint(1);
            sets2break=thetas_sets(thetas_sets(:,2)==maxrange_idx,:);
            largeset=sets2break(sets2break(:,1)>fenge,:);
            smallset=sets2break(sets2break(:,1)<fenge,:);
            range_breakset=2*pi-min(largeset(:,1))+max(smallset(:,1));
            range_per_remainset=range_per_set(range_per_set(:,1)~=maxrange_idx,:);
            range_per_set=[range_per_remainset; 
            [maxrange_idx,range_breakset,2]
            ];
            max_id=range_per_set(range_per_set(:,2)==max(range_per_set(:,2)),1);
            long_ip=p_proj_sel(p_proj_sel(:,3)==max_id,:);
            arclength=max([range_per_remainset(:,2)',range_breakset])*radii;
        else% no culster pass the 0 degree
            arclength=max(range_per_set(:,2))*radii; 
            [~,max_id]=max(range_per_set(:,2));
            long_ip=p_proj_sel(p_proj_sel(:,3)==max_id,:);
        end
    end
    br=max(range_per_set(:,2));
end
