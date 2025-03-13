function [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel]=myRanSac(data,sigma,itermax,max_r)
r=0;p=0;BR=0;
ip_on_circle_sel=[];
longest_ip_on_circle_sel=[];
a = data;
iter = 0;
m=length(a);
if m<6
    fprintf('point is not enough')
    r=0;p=0;
    return
end

largestarclength = 0;


while iter<itermax
    ran = randperm(m,3)';
    b = a(ran,:);
    [r1,p1] = ThreePoint2Circle(b(1,1:2), b(2,1:2), b(3,1:2));
    if r1>max_r*1.1|| r1<0.01 %limit the min and max
        iter=iter+1;
        continue
    end

    dis = sqrt(sum((a(:,1:2)-p1).^2,2));
    res = abs(dis - r1);
    d = a(res<sigma,:);

    if(isempty(d))
        iter=iter+1;
        continue
    end

    [arclength,br,p_proj_sel,long_ip]=findStartandEndPerSet(d,r1,p1,15);
    if arclength==-1
        iter=iter+1;
        continue
    end

    if (arclength>=largestarclength)
        largestarclength = arclength;
        r = r1; p = p1; BR=br;ip_on_circle_sel = p_proj_sel; longest_ip_on_circle_sel = long_ip;
    end
    iter = iter + 1;
end






