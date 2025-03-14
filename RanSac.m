function [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel]=RanSac(data,sigma,itermax,max_r)

r=0;p=0;BR=0;
ip_on_circle_sel=[];
longest_ip_on_circle_sel=[];
a = data;
iter = 0;
m=length(a);
try
    if m<6
        fprintf('number of points is too small')
        r=0;p=0;
        return
    end
    lengest = 1;
    while iter<itermax
        ran = randperm(m,3)';
        b = a(ran,:);
        [r1,p1] = ThreePoint2Circle(b(1,1:2), b(2,1:2), b(3,1:2));
        if r1>max_r*1.2|| r1<0.01 %limit the min and max
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
        len = length(d(:,1));
        
%         if arclength==-1
%             iter=iter+1;
%             continue
%         end

        if (len>=lengest)
            lengest = len;
            [~,br,p_proj_sel,long_ip]=findStartandEndPerSet(d,r1,p1,15);
%             points=d;radii=r1;center=p1;angle_eps=15;
             r = r1; p = p1; BR=br;ip_on_circle_sel = p_proj_sel; longest_ip_on_circle_sel = long_ip;
        end
        iter = iter + 1;
    end
catch
    fprintf('number of points is too small')
    r=0;p=0;
end




