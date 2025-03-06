function [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel]=RanSac(data,sigma,itermax,max_r)
%初始赋值
r=0;p=0;BR=0;%初始赋值
ip_on_circle_sel=[];
longest_ip_on_circle_sel=[];
a = data;
% RANSCA参数
% 迭代次数
iter = 0;
% 查看圆数据的大小
m=length(a);
try
    if m<6
        fprintf('点数太少')
        r=0;p=0;
        return
    end
    % 误差参数
    lengest = 1;
    % 开始循环迭代
    while iter<itermax
        % 随机挑选三个点，三个点不重复
        % 拟合圆最少需要三个点，拟合直线最少需要两个
        % ran为索引编号
        ran = randperm(m,3)';
        % b为索引得到的点
        b = a(ran,:);
        % 根据随机得到的三个点，计算圆的半径和圆心
        [r1,p1] = ThreePoint2Circle(b(1,1:2), b(2,1:2), b(3,1:2));
        if r1>max_r*1.2|| r1<0.01%设置上下限
            iter=iter+1;
            continue
        end
        % 选择除了随机得到的三个点外的其他点
%         c = setdiff(a,b,"rows");
    %   构成圆的点也算进内点
        c=a;
        % 计算每个点到圆心的距离dis
        dis = sqrt(sum((c(:,1:2)-p1).^2,2));
        % 计算 dis和拟合圆的误差
        res = abs(dis - r1);
        % 选择小于误差的点，进入到内点中
        d = c(res<sigma,:);
        if(isempty(d))
            iter=iter+1;
            continue
        end
        len = length(d(:,1));
        [arclength,br,p_proj_sel,long_ip]=findStartandEndPerSet(d,r1,p1,15);
        if arclength==-1
            iter=iter+1;
            continue
        end

        % 存最优的圆
        if (len>=lengest)
            lengest = len;
             r = r1; p = p1; BR=br;ip_on_circle_sel = p_proj_sel; longest_ip_on_circle_sel = long_ip;
        end
        iter = iter + 1;
    end
catch
    fprintf('点数太少')
    r=0;p=0;
end




