%% import data

clear 
close all 
[fileName,filePath]=uigetfile('*.csv','Input Data-File'); 

fullPath = fullfile(filePath, fileName);

stem = readtable(fullPath); 
%information of tree and plot
plotno='yml02';treeno=3;lidarno=6;fagen=0.08;% fagen is stump height. if collar height is unknown, set it to be 0.

xyz = stem(:,1:3);
xyz = table2array(xyz);
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 'filled'); 

% normalization
chazhi=0;
xyz(:,3)=xyz(:,3)-min(xyz(:,3))+chazhi;
xyz(:,2)=xyz(:,2)-min(xyz(:,2));
xyz(:,1)=xyz(:,1)-min(xyz(:,1));
ptCloud = pointCloud(xyz);
% downsampling
ptCloud_d=pcdownsample(ptCloud,'gridAverage',0.005);
xyz=ptCloud_d.Location;
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 'filled');
figure
pcshow(ptCloud_d); 
H=max(xyz(:,3));

%%  calculate flatness and verticality
tic
cells = cell(length(xyz),7);
for i = 1:length(xyz)
    [index,dist] = findNeighborsInRadius(ptCloud_d,xyz(i,:),0.05);
    cells{i,1} = index;
    cells{i,2} = xyz(index,:);
    cells{i,3} = mean(dist);
end
toc

tic
for i = 1:length(xyz)
    if length(cells{i,1})>4
        mat = cell2mat(cells(i,2));
        covs = cov(mat); 
        [v,d] = eig(covs);

        d_list = diag(d);
        flatness = 1 - d_list(1)/sum(d_list);
        cells{i,4} = flatness;
        normal = v(:,1);
        cos_va = abs(dot(normal,[0,0,1]));
        rad = asin(cos_va);
        angel = rad2deg(rad);
        cells{i,5} = angel;
    else
        continue
    end
end
toc

C=cells(:,4);
C(cellfun('isempty', C)) = {0.7};
C = cell2mat(C);
condition1 = C<0.9;

C=cells(:,5);
C(cellfun('isempty', C)) = {90};
C = cell2mat(C);
condition2 = C>15;



condition = condition2;
% condition = condition1 | condition2;

mat = xyz(condition,:);
mat2 = xyz(~condition,:);
figure
pcshow(mat);
figure
pcshow(mat2);

%% GMM clustering 
model = fitgmdist(mat2, 2);
idx = cluster(model,mat2);



mat2(:,4)=idx;
mat_candi1=mat2(mat2(:,4) == 1,1:3);
mat_candi2=mat2(mat2(:,4) == 2,1:3);
mat_sel=mat_candi1;
if std(mat_candi1(:,1))>std(mat_candi2(:,1))
    mat_sel=mat2(mat2(:,4) == 2,1:3);
end



figure
scatter3(mat2(:,1), mat2(:,2), mat2(:,3), 6 ,'filled');
hold on
scatter3(mat_sel(:,1), mat_sel(:,2), mat_sel(:,3));

figure
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 6 ,'filled');
hold on
scatter3(mat_sel(:,1), mat_sel(:,2), mat_sel(:,3));

figure
pcshow(ptCloud_d);
hold on
scatter3(mat_sel(:,1), mat_sel(:,2), mat_sel(:,3),4,'filled','red');




covs = cov(mat_sel(:,1:3)); 
[v,~] = eig(covs);

orient = v(:,3);
cos_va = abs(dot(orient,[0,0,1]));
rad = asin(cos_va);
angel = 90-rad2deg(rad);
fprintf(['lean angle of stem：',num2str(angel)]); 

%% DBSACN Clustering :IF the GMM is good, jump this section
epsilon = 0.1;  
minPts = 100;     
tic
[idx, coreIndices] = dbscan(mat2, epsilon, minPts);
toc


figure
scatter3(mat2(:,1), mat2(:,2), mat2(:,3), 6, idx, 'filled');

% calculate H-W ratio
mat2(:,4)=idx;
HWs=zeros(1,max(idx));
for i=1:max(idx)
    sub=mat2(mat2(:,4) == i,:);
    xvar=var(sub(:,1));
    yvar=var(sub(:,2));
    zsd=std(sub(:,3));
    HW=zsd/(sqrt(xvar+yvar));
    HWs(i)=HW;
end
compare=1.8;
sets = find(HWs>compare);
condi = ismember(mat2(:,4),sets);
mat_sel=mat2(condi,1:3);

scatter3(mat2(:,1), mat2(:,2), mat2(:,3), 6, idx, 'filled');
hold on
scatter3(mat_sel(:,1), mat_sel(:,2), mat_sel(:,3));


covs = cov(mat_sel(:,1:3)); 
[v,~] = eig(covs);

orient = v(:,3);
cos_va = abs(dot(orient,[0,0,1]));
rad = asin(cos_va);
angel = 90-rad2deg(rad);
fprintf(['lean angle of stem：',num2str(angel)]); 
%% close all figures
close all 
%% diameter along the trunk
alg='ran'; %ran for ransac, plus for modified-ransac, ols for Pratt.
fenli=mat_sel;
theta = linspace(0, 2*pi, 100);
max_r=0.20;% max diameter   

hs = fagen:0.1:H;
sigma=0.004;
ratios=zeros(length(hs),1);
ranges=zeros(length(hs),1);
ds = zeros(length(hs),1);
ps = zeros(length(hs),2);
thick=0.02; %slice thickness
for i=1:length(hs)
    if hs(i)>6
        thick = 0.05;
    end
    if mod(i-1,20) == 0
        figure
    end
    span = [hs(i)-thick/2,hs(i)+thick/2];
    sub_fenli = fenli(span(1)<fenli(:,3)&span(2)>fenli(:,3),:);
    crosss = sub_fenli(:,1:2);
    crosss = unique(crosss,'rows');
    if (isempty(crosss))
        continue
    end

    
    if strcmp(alg,'plus')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = myRanSac(crosss,sigma,150,max_r);
    end
    if strcmp(alg,'ran')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = RanSac(crosss,sigma,150,max_r);
    end
    if strcmp(alg,'ols')
        [r,p,BR,ip_on_circle_sel,longest_ip_on_circle_sel] = CircleFitByPratt3(crosss,sigma);
    end

    if(isempty(ip_on_circle_sel))
        continue
    end
    ip_on_circle_sel_t=array2table(ip_on_circle_sel,'VariableNames',{'x','y','set'});
    k=mod(i,20);
    if k==0
        k=20;
    end
    subplot(4,5,k);
    scatter(crosss(:,1),crosss(:,2),20,'filled');
    axis equal;
    hold on
    scatter(ip_on_circle_sel(:,1), ip_on_circle_sel(:,2),10, 'MarkerEdgeColor','red');
    scatter(longest_ip_on_circle_sel(:,1), longest_ip_on_circle_sel(:,2),30, 'MarkerEdgeColor','blue');
    if r==0
        fprintf([',no',num2str(i),', height',num2str(hs(i)),'fail to fit\n']);
    end
    if r~=0
        xfit = p(1) + r * cos(theta);
        yfit = p(2) + r * sin(theta);
        plot(xfit, yfit, '--', 'Color', 'red', 'LineWidth', 1.5);
        title(['no:',num2str(i),', h=',num2str(round(hs(i),3)),', d=',num2str(round(r*200,3))]);
        hold off
        ds(i)=r*200;
        ps(i,:)=p;
    else
        title(['no:',num2str(i),', h=',num2str(hs(i)),', d=',num2str(r*200)]);
        hold off
    end
       %constrain the max radii
    max_r=median([ds(ds~=0);50])/200;
    ranges(i)=BR;
end
    hds=[hs',ds,ranges];
    hds_t=array2table(hds,"VariableNames",{'h','d','range'});

%% slice selection method 1
hds_t_ex0=hds_t(hds_t.d>0,:);
hds_t_ex0_higher5=hds_t_ex0(hds_t_ex0.h>0.6,:);
heights=height(hds_t_ex0_higher5);
hds_t_ex0_higher5.save=zeros(height(hds_t_ex0_higher5),1);
interval=10;


for i=1:ceil(heights/interval)
    hds_sub=hds_t_ex0_higher5(i*interval-interval+1:min(heights,i*interval),:);  
    sortedDf = sortrows(hds_sub,"range",'descend');
    choosed_h=sortedDf.h(1); 
    rowIndices = find(hds_t_ex0_higher5.h == choosed_h);
    hds_t_ex0_higher5.save(rowIndices)=1;
end
hds_t_ex0_higher5_sel=hds_t_ex0_higher5(hds_t_ex0_higher5.save==1,:);
hds_t_ex0_higher5_sel.save=[];
maxrange=[hds_t_ex0(hds_t_ex0.h<=0.6,:);hds_t_ex0_higher5_sel;{H,0,0}];% add (H,0)
maxrange.save=ones(height(maxrange),1);

fitdata=maxrange(maxrange.save==1,:);
figure
scatter(fitdata.h,fitdata.d, [],fitdata.range);
lightBlue = [0.7, 0.85, 1]; 
blue = [0, 0, 1];           
N = 256; 
customColormap = [linspace(lightBlue(1), blue(1), N)', ...
                  linspace(lightBlue(2), blue(2), N)', ...
                  linspace(lightBlue(3), blue(3), N)'];
colormap(customColormap); 


hold on
ft = fittype('smoothingspline');

opts = fitoptions(ft);
% opts.Weights = fitdata.ratio;
opts.SmoothingParam = 0.4; % [0(smooth), 1(unsmooth)]

[fitresult, ~] = fit(fitdata.h, fitdata.d, ft, opts);

plot(fitresult);
ylim([0, 40]); 
title('cubic curve by method 1');xlabel('h');ylabel('d');
figure
scatter(hds_t_ex0.h,hds_t_ex0.d,[],hds_t_ex0.range);
title('all h-d pairs');xlabel('h');ylabel('d');

%% slice selection method 2
hds_t_ex0=hds_t(hds_t.d>0,:);
hds_t_ex0.save=zeros(height(hds_t_ex0),1);
low_constrain=0.2;high_constrain=0.1;
hds_t_ex0.save(1:3)=[1,1,1];
for iter=4:height(hds_t_ex0)
    prev3=hds_t_ex0.d(iter-3);prev2=hds_t_ex0.d(iter-2);prev1=hds_t_ex0.d(iter-1);
    current_d=hds_t_ex0.d(iter);
    current_h=hds_t_ex0.h(iter);
    prev=median([prev3,prev2,prev1]);
    if current_h<1.3
        if abs(current_d-prev)<prev*low_constrain
            hds_t_ex0.save(iter)=1;
        else
            hds_t_ex0.save(iter)=0;
        end
    else
        if abs(current_d-prev)<prev*high_constrain
            hds_t_ex0.save(iter)=1;
        else
            hds_t_ex0.save(iter)=0;
        end
    end
end


hds_t_ex0(end+1,:)={H,0,1,1};
hds_t_ex0_sel=hds_t_ex0(hds_t_ex0.save==1,:);
fitdata=hds_t_ex0_sel;

figure
scatter(fitdata.h,fitdata.d, [],fitdata.range);
lightBlue = [0.7, 0.85, 1]; 
blue = [0, 0, 1];           
N = 256; 
customColormap = [linspace(lightBlue(1), blue(1), N)', ...
                  linspace(lightBlue(2), blue(2), N)', ...
                  linspace(lightBlue(3), blue(3), N)'];
colormap(customColormap); 

hold on
ft = fittype('smoothingspline');

opts = fitoptions(ft);
% opts.Weights = fitdata.ratio;
opts.SmoothingParam = 0.6; 

[fitresult, ~] = fit(fitdata.h, fitdata.d, ft, opts);

plot(fitresult);
ylim([0, 40]); 
title('cubic curve by method 2');xlabel('h');ylabel('d');