close;
clear;
clc;
% X=[];
% A1=csvread('C:\Users\EDA user\Desktop\pc\PointCloud0.csv');
% A2=csvread('C:\Users\EDA user\Desktop\pc\PointCloud1.csv');
% A3=csvread('C:\Users\EDA user\Desktop\pc\PointCloud2.csv');
% A4=csvread('C:\Users\EDA user\Desktop\pc\PointCloud3.csv');
% A5=csvread('C:\Users\EDA user\Desktop\pc\PointCloud4.csv');
% A6=csvread('C:\Users\EDA user\Desktop\pc\PointCloud5.csv');
% A7=csvread('C:\Users\EDA user\Desktop\pc\PointCloud6.csv');
% A8=csvread('C:\Users\EDA user\Desktop\pc\PointCloud7.csv');
% A9=csvread('C:\Users\EDA user\Desktop\pc\PointCloud8.csv');
% A10=csvread('C:\Users\EDA user\Desktop\pc\PointCloud9.csv');
% A11=csvread('C:\Users\EDA user\Desktop\pc\PointCloud10.csv');
% A12=csvread('C:\Users\EDA user\Desktop\pc\PointCloud11.csv');
% A13=csvread('C:\Users\EDA user\Desktop\pc\PointCloud12.csv');
% A14=csvread('C:\Users\EDA user\Desktop\pc\PointCloud13.csv');
% A15=csvread('C:\Users\EDA user\Desktop\pc\PointCloud14.csv');
% A16=csvread('C:\Users\EDA user\Desktop\pc\PointCloud15.csv');
% A17=csvread('C:\Users\EDA user\Desktop\pc\PointCloud16.csv');
% A18=csvread('C:\Users\EDA user\Desktop\pc\PointCloud17.csv');
% A19=csvread('C:\Users\EDA user\Desktop\pc\PointCloud18.csv');
% A20=csvread('C:\Users\EDA user\Desktop\pc\PointCloud19.csv');
% A21=csvread('C:\Users\EDA user\Desktop\pc\PointCloud20.csv');
% A22=csvread('C:\Users\EDA user\Desktop\pc\PointCloud21.csv');
% A23=csvread('C:\Users\EDA user\Desktop\pc\PointCloud22.csv');
% A24=csvread('C:\Users\EDA user\Desktop\pc\PointCloud23.csv');
% A25=csvread('C:\Users\EDA user\Desktop\pc\PointCloud24.csv');
% A26=csvread('C:\Users\EDA user\Desktop\pc\PointCloud25.csv');
% A27=csvread('C:\Users\EDA user\Desktop\pc\PointCloud26.csv');
% A28=csvread('C:\Users\EDA user\Desktop\pc\PointCloud27.csv');
% A29=csvread('C:\Users\EDA user\Desktop\pc\PointCloud28.csv');
% A30=csvread('C:\Users\EDA user\Desktop\pc\PointCloud29.csv');
% A31=csvread('C:\Users\EDA user\Desktop\pc\PointCloud30.csv');
% 
% X1=[A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;A13;A14;A15;A16;A17;A18;A19;A20;A21;A22;A23;A24;A25;A26;A27;A28;A29;A30;A31];
% for i=1:571350
%     X(i,1:3)=X1(10*i,1:3);
% end
load ('data.mat')
k=10;
figure
hold on
h1 = scatter3(data(:,1),data(:,2),data(:,3),'.','b');
 %GMModel = fitgmdist(data,k);
 X=data';
[z1,model,llh] = mixGaussEm(X,k);


% [dim,Num]=size(data);
% max_iter=10;%最大迭代次数
% min_improve=1e-4;% 提升的精度
% Ngauss=3;%混合高斯函数个数
% Pw=zeros(1,Ngauss);%保存权重
% mu= zeros(dim,Ngauss);%保存每个高斯分类的均值,每一列为一个高斯分量
% sigma= zeros(dim,dim,Ngauss);%保存高斯分类的协方差矩阵
% fprintf('采用K均值算法对各个高斯分量进行初始化\n');
% [cost,cm,cv,cc,cs,map] = vq_flat(data, Ngauss);%聚类过程  map:样本所对应的聚类中心
% mu=cm;%均值初始化
% for j=1:Ngauss
%    gauss_labels=find(map==j);%找出每个类对应的标签
%    Pw(j)= length(gauss_labels)/length(map);%类别为1的样本个数占总样本的个数 
%    sigma(:,:,j)  = diag(std(data(:,gauss_labels),0,2)); %求行向量的方差，只取对角线，其他特征独立，并将其赋值给对角线
% end
% 
% last_loglik = -Inf;%上次的概率
% % 采用EM算法估计GMM的各个参数
% if Ngauss==1,%一个高斯函数不需要用EM进行估计
%     sigma(:,:,1)  = sqrtm(cov(data',1));
%     mu(:,1)       = mean(data,2);
% else
%      sigma_i  = squeeze(sigma(:,:,:));
%      
%      iter= 0;
%      for iter = 1:max_iter
%           %E 步骤
%           %求每一样样本对应于GMM函数的输出以及每个高斯分量的输出，
%           sigma_old=sigma_i;
%           %E步骤。。。。。
%           for i=1:Ngauss
%           P(:,i)= Pw(i) * p_single(data, squeeze(mu(:,i)), squeeze(sigma_i(:,:,i)));%每一个样本对应每一个高斯分量的输出
%           end
%           s=sum(P,2);%
%         for j=1:Num
%             P(j,:)=P(j,:)/s(j);
%         end
%        %%%Max步骤
%         Pw(1:Ngauss) = 1/Num*sum(P);%权重的估计
%         %均值的估计
%         for i=1:Ngauss
%             sum1=0;
%             for j=1:Num
%              sum1=sum1+P(j,i).*data(:,j);
%             end
%           mu(:,i)=sum1./sum(P(:,i));
%         end
%        
%         %方差估计按照公式类似
%          %sigma_i
%          if((sum(sum(sum(abs(sigma_i- sigma_old))))<min_improve))
%              break;
%         end
%         
%         
%      end
%     
%      
% end
% 
% 


for n=1:8
 mu =model.mu(1:3,n);
 sigma=model.Sigma(1:3,1:3,n);
[V,D]=eig(sigma);
r =3;
xhalf = linspace(sqrt(r^2*D(1,1)),0,10);
Ninthalf = round(10/2);
zsect = zeros(10,Ninthalf);
ysect = zeros(10,Ninthalf);
for ti = 1:10
r2d = r^2 - xhalf(ti).^2/D(1,1);
ysect(ti,:) = linspace(0,sqrt(r2d*D(2,2)),Ninthalf);

zsect(ti,:) = sqrt((r2d - ysect(ti,:).^2/D(2,2) )*D(3,3));
xsect(ti,1:Ninthalf) = xhalf(ti);
end
zsect = real(zsect);
%x&gt;0,Z&gt;0
xsect = [xsect,xsect];
ysect = [ysect,fliplr(ysect)];
zsect = [zsect,-fliplr(zsect)];
%x&gt;0
xsect = [xsect,xsect];
ysect = [ysect,-fliplr(ysect)];
zsect = [zsect,fliplr(zsect)];
% make it a whole
xsect = [xsect;-flipdim(xsect,1)];
ysect = [ysect;flipdim(ysect,1)];
zsect = [zsect;flipdim(zsect,1)];
% rotate
[lr,lc] = size(xsect);
for ti = 1:lr
for tj = 1:lc
newcodi = [xsect(ti,tj),ysect(ti,tj),zsect(ti,tj)]*inv(V);
xsect(ti,tj) = newcodi(1);
ysect(ti,tj) = newcodi(2);
zsect(ti,tj) = newcodi(3);
end
end
% shift
xsect = xsect+mu(1);
ysect = ysect+mu(2);
zsect = zsect+mu(3);
xsect1=real(xsect);
ysect1=real(ysect);
zsect1=real(zsect);
surf(xsect1,ysect1,zsect1);
alpha(0.3)
xsect=[];
xsect1=[];
ysect=[];
ysect1=[];
zsect=[];
zsect1=[];

title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
end
hold off