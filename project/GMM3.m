data1=textread('C:\Users\EDA user\Desktop\GMM1\1_pointcloud.txt');
data2=textread('C:\Users\EDA user\Desktop\GMM1\30_pointcloud.txt');
data3=textread('C:\Users\EDA user\Desktop\GMM1\68_pointcloud.txt');


i=20;
L=[];

  
 pose1=textread('C:\Users\EDA user\Desktop\GMM1\1_T.txt');   
pose2=textread('C:\Users\EDA user\Desktop\GMM1\30_T.txt');
pose3=textread('C:\Users\EDA user\Desktop\GMM1\68_T.txt');

R1=pose1(1:3,1:3);
R2=pose2(1:3,1:3);
R3=pose3(1:3,1:3);

 T1=pose1(1:3,4);
 T2=pose2(1:3,4);
 T3=pose3(1:3,4);

  
  
  
GMModel1 = fitgmdist(data1,i,'RegularizationValue',0.0001);

mu1=GMModel1.mu';
mu1=R1'*mu1;
mu1=mu1-R1'*T1;
mu1=mu1';

Sigma1=GMModel1.Sigma;

for i1=1:i
Sigma1(:,:,i1)=R1'*Sigma1(:,:,i1)*R1;
end


GMModel2 = fitgmdist(data2,i,'RegularizationValue',0.0001);

mu2=GMModel2.mu';
mu2=R2'*mu2;
mu2=mu2-R2'*T2;
mu2=mu2';
Sigma2=GMModel2.Sigma;

for i1=1:i
Sigma2(:,:,i1)=R2'*Sigma2(:,:,i1)*R2;
end
GMModel3 = fitgmdist(data3,i,'RegularizationValue',0.0001);

mu3=GMModel3.mu';
mu3=R3'*mu3;
mu3=mu3-R3'*T3;
mu3=mu3';
Sigma3=GMModel3.Sigma;

for i1=1:i
Sigma3(:,:,i1)=R3'*Sigma3(:,:,i1)*R3;
end


mu=[mu3;mu2;mu1];
Sigma=cat(3,Sigma3,Sigma2,Sigma1);
size1=size(Sigma);
sizek=size1(1,3);
for v=1:sizek;
    zh(1:3,1:3,v)=Sigma(:,:,v);
    zh(4,1:3,v)=mu(v,1:3);
end

C=nchoosek(1:sizek,2) ;
size2=size(C);
size3=size2(1,1);
for g=1:size3;
    
    kl = kldivGaussian(zh(4,1:3,C(g,1)), zh(1:3,1:3,C(g,1)),zh(4,1:3,C(g,2)), zh(1:3,1:3,C(g,2)));
    L(g,1)=kl;
    
    
    
end



figure(1)
 hold on
  h1 = scatter3(data1(:,1),data1(:,2),data1(:,3),'.','b');
figure(2)
hold on
  h2 = scatter3(data2(:,1),data2(:,2),data2(:,3),'.','b');

  


  

 
 
  
  
 

  
  
 
  
  
for n=1:i*3
Mu =mu(n,1:3);
SIGMA=Sigma(1:3,1:3,n);
[V,D]=eig(SIGMA);
r =1.7;
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
xsect = xsect+Mu(1);
ysect = ysect+Mu(2);
zsect = zsect+Mu(3);
xsect1=real(xsect);
ysect1=real(ysect);
zsect1=real(zsect);
figure(3)
hold on
h3=surf(xsect1,ysect1,zsect1);
alpha(0.3)
xsect=[];
xsect1=[];
ysect=[];
ysect1=[];
zsect=[];
zsect1=[];

title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
end

  