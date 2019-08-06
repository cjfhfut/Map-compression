data1=textread('C:\Users\EDA user\Desktop\GMM\1_pointcloud.txt');
data2=textread('C:\Users\EDA user\Desktop\GMM\6_pointcloud.txt');
data3=textread('C:\Users\EDA user\Desktop\GMM\11_pointcloud.txt');
data4=textread('C:\Users\EDA user\Desktop\GMM\16_pointcloud.txt');
data5=textread('C:\Users\EDA user\Desktop\GMM\21_pointcloud.txt');
data6=textread('C:\Users\EDA user\Desktop\GMM\26_pointcloud.txt');
data7=textread('C:\Users\EDA user\Desktop\GMM\31_pointcloud.txt');
data1(:,4)=[];
data2(:,4)=[];
data3(:,4)=[];
data4(:,4)=[];
data5(:,4)=[];
data6(:,4)=[];
data7(:,4)=[];


i=10;
L=[];

  
  pose1=textread('C:\Users\EDA user\Desktop\GMM\1_pose.txt');
  pose2=textread('C:\Users\EDA user\Desktop\GMM\6_pose.txt');
  pose3=textread('C:\Users\EDA user\Desktop\GMM\11_pose.txt');
  pose4=textread('C:\Users\EDA user\Desktop\GMM\16_pose.txt');
  pose5=textread('C:\Users\EDA user\Desktop\GMM\21_pose.txt');
  pose6=textread('C:\Users\EDA user\Desktop\GMM\26_pose.txt');
  pose7=textread('C:\Users\EDA user\Desktop\GMM\31_pose.txt');
            

    q1=pose1(1,4:7);
     q2=pose2(1,4:7);
     q3=pose3(1,4:7);
     q4=pose4(1,4:7);
     q5=pose5(1,4:7);
     q6=pose6(1,4:7);
     q7=pose7(1,4:7);
    
%     R1=[ 2*q1(1,1).^2-1+2*q1(1,2)^2    2*(q1(1,2)*q1(1,3)-q1(1,1)*q1(1,4)) 2*(q1(1,2)*q1(1,4)+q1(1,1)*q1(1,3));
%     2*(q1(1,2)*q1(1,3)+q1(1,1)*q1(1,4)) 2*q1(1,1)^2-1+2*q1(1,3)^2     2*(q1(1,3)*q1(1,4)-q1(1,1)*q1(1,2));
%     2*(q1(1,2)*q1(1,4)-q1(1,1)*q1(1,3)) 2*(q1(1,3)*q1(1,4)+q1(1,1)*q1(1,2)) 2*q1(1,1)^2-1+2*q1(1,4)^2];
% 
%     R2=[ 2*q2(1,1).^2-1+2*q2(1,2)^2    2*(q2(1,2)*q2(1,3)-q2(1,1)*q2(1,4)) 2*(q2(1,2)*q2(1,4)+q2(1,1)*q2(1,3));
%     2*(q2(1,2)*q2(1,3)+q2(1,1)*q2(1,4))  2*q2(1,1)^2-1+2*q2(1,3)^2     2*(q2(1,3)*q2(1,4)-q2(1,1)*q2(1,2));
%     2*(q2(1,2)*q2(1,4)-q2(1,1)*q2(1,3)) 2*(q2(1,3)*q2(1,4)+q2(1,1)*q2(1,2)) 2*q2(1,1)^2-1+2*q2(1,4)^2];
R1=quat2rot(q1);
R2=quat2rot(q2);
R3=quat2rot(q3);
R4=quat2rot(q4);
R5=quat2rot(q5);
R6=quat2rot(q6);
R7=quat2rot(q7);

  
R71=R7*inv(R1);
R72=R7*inv(R2);
R73=R7*inv(R3);
R74=R7*inv(R4);
R75=R7*inv(R5);
R76=R7*inv(R6);



     T1=0; 
     T2=0; 
     T3=0; 
     T4=0; 
     T5=0; 
     T6=0; 
     
GMModel1 = fitgmdist(data1,i,'RegularizationValue',0.0001);

mu1=GMModel1.mu';
mu1=R71*mu1;
mu1=mu1';
mu1=mu1+T1;
Sigma1=GMModel1.Sigma;

for i1=1:i
Sigma1(:,:,i1)=R71*Sigma1(:,:,i1)*R71';
end
GMModel2 = fitgmdist(data2,i,'RegularizationValue',0.0001);

mu2=GMModel2.mu';
mu2=R72*mu2;
mu2=mu2';
mu2=mu2+T2;
Sigma2=GMModel2.Sigma;

for i1=1:i
Sigma2(:,:,i1)=R72*Sigma2(:,:,i1)*R72';
end
GMModel3 = fitgmdist(data3,i,'RegularizationValue',0.0001);

mu3=GMModel3.mu';
mu3=R73*mu3;
mu3=mu3';
mu3=mu3+T3;
Sigma3=GMModel3.Sigma;

for i1=1:i
Sigma3(:,:,i1)=R73*Sigma3(:,:,i1)*R73';
end
GMModel4 = fitgmdist(data4,i,'RegularizationValue',0.0001);

mu4=GMModel4.mu';
mu4=R74*mu4;
mu4=mu4';
mu4=mu4+T4;
Sigma4=GMModel4.Sigma;

for i1=1:i
Sigma4(:,:,i1)=R74*Sigma4(:,:,i1)*R74';
end
GMModel5 = fitgmdist(data5,i,'RegularizationValue',0.0001);

mu5=GMModel5.mu';
mu5=R75*mu5;
mu5=mu5';
mu5=mu5+T5;
Sigma5=GMModel5.Sigma;

for i1=1:i
Sigma5(:,:,i1)=R75*Sigma5(:,:,i1)*R75';
end
GMModel6 = fitgmdist(data6,i,'RegularizationValue',0.0001);

mu6=GMModel6.mu';
mu6=R76*mu6;
mu6=mu6';
mu6=mu6+T6;
Sigma6=GMModel6.Sigma;

for i1=1:i
Sigma6(:,:,i1)=R76*Sigma6(:,:,i1)*R76';
end










GMModel7 = fitgmdist(data7,i,'RegularizationValue',0.0001);

mu7=GMModel7.mu;
Sigma7=GMModel7.Sigma;


mu=[mu7;mu6;mu5;mu4;mu3;mu2;mu1];
Sigma=cat(3,Sigma7,Sigma6,Sigma5,Sigma4,Sigma3,Sigma2,Sigma1);
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



figure(11)
 hold on
  h1 = scatter3(data1(:,1),data1(:,2),data1(:,3),'.','b');
figure(12)
hold on
  h2 = scatter3(data2(:,1),data2(:,2),data2(:,3),'.','b');
  figure(13)
hold on
  h3 = scatter3(data3(:,1),data3(:,2),data3(:,3),'.','b');
  figure(14)
hold on
  h4 = scatter3(data4(:,1),data4(:,2),data4(:,3),'.','b');
  figure(15)
hold on
  h5 = scatter3(data5(:,1),data5(:,2),data5(:,3),'.','b');
  figure(16)
hold on
  h6 = scatter3(data6(:,1),data6(:,2),data6(:,3),'.','b');

  

  
  
  for n=1:i
Mu =mu1(n,1:3);
SIGMA=Sigma1(1:3,1:3,n);
[V,D]=eig(SIGMA);
r =2;
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
figure(5)
hold on
h5=surf(xsect1,ysect1,zsect1);
alpha(0.3)
xsect=[];
xsect1=[];
ysect=[];
ysect1=[];
zsect=[];
zsect1=[];

title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
  end

  
  
  for n=1:i
Mu =mu2(n,1:3);
SIGMA=Sigma2(1:3,1:3,n);
[V,D]=eig(SIGMA);
r =2;
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
figure(4)
hold on
h4=surf(xsect1,ysect1,zsect1);
alpha(0.3)
xsect=[];
xsect1=[];
ysect=[];
ysect1=[];
zsect=[];
zsect1=[];

title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
end
  
  
for n=1:i*7
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




hold off