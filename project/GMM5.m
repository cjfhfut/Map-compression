filepath=['C:\Users\EDA user\Desktop\GMM1'];
i=20;
c=0;
MU=[];
SIGMA=[];
L=[];
figure(1);
hold on
for k=1:5:21
data = textread(...
        [filepath '\' num2str(k) '_pointcloud.txt']);    
    pose = textread(...
        [filepath '\' num2str(k) '_T.txt']);  
    R=pose(1:3,1:3);
    T=pose(1:3,4);
    
     data=(R'*data'-R'*T)';
    
    
 GMModel1 = fitgmdist(data,i,'RegularizationValue',0.0001);

mu1=GMModel1.mu;


Sigma1=GMModel1.Sigma;


MU=[MU;mu1];
SIGMA=cat(3,SIGMA,Sigma1);


for n=1:i
Mu =mu1(n,1:3);
Sigma=Sigma1(1:3,1:3,n);
[V,D]=eig(Sigma);
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
    
size1=size(SIGMA);
sizek=size1(1,3);
for v=1:sizek;
    zh(1:3,1:3,v)=SIGMA(:,:,v);
    zh(4,1:3,v)=MU(v,1:3);
end

C=nchoosek(1:sizek,2) ;
size2=size(C);
size3=size2(1,1);
for g=1:size3;
    
    kl = kldivGaussian(zh(4,1:3,C(g,1)), zh(1:3,1:3,C(g,1)),zh(4,1:3,C(g,2)), zh(1:3,1:3,C(g,2)));
    L(g,1)=kl;
    
    
    
end


    
end