
A1=textread('C:\Users\EDA user\Desktop\GMM\21_pointcloud.txt');
A1(:,4)=[];


for i=1:2813
    data21(i,1:3)=A1(5*i,1:3);
end