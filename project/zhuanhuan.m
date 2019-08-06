data1=textread('C:\Users\EDA user\Desktop\GMM\1_pointcloud.txt');
figure(1)
hold on
 scatter3(data1(:,1),data1(:,2),data1(:,3),'.','b');



data2=textread('C:\Users\EDA user\Desktop\GMM\21_pointcloud.txt');
 figure(2)
hold on
 scatter3(data2(:,1),data2(:,2),data2(:,3),'.','b');
data3=[];
data4=[];
data5=[];
data6=[];
data1(:,4)=data1(:,1);
data1(:,2)=-data1(:,2);
data1(:,3)=-data1(:,3);
data2(:,4)=data2(:,1);
data2(:,2)=-data2(:,2);
data2(:,3)=-data2(:,3);

data3(:,1)=data1(:,2);
data3(:,2)=data1(:,3);
data3(:,3)=data1(:,4);

data4(:,1)=data2(:,2);
data4(:,2)=data2(:,3);
data4(:,3)=data2(:,4);

pose1=textread('C:\Users\EDA user\Desktop\GMM\1_pose.txt');
pose2=textread('C:\Users\EDA user\Desktop\GMM\21_pose.txt');

q1=[ -0.0221145442033, -3.11558682488e-05 , -0.000371488717354 ,0.999755374059];
q2=[ -0.00114818308711, 0.643118032594 ,0.000218562284025 ,0.765766171921];



R1=quat2rot(q1);
R2=quat2rot(q2);

T1=[ 0.001457824721, 0, 4.2644405365];
T2=[ 2.81376481056, 0, 14.9369363785];



data3=R1*data3';
data3=data3';
data3=data3+T1;





data4=R2*data4';
data4=data4';
data4=data4+T2;


data5(:,1)=-data3(:,2);
data5(:,2)=-data3(:,3);
data5(:,3)=data3(:,1);


data6(:,1)=-data4(:,2);
data6(:,2)=-data4(:,3);
data6(:,3)=data4(:,1);


figure(3)
hold on
 scatter3(-data5(:,1),data5(:,2),data5(:,3),'.','b');

 scatter3(-data6(:,1),data6(:,2),data6(:,3),'.','r');
 
 figure(4)
hold on
 scatter3(-data6(:,1),data6(:,2),-data6(:,3),'.','r');
 