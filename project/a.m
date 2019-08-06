close;
clear;
clc;
% mu1 = [1 2];
% Sigma1 = [2 0; 0 0.5];
% mu2 = [-3 -5];
% Sigma2 = [1 0;0 1];
% rng(1); % For reproducibility
% X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000)]; % 2000 x 2
X=csvread('C:\Users\EDA user\Desktop\pc\PointCloud11.csv');
GMModel = fitgmdist(X,1); %fit GMM distribution

figure
y = [zeros(10000,1);ones(10000,1)];
h = gscatter(X(1:20000,1),X(1:20000,2),y);
hold on
gmPDF = @(x1,x2)pdf(GMModel,[x1(:) x2(:)]);
g = gca;
fcontour(gmPDF)
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off
