% plot 3D Fitted Gaussian Mixture 
close;
clear;
clc;
mu1 = [1 2 3];
Sigma1 = [1 0 0;0 1 0; 0 0 1];
mu2 = [-3 -5 -4];
Sigma2 = [1 0 0;0 1 0; 0 0 1];
rng(1); % For reproducibility
X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000)]; % 2000 x 2
GMModel = fitgmdist(X,2); %fit GMM distribution
figure
% y = [zeros(1000,1);ones(1000,1)];
h = scatter3(X(:,1),X(:,2),X(:,3),'.');
hold on
gmPDF = @(x1,x2,x3)reshape(pdf(GMModel,[x1(:) x2(:) x3(:)]),size(x1));
g = gca;
contour3(gmPDF)
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off
