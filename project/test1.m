Mu = [2 3 5]';
Sigma = [0.9 0.4 0.5;0.4 0.2 0.2 ;0.5 0.1 0.8];
r =2;
Nint=9;
[V,D] = eig(Sigma);



xhalf = linspace(sqrt(r^2*D(1,1)),0,Nint);
Ninthalf = round(Nint/2);
zsect = zeros(Nint,Ninthalf);
ysect = zeros(Nint,Ninthalf);
for ti = 1:Nint
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
surf(xsect1,ysect1,zsect1);
alpha(0.3)
xsect=[];

