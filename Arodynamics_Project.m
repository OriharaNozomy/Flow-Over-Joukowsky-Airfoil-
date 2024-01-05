%% Name and ID
%Name: Aya Abdelrahman Atya
%ID: 202000247

clear all

%% Givens
vfree = 100; % freestream velocity
c = 1; % c length
tm_c = 0.08; % max thickness to c
Cm_c = 0.04; % Max chamber to c
alpha = 7 * pi / 180; % angle of attack;

%% Airfoil parameters
b = c / 4;
e = tm_c / 1.3;
beta = 2 * Cm_c;
a = b * (1 + e) / cos(beta);
x0 = -b * e;
y0 = a * beta;

%% Drawing the airfoil
theta1 = linspace(0,2 * pi,50);
x1 = 2 * b * cos(theta1);
y1 = 2 * b * e .* (1 - cos(theta1)) .* sin(theta1) + 2 * b * beta .* sin(theta1).^2;
figure
plot(x1, y1,'r');
grid on;
xlim([-0.6 0.6]);
ylim([-0.25 0.25]);

%% Drawing the streamlines
rd=linspace(a,50,50);
thetad=linspace(0,2*pi,50);
Gamma=4*pi*vfree*a*sin(alpha+beta); %% vortex Strength Caculation
[Thetad,Rd]= meshgrid(thetad,rd);
Z=Rd.*exp(1i*Thetad)+x0+1i*y0; %z=z0+z'
Z1=Z+b^2./Z;%z1=z+b2/z
W = vfree .*exp(-1i.*alpha).*Z+vfree*a^2*exp( 1i *alpha)./Z+1i*log(Z)*Gamma/2/pi; %douplet+vortex+uniform over cylinder
psi=imag(W); %% Lines of stream line
psi1=psi+b^2./psi; %% cnvorting the stream lines from circular plan z to airfoil plan z1
figure
plot(x1, y1,'r');
hold on
contour(real(Z1),imag(Z1),psi1,'levelstep',2);
xlim([-0.6 0.6]);
ylim([-0.25 0.25]);

%% Drawing the velocity 
%Velocity over the upper surface
thetaup =linspace(0.0001,pi,50);
x1up = 2 * b * cos(thetaup);
y1up = 2 * b * e .* (1 - cos(thetaup)) .* sin(thetaup) + 2 * b * beta .* sin(thetaup).^2;
rup=b.*(1+e.*(1-cos(thetaup))+beta.*sin(thetaup));
xup=rup.*cos(thetaup);
yup=rup.*sin(thetaup);
xdup=xup-x0;
ydup=yup-y0;
thetaupdash=atan2(ydup,xdup);
vthetaupd = -vfree .* (2 .*sin(thetaupdash- alpha) + 2 .* sin(alpha + beta));
Vup=(vthetaupd.^2./(1-2*(b^2./rup.^2).*cos(2*thetaup)+(b^4./rup.^4))).^0.5;
figure
plot(x1up,Vup);
hold on

%Velocity over the lower surface
thetalow = linspace(pi,2*pi+0.0001,50);
x1low = 2 * b * cos(thetalow);
y1low = 2 * b * e .* (1 - cos(thetalow)) .* sin(thetalow) + 2 * b * beta .* sin(thetalow).^2;
rlow=b.*(1+e.*(1-cos(thetalow))+beta.*sin(thetalow));
xlow=rlow.*cos(thetalow);
ylow=rlow.*sin(thetalow);
xdlow=xlow-x0;
ydlow=ylow-y0;
thetalowdash=atan2(ydlow,xdlow);
vthetadlow = -vfree .* (2 .*sin(thetalowdash- alpha) + 2 .* sin(alpha + beta));
Vlow=(vthetadlow.^2./(1-2*(b^2./rlow.^2).*cos(2*thetalow)+(b^4./rlow.^4))).^0.5;
plot(x1low,Vlow);

%% Plotting pressure coeffecient

Cplow=1-(Vlow./vfree).^2; %Cp on the upper surface
Cpup=1-(Vup./vfree).^2; %Cp on the lower surface
figure
plot(x1up,Cpup);
hold on 
plot(x1low,Cplow);



%% Variation of angle of attack with lift coeffecient
AOA=-5:1:10;
Cl=2*pi*(1+e)*sin(AOA*pi/180+beta);
figure
plot(AOA,Cl);

%% Variation of moment coeffecient at c/4 with angle of attack

%M_c/4=My-Mx=integral(pu-pl)*(x1-c/4)*dx1-integral(pu-pl)*y1*dy1
%Write x1 and y1 in terms of theta
AOA=linspace(-5.*pi./180,10.*pi./180,50);

for i=1:length(AOA)
    vthetaupd = -vfree .* (2 .*sin(thetaupdash- AOA(i)) + 2 .* sin(AOA(i) + beta)); 
    Vup=(vthetaupd.^2./(1-2*(b^2./rup.^2).*cos(2*thetaup)+(b^4./rup.^4))).^0.5;

    vthetadlow = -vfree .* (2 .*sin(thetalowdash- AOA(i)) + 2 .* sin(AOA(i) + beta));
    Vlow=(vthetadlow.^2./(1-2*(b^2./rlow.^2).*cos(2*thetalow)+(b^4./rlow.^4))).^0.5;
    
    %I1=integral(pu*(x1-c/4)*dx1)-integral(pu*y1*dy1)
    I1=Vup.^2.*(2.*a.*cos(thetaup)-c/4).*(2.*a.*sin(thetaup))...
        -Vup.^2.*(2.*b.*e.*(1-cos(thetaup)).*sin(thetaup)+2.*b.*beta.*sin(thetaup).^2)...
        .*(2.*b.*e.*((1+sin(thetaup)).*sin(thetaup)+(1-cos(thetaup)).*cos(thetaup))...
        +4.*beta.*b.*sin(thetaup).*cos(thetaup));

    %I1=integral(-pu*(x1-c/4)*dx1)+integral(pu*y1*dy1)
    I2=-Vlow.^2.*(2.*a.*cos(thetalow)-c/4).*(2.*a.*sin(thetalow))+...
        Vlow.^2.*(2.*b.*e.*(1-cos(thetalow)).*sin(thetalow)+2.*b.*beta.*sin(thetalow).^2)...
        .*(2.*b.*e.*((1+sin(thetalow)).*sin(thetalow)+(1-cos(thetalow)).*cos(thetalow))...
        +4.*beta.*b.*sin(thetalow).*cos(thetalow));

    Q1(i)=trapz(I1,thetaup)/(0.5*1.225*vfree^2*c^2);
    Q2(i)=trapz(I2,thetalow)/(0.5*1.225*vfree^2*c^2);
    Cm(i)=(Q1(i)+Q2(i));

end
figure
plot(AOA,Cm);

%% Put the coordinates on text file for ANSYS

fileID = fopen('airfoil.txt', 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Unable to open file for writing');
end

% Write the result to the file
zcoordinate=zeros(size(x1));
[r,c]=size(x1);
group=ones(size(x1));

points=1:1:c;
fprintf(fileID, '%f\t%f\t%f\t%f\t%f\n', [group;points;x1; y1;zcoordinate]);
fprintf(fileID,'%f\t%f\t%f\t%f\t%f\n',[2;1;x1(c);y1(c);0]);
fprintf(fileID,'%f\t%f\t%f\t%f\t%f\n',[2;2;x1(1);y1(1);0]);
% Close the file
fclose(fileID);






