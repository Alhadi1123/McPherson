clc
clear
close all
s=2;
yy=0.5;
Rwheel=0.29;
 
ya0=0+yy;
za0=0;
yb0=0.2490+yy;  
zb0=-0.0608;
yc0=0.3721+yy;
zc0=0.0275;
yd0=0.1074+yy;  
zd0=0.5825;
 
za0=za0-(zc0-Rwheel);
zb0=zb0-(zc0-Rwheel);
zd0=zd0-(zc0-Rwheel);
zc0=zc0-(zc0-Rwheel);
 
theta0=atan((za0-zb0)/(yb0-ya0));
L01=sqrt((ya0-yb0)^2+(za0-zb0)^2);
L02=sqrt((yd0-yb0)^2+(zd0-zb0)^2);
L03=L02/3;
L05=0.2*L02;
L04=L02-L03-L05;
wheel_axis_dist=(zc0-zb0)/(zd0-zb0)*L02/L03;
 
ms=453;
mu=71;
ks=17658;
bs=1950;
Rwheel=0.29;
wheelWidth=Rwheel;
I=0.021;
kt1=183887;
bt1=2500;
kt2=50000;
bt2=2500;
 
k1=L02/(yb0*L02+yb0*zb0-yc0*zb0);
k2=k1*zb0/L02;
k3=k2*(zc0-zb0-L02);
k4=k3-Rwheel*k2;
k5=yd0*zb0-zd0*yb0;
k6s=ms+mu*k3+I*k2^2-((mu*k3+I*k2^2)^2)/(mu+mu*k3+I*k2^2);
k6u=mu+mu*k3+I*k2^2-((mu*k3+I*k2^2)^2)/(ms+mu*k3+I*k2^2);
k7s=(mu*k3+I*k2^2)/(mu+mu*k3+I*k2^2);
k7u=(mu*k3+I*k2^2)/(ms+mu*k3+I*k2^2);
 
tw=5;
h=0.001;

Zr(1:0.2/h)=0;
Zr(0.2/h+1:1/h+0.2/h+1)=0.2*sin(0:3.14*h/1:3.14);
Zr(end+1:tw/h+1)=0;

% Zr=0:0.4*h:0.4*2;             %choose between those for the input signals
% Zr(end+1:tw/h+1)=Zr(end);

% Zr(1)=0;
% Zr(end+1:tw/h+1)=0.2;
dZr=diff(Zr)/h;
 
 
t(1)=0;
Y(1,1)=0;   %zs
Y(2,1)=0;   %zu
Y(3,1)=0;   %s
Y(4,1)=0;   %u
 
n=4;
for i=1:tw/h
    zr=Zr(i);
    dzr=dZr(i);
    zs=Y(1,i);
    zu=Y(2,i);
    s=Y(3,i);
    u=Y(4,i);
    m1=[s,u,(1/k6s)*(((ks*k1*k5*(1-L02/(sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(1-k7s)-k7s*(kt1*(zu-zr)+bt1*(u-dzr))) ...
        1/k6u*(((ks*k1*k5*(1-(L02/sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(k7u-1)-(kt1*(zu-zr)+bt1*(u-dzr)))];
    
    zs=Y(1,i)+0.5*m1(1)*h;
    zu=Y(2,i)+0.5*m1(2)*h;
    s=Y(3,i)+0.5*m1(3)*h;
    u=Y(4,i)+0.5*m1(4)*h;
    m2=[s,u,(1/k6s)*(((ks*k1*k5*(1-L02/(sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(1-k7s)-k7s*(kt1*(zu-zr)+bt1*(u-dzr))) ...
        1/k6u*(((ks*k1*k5*(1-(L02/sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(k7u-1)-(kt1*(zu-zr)+bt1*(u-dzr)))];
 
    zs=Y(1,i)+0.5*m2(1)*h;
    zu=Y(2,i)+0.5*m2(2)*h;
    s=Y(3,i)+0.5*m2(3)*h;
    u=Y(4,i)+0.5*m2(4)*h;
    m3=[s,u,(1/k6s)*(((ks*k1*k5*(1-L02/(sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(1-k7s)-k7s*(kt1*(zu-zr)+bt1*(u-dzr))) ...
        1/k6u*(((ks*k1*k5*(1-(L02/sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(k7u-1)-(kt1*(zu-zr)+bt1*(u-dzr)))];
 
    zs=Y(1,i)+m3(1)*h;
    zu=Y(2,i)+m3(2)*h;
    s=Y(3,i)+m3(3)*h;
    u=Y(4,i)+m3(4)*h;
    m4=[s,u,(1/k6s)*(((ks*k1*k5*(1-L02/(sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(1-k7s)-k7s*(kt1*(zu-zr)+bt1*(u-dzr))) ...
        1/k6u*(((ks*k1*k5*(1-(L02/sqrt(2*k1*k5*(zu-zs)+L02^2))))+kt2*(k4)^2*(zu-zs)+bs*((k1*k5)^2)*(u-s)/(2*k1*k5*(zu-zs)+L02^2)+bt2*k4^2*(u-s))*(k7u-1)-(kt1*(zu-zr)+bt1*(u-dzr)))];
    Y(:,i+1)=Y(:,i)+h*(m1'+2*m2'+2*m3'+m4')/6;
    
    t(i+1)=t(i)+h;    
end
f1=figure('Name','System Response','NumberTitle','off');
plot(t,Y(1,:),'r','DisplayName','Zs');
grid
hold on 
plot(t,Y(2,:),'k','DisplayName','Zu');
plot(t,Zr,'b','DisplayName','Zr');
 
tt=linspace(0,2*pi,100);
sins=sin(tt);
coss=cos(tt);
 
Rpin=0.02;
Lthick=Rpin/2;
 
ya=ya0;
yd=yd0;
jj=0;
i=1;
f2=figure('Name','Half Car','NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
tic

myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 25;  %can adjust this, 5 - 10 works well for me
open(myVideo)

for i = 1:40:length(t)
    zs=Y(1,i);
    zu=Y(2,i);
    zr=Zr(i);
    theta=k1*(zu-zs);
    phi=k2*(zu-zs);
    
    za=za0+zs;
    zd=zd0+zs;
    yb=yb0*cos(theta)-zb0*sin(theta);
    zb=zb0*cos(theta)+yb0*sin(theta)+zs;
    yc=yb+(yc0-yb0)*cos(phi)+(zc0-zb0)*sin(phi);
    zc=zc0+zu;
    L1=sqrt((ya-yb)^2+(za-zb)^2);
    
    L2=sqrt((zd-zb)^2+(yd-yb)^2);
    L4=L2-L03-L05;
    psy=atan((zd-zb)/(yb-yd));
    ye=yd+L05*cos(psy);
    ze=zd-L05*sin(psy);
    yf=yb-L03*cos(psy);
    zf=zb+L03*sin(psy);
    Y_wheel_axis=yb-wheel_axis_dist*(yb-yf);
    Z_wheel_axis=zb+wheel_axis_dist*(zf-zb);
    theta=theta0-theta;
    YL1=[ya+Lthick*sin(theta) yb+Lthick*sin(theta) yb-Lthick*sin(theta) ya-Lthick*sin(theta) ya+Lthick*sin(theta)];
    ZL1=[za+Lthick*cos(theta) zb+Lthick*cos(theta) zb-Lthick*cos(theta) za-Lthick*cos(theta) za+Lthick*cos(theta)];
    YL2=[yb+Lthick*sin(psy) yf+Lthick*sin(psy) yf-Lthick*sin(psy) yb-Lthick*sin(psy) yb+Lthick*sin(psy)];
    ZL2=[zb+Lthick*cos(psy) zf+Lthick*cos(psy) zf-Lthick*cos(psy) zb-Lthick*cos(psy) zb+Lthick*cos(psy)];
    YL3=[yd+Lthick*sin(psy) ye+Lthick*sin(psy) ye-Lthick*sin(psy) yd-Lthick*sin(psy) yd+Lthick*sin(psy)];
    ZL3=[zd+Lthick*cos(psy) ze+Lthick*cos(psy) ze+Lthick*cos(psy) zd-Lthick*cos(psy) zd+Lthick*cos(psy)];
    YL4=[Y_wheel_axis-Lthick*cos(psy) yc+Lthick*sin(phi) yc-Lthick*sin(phi) Y_wheel_axis+Lthick*cos(psy) Y_wheel_axis-Lthick*cos(psy)];
    ZL4=[Z_wheel_axis+Lthick*sin(psy) zc+Lthick*cos(phi) zc-Lthick*cos(phi) Z_wheel_axis-Lthick*sin(psy) Z_wheel_axis+Lthick*sin(psy)];
    
    Ychassis=[ya+2*Rpin*sins(1:end/2) ya  yd-2*Rpin*coss(1:end/2) yc0 yc0 yc0+0.9*wheelWidth yc0+0.9*wheelWidth];
    Zchassis=[za-2*Rpin*coss(1:end/2) zd  zd-2*Rpin*sins(1:end/2) zd zd-(zd0-zc0)/4 zd-(zd0-zc0)/4 zd+0.75*Rwheel];
    Ychassis1=[Ychassis 0.2 0.2];
    Zchassis1=[Zchassis Zchassis(end) Zchassis(1)];
    Ychassis(end+1:2*end)=-flip(Ychassis);
    Zchassis(end+1:2*end)=flip(Zchassis);
    Ywheel=[yc+Rwheel*sin(phi) yc+wheelWidth*cos(phi)+Rwheel*sin(phi) yc+wheelWidth*cos(phi)-Rwheel*sin(phi) yc-Rwheel*sin(phi) yc+Rwheel*sin(phi)];
    Zwheel=[zc+Rwheel*cos(phi) zc+Rwheel*cos(phi)-wheelWidth*sin(phi) zc-Rwheel*cos(phi)-wheelWidth*sin(phi) zc-Rwheel*cos(phi) zc+Rwheel*cos(phi)];
    
    n=10;           %spring coils
    dy=L4/(n+1);    %single spring coil length
    dx=L04/2;       %strut width
    yy1=yf+dx/2*sin(pi-psy)+dy/2*cos(pi-psy);
    zz1=zf+dy/2*sin(pi-psy)-dx/2*cos(pi-psy);
    
    figure(f2)
    fill([yf yy1],[zf zz1],'k',[-yf -yy1],[zf zz1],'k','LineWidth',2)
    hold on
    for j=1:n
        fill([yy1 yy1+((-1)^j)*dx*sin(pi-psy)+dy*cos(pi-psy)],[zz1 zz1+dy*sin(pi-psy)-(-1)^j*dx*cos(pi-psy)],'r',...
            [-yy1 -(yy1+((-1)^j)*dx*sin(pi-psy)+dy*cos(pi-psy))],[zz1 zz1+dy*sin(pi-psy)-(-1)^j*dx*cos(pi-psy)],'r','LineWidth',2)
        yy1=yy1+(-1)^j*dx*sin(pi-psy)+dy*cos(pi-psy);
        zz1=zz1+dy*sin(pi-psy)-(-1)^j*dx*cos(pi-psy);
    end
    fill([yy1,ye],[zz1,ze],'k',[-yy1,-ye],[zz1,ze],'k','LineWidth',2) %end of strut
    
    fill([-4,4,4,-4,-4],[zr,zr,0,0,zr],'k','FaceAlpha',0.7,'LineWidth',0.1)    %road
    fill(YL1,ZL1,'k',YL2,ZL2,'k',YL3,ZL3,'k',YL4,ZL4,'k',Ywheel,Zwheel,'k',Ychassis,Zchassis,'',yb+Rpin*coss,zb+Rpin*sins,'g',...
        -YL1,ZL1,'k',-YL2,ZL2,'k',-YL3,ZL3,'k',-YL4,ZL4,'k',-Ywheel,Zwheel,'k',-yb+Rpin*coss,zb+Rpin*sins,'g',...
        yd+Rpin*coss,zd+Rpin*sins,'g',...
        ye+Rpin*coss,ze+Rpin*sins,'b',...
        yf+Rpin*coss,zf+Rpin*sins,'b',...
        ya+Rpin*coss,za+Rpin*sins,'g',...
        -yd+Rpin*coss,zd+Rpin*sins,'g',...
        -ye+Rpin*coss,ze+Rpin*sins,'b',...
        -yf+Rpin*coss,zf+Rpin*sins,'b',...
        -ya+Rpin*coss,za+Rpin*sins,'g','LineWidth',0.1)
    fill([Ychassis(105) Ychassis(105)*(1-1/3) Ychassis(105)*(1/3-1) -Ychassis(105)],...
        [Zchassis(105) Zchassis(105)+0.5 (Zchassis(105)+0.5) Zchassis(105)],'b','LineWidth',0.1,'FaceAlpha',0.3)
    axis([-1.5 1.5 0 3])
    axis equal
    drawnow
    hold off
    

    frame = getframe(f2);
    writeVideo(myVideo, frame);

end
close(myVideo)
