clear all;
clc;

dt=1;

m11=120000;m12=m11;m13=m11;%���Է�������1��2��3�Ҵ���x���ϵĹ��Է������
m21=217900;m22=m21;m23=m21;%y
m31=63600000;m32=m31;m33=m31;%z
d11=21500;d12=d11;d13=d11;%�����������1��2��3�Ҵ���x���ϵ�����������
d21=117000;d22=d21;d23=d21;%y
d31=8020000;d32=d31;d33=d31;%z

a1=m11/m21;b1=d21/m21;%�������ͱ���
a2=m12/m22;b2=d22/m22;
a3=m13/m23;b3=d23/m23;

%����������ŵ���������ɸ�kֵ

k1=0.03;k2=0.03;k3=0;k4=0;k5=0.3;k6=0.2;k7=0.05;%���ܲο���Ǭ˶ʿ����

%�캽�߳�ʼλ��

X1=zeros(1,101);Y1=zeros(1,101);
X1(1,1)=0;Y1(1,1)=0;
Phi1=zeros(1,101);Phi1(1,1)=0;
v1=0;

%������1��ʼλ�ü���ʼ״̬

X2=zeros(1,100);Y2=zeros(1,100);
X2(1,1)=0;Y2(1,1)=5;%��ʼ����
Phi2=zeros(1,101);Phi2(1,1)=0.1;%��ʼ�����
u2=zeros(1,101);r2=zeros(1,101);v2=zeros(1,101);
u2(1,1)=0;r2(1,1)=0;v2(1,1)=0;

%������2��ʼλ�ü���ʼ״̬

X3=zeros(1,100);Y3=zeros(1,100);
X3(1,1)=0;Y3(1,1)=-5;%��ʼ����
Phi3=zeros(1,101);Phi3(1,1)=0.1;%��ʼ�����
u3=zeros(1,101);r3=zeros(1,101);v3=zeros(1,101);
u3(1,1)=0;r3(1,1)=0;v3(1,1)=0;

%��Ӹ���״̬

L12d=5;L13d=5;%��������
Phi12d=pi/2;Phi13d=-pi/2;%������ԽǶ�

%��������

ePhi12=zeros(1,101);ePhi12=Phi1-Phi2;
lx2=zeros(1,101);ly2=zeros(1,101);

ePhi13=zeros(1,101);ePhi13=Phi1-Phi3;
lx3=zeros(1,101);ly3=zeros(1,101);

u2a=zeros(1,101);r2a=zeros(1,101);
du2a=zeros(1,101);dr2a=zeros(1,101);

u3a=zeros(1,101);r3a=zeros(1,101);
du3a=zeros(1,101);dr3a=zeros(1,101);

z12=zeros(1,101);z22=zeros(1,101);
ex2=zeros(1,101);ey2=zeros(1,101);ve2=zeros(1,101);

z13=zeros(1,101);z23=zeros(1,101);
ex3=zeros(1,101);ey3=zeros(1,101);ve3=zeros(1,101);

f12=zeros(1,101);f22=zeros(1,101);
f13=zeros(1,101);f23=zeros(1,101);

delta2=zeros(1,101);
delta3=zeros(1,101);

taou2=zeros(1,101);taor2=zeros(1,101);
taou3=zeros(1,101);taor3=zeros(1,101);

%��ʼѭ��

for t=2:100
    if t<30
     u1=1;    %leader�����ٶ� 
     r1(t)=0.1;  %leader�Ľ��ٶ�  
    end
    if t>=30&&t<80
     u1=1;    %leader�����ٶ� 
     r1(t)=0;    %leader�Ľ��ٶ� 
    end
    if t>=80&&t<101
     u1=1;    %leader�����ٶ� 
     r1(t)=0.05; %leader�Ľ��ٶ� 
    end
    
    %%%�캽����ѧģ��%%%
    
    Phi1(t)=dt*r1(t)+Phi1(t-1);
    X1(t)=dt*(u1*cos(Phi1(t-1))-v1*sin(Phi1(t-1)))+X1(t-1);
    Y1(t)=dt*(u1*sin(Phi1(t-1))+v1*cos(Phi1(t-1)))+Y1(t-1);
    x1(t)=X1(t);
    y1(t)=Y1(t);
    phi(t)=Phi1(t);
    
    %%%����leader-follower�Ĵ���ģ��%%%
    
    %%%%%%%%%%%������1%%%%%%%%%%%
    
    Phi2(t)=dt*r2(t-1)+Phi2(t-1);
    ePhi12(t)=dt*(r1(t-1)-r2(t-1))+ePhi12(t-1);%�캽�ߺ͸�����1�ĺ�������
    
    lx2(t)=-(X1(t-1)-X2(t-1))*cos(Phi2(t-1))-(Y1(t-1)-Y2(t-1))*sin(Phi2(t-1));%�캽�ߺ͸�����1����Ծ���
    ly2(t)=(X1(t-1)-X2(t-1))*sin(Phi2(t-1))-(Y1(t-1)-Y2(t-1))*cos(Phi2(t-1));
    
    theta12d=pi-Phi12d;
    
    lx2d=L12d*cos(theta12d);%�캽�ߺ͸�����1��������Ծ���
    ly2d=L12d*sin(theta12d);
    
    f12(t)=-u1+ly2d*r1(t-1);%ɾȥlx2d�ĵ���
    f22(t)=-v1-lx2d*r1(t-1);%ɾȥly2d�ĵ���
    ex2(t)=dt*(u2(t-1)*cos(ePhi12(t-1))+v2(t-1)*sin(ePhi12(t-1))+ey2(t-1)*r1(t-1)+f12(t-1))+ex2(t-1);
    ey2(t)=dt*(-u2(t-1)*sin(ePhi12(t-1))+v2(t-1)*cos(ePhi12(t-1))-ex2(t-1)*r1(t-1)+f22(t-1))+ey2(t-1);
    
    u2(t)=dt*((m22/m12)*v2(t-1)*r2(t-1)-(d12/m12)*u2(t-1)+(1/m12)*taou2(t-1))+u2(t-1);
    v2(t)=dt*(a2*u2(t-1)*r2(t-1)-b2*v2(t-1))+v2(t-1);
    r2(t)=dt*(((m12-m22)/m32)*u2(t-1)*v2(t-1)-(d32/m32)*r2(t-1)+(1/m32)*taor2(t-1))+r2(t-1);
    
    z12(t)=ex2(t)*cos(ePhi12(t))-ey2(t)*sin(ePhi12(t));%����任
    z22(t)=ex2(t)*sin(ePhi12(t))+ey2(t)*cos(ePhi12(t));
    
    u2a(t)=-k1*z12(t)-f12(t)*cos(ePhi12(t))+f22(t)*sin(ePhi12(t));%�������������������
    v2a(t)=-k2*z22(t)-f12(t)*sin(ePhi12(t))-f22(t)*cos(ePhi12(t));
    
    ue2(t)=u2(t)-u2a(t);%ǰ���ٶ�����Ư�ٶ�����ҡ���ٶ����
    ve2(t)=v2(t)-v2a(t);
    
    delta2(t)=(ve2(t)*(k4+k3*k5)*(a1-1)*(f12(t)+f22(t)))/((a1-1)*ve2(t)*(f12(t)+f22(t))+k3);
    r2a(t)=r1(t-1)+k4*sin(ePhi12(t))-k5*(a1-1)*ve2(t)*(f12(t)*cos(ePhi12(t))-f22(t)*sin(ePhi12(t)))+delta2(t);
    
    re2(t)=r2(t)-r2a(t);
    
    du2a(t)=(u2a(t)-u2a(t-1))/dt;
    dr2a(t)=(r2a(t)-r2a(t-1))/dt;
    taou2(t)=m12*(-k6*ue2(t-1)-(m22/m12)*v2(t-1)*r2(t-1)+(d12/m12)*u2(t-1)+du2a(t-1));%��������-ǰ������
    taor2(t)=m32*(-k7*re2(t-1)-((m12-m22)/m32)*u2(t-1)*v2(t-1)+(d32/m32)*r2(t-1)+dr2a(t-1));%��������-��ҡ������
    
    %λ�����
    X2(t)=dt*(u2(t-1)*cos(Phi2(t-1))-v2(t-1)*sin(Phi2(t-1)))+X2(t-1);
    Y2(t)=dt*(u2(t-1)*sin(Phi2(t-1))+v2(t-1)*cos(Phi2(t-1)))+Y2(t-1);
    
    %%%%%%%%%%%������3%%%%%%%%%%%
    
    Phi3(t)=dt*r3(t-1)+Phi3(t-1);
    ePhi13(t)=dt*(r1(t-1)-r3(t-1))+ePhi13(t-1);%�캽�ߺ͸�����2�ĺ�������
    
    lx3(t)=-(X1(t-1)-X3(t-1))*cos(Phi3(t-1))-(Y1(t-1)-Y3(t-1))*sin(Phi3(t-1));%�캽�ߺ͸�����2����Ծ���
    ly3(t)=(X1(t-1)-X3(t-1))*sin(Phi3(t-1))-(Y1(t-1)-Y3(t-1))*cos(Phi3(t-1));
    
    theta13d=pi-Phi13d;
    
    lx3d=L13d*cos(theta13d);%�캽�ߺ͸�����2��������Ծ���
    ly3d=L13d*sin(theta13d);

    f13(t)=-u1+ly3d*r1(t-1);%ɾȥlx3d�ĵ���
    f23(t)=-v1-lx3d*r1(t-1);%ɾȥlx3d�ĵ���
    ex3(t)=dt*(u3(t-1)*cos(ePhi13(t-1))+v3(t-1)*sin(ePhi13(t-1))+ey3(t-1)*r1(t-1)+f13(t))+ex3(t-1);
    ey3(t)=dt*(-u3(t-1)*sin(ePhi13(t-1))+v3(t-1)*cos(ePhi13(t-1))-ex3(t-1)*r1(t-1)+f23(t))+ey3(t-1);
    
    u3(t)=dt*((m23/m13)*v3(t-1)*r3(t-1)-(d13/m13)*u3(t-1)+(1/m13)*taou3(t-1))+u3(t-1);
    v3(t)=dt*(a3*u3(t-1)*r3(t-1)-b3*v3(t-1))+v3(t-1);
    r3(t)=dt*(((m13-m23)/m33)*u3(t-1)*v3(t-1)-(d33/m33)*r3(t-1)+(1/m33)*taor3(t-1))+r3(t-1);
    
    z13(t)=ex3(t)*cos(ePhi13(t))-ey3(t)*sin(ePhi13(t));%����任
    z23(t)=ex3(t)*sin(ePhi13(t))+ey3(t)*cos(ePhi13(t));
    
    
    u3a(t)=-k1*z13(t)-f13(t)*cos(ePhi13(t))+f23(t)*sin(ePhi13(t));%�������������������
    v3a(t)=-k2*z23(t)-f13(t)*sin(ePhi13(t))-f23(t)*cos(ePhi13(t));
    
    ue3(t)=u3(t)-u3a(t);%ǰ���ٶ�����Ư�ٶ�����ҡ���ٶ����
    ve3(t)=v3(t)-v3a(t);
    
    delta3(t)=(ve3(t)*(k4+k3*k5)*(a1-1)*(f13(t)+f23(t)))/((a1-1)*ve3(t)*(f13(t)+f23(t))+k3);
    r3a(t)=r1(t)+k4*sin(ePhi13(t))-k5*(a1-1)*ve3(t)*(f13(t)*cos(ePhi13(t))-f23(t)*sin(ePhi13(t)))+delta3(t);
    
    re3(t)=r3(t)-r3a(t);
    
    du3a(t)=(u3a(t)-u3a(t-1))/dt;
    dr3a(t)=(r3a(t)-r3a(t-1))/dt;
    taou3(t)=m13*(-k6*ue3(t-1)-(m23/m13)*v3(t-1)*r3(t-1)+(d13/m13)*u3(t-1)+du3a(t-1));%��������-ǰ������
    taor3(t)=m33*(-k7*re3(t-1)-((m13-m23)/m33)*u3(t-1)*v3(t-1)+(d33/m33)*r3(t-1)+dr3a(t-1));%��������-��ҡ������
    
    %λ�����
    X3(t)=dt*(u3(t-1)*cos(Phi3(t-1))-v3(t-1)*sin(Phi3(t-1)))+X3(t-1);
    Y3(t)=dt*(u3(t-1)*sin(Phi3(t-1))+v3(t-1)*cos(Phi3(t-1)))+Y3(t-1);
    
end

figure(1);
plot(x1,y1,'r-',X2,Y2,'b--',X3,Y3,'k-.');
%axis([-80,30,-15,50]);
legend('�캽���˶��켣','������1�˶��켣','������2�˶��켣');

figure(2);
plot(ePhi12,'r-');hold on;
plot(ePhi13,'b--');
axis([0,100,-1,1]);
legend('�캽�ߺ͸�����1�ĺ�������','�캽�ߺ͸�����2�ĺ�������');
