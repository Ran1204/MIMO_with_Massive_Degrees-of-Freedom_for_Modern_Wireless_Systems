close all;
clear all
clc
Nt=4;
Nr=4;
Nc=31;
Ki=2;
T_c = 1*10^-4;
c = 3*10^8;
L = 200;
code=randi([1,4],1,L);
SNR=50;
%steering vector of target 1
theta=[30,150];
r=0:Nr-1;r=90*r';
Sik=exp(-1j*r*sin(theta(1)));
Sik_u=exp(-1j*r*sin(theta(2)));
hik=kron(conj(Sik_u),Sik);
for n=1:L
    xx(:,n)=hik*randi([1,8],1,1);%+0.05*randn(64,1);
     xx(:,n)=awgn(xx(:,n),SNR,'measured');
end
 xt{1}=xx;
% steering vector of target 2
theta=[60,30];
r=0:Nr-1;r=90*r';
Sik=exp(-1j*r*sin(theta(1)));
Sik_u=exp(-1j*r*sin(theta(2)));
hik=kron(conj(Sik_u),Sik);
for n=1:L
    xx(:,n)=hik*randi([1,8],1,1);%+0.05*randn(64,1);
    xx(:,n)=awgn(xx(:,n),SNR,'measured');
end
xt{2}=xx;
% steering vector of target 3
theta=[120,90];
r=0:Nr-1;r=90*r';
Sik=exp(-1j*r*sin(theta(1)));
Sik_u=exp(-1j*r*sin(theta(2)));
hik=kron(conj(Sik_u),Sik);
for n=1:L
    xx(:,n)=hik*randi([1,8],1,1);%+0.05*randn(64,1);
    xx(:,n)=awgn(xx(:,n),SNR,'measured');
end
xt{3}=xx;

J_up = zeros(1,2*Nc);
J_down = horzcat(eye(2*Nc-1),zeros(2*Nc-1,1));
J = vertcat(J_up,J_down);

C_i = randsrc(31,8,[-1,1]);
O = zeros(31,8);
temp = vertcat(C_i,O);
N_hat=Nt;

F = rand(2*Nc,N_hat); %doppler frequency (random from 0-1)
F_hat = cell(2*Nc,N_hat); %formula 9
self_plus_array = [0:1:2*Nc-1].';

for i=1:2*Nc
    for k=1:N_hat
        F_hat{i,k} = exp(1i*2*pi*F(i,k)*self_plus_array*T_c);
    end
end

tau = randsrc(7,1,[0,30]);

one_N_hat_T = ones(1,N_hat);
l_i_k = zeros(2*Nc,N_hat);


T_i_k = cell(2*Nc,N_hat);


N_n = randn(Nr,2*Nc);

X = zeros(Nr,2*Nc);
X_ISI = zeros(Nr,2*Nc);
X_MAI = zeros(Nr,2*Nc);


Ci=sign(randn(Nc,Nt));

J=[zeros(1,2*Nc-1),0;
   eye(2*Nc-1),zeros(2*Nc-1,1)];

li=[20,40];

fCi=sign(randn(Nc,Nt));

J=[zeros(1,2*Nc-1),0;
   eye(2*Nc-1),zeros(2*Nc-1,1)];

li=[20,40];

for m=1:Nt %the  channel model
    cim=Ci;
    cim(:,m)=[];
    
    JJ=[];
    for i=1:Ki
        JJ=[JJ,J.^li(i)];
    end
    Bim=kron(ones(1,Ki),J*[cim;zeros(Nc,Nt-1)]);
    Bi{m}=kron(ones(1,Ki),cim);
    PBim{m}=eye(2*Nc)-Bim*inv(Bim'*Bim+eye(6))*Bim';
    PBim_pro{m}=kron(eye(Nr),PBim{m});
end

PB1=zeros(2*Nr*Nc*Nt,2*Nr*Nc*Nt);
for i=1:Nt
    aa=(i-1)*2*Nc*Nr+1;
    PB1(aa:aa+2*Nc*Nr-1,aa:aa+2*Nc*Nr-1)=PBim_pro{i};
end


for ti=1:3
x=xt{ti};
Rx=x*x'./L;
Rx_inv=inv(Rx);
    
    
THETA1=0:1:180;
THETA2=0:1:180;
obj_fun=[];
n=1;m=1;
for theta1=THETA1
    n=1;
    for theta2=THETA2
        sik=exp(-1j*r*sin(theta1));
        sik_u=exp(-1j*r*sin(theta2));
        h{m,n}=kron(conj(sik_u),sik);
        obj_fun(m,n)=h{m,n}'*sum(x.').'/(h{m,n}'*Rx_inv*h{m,n});
        n=n+1;
    end
    m=m+1;
end

obj_t{ti}=abs(obj_fun);

end
mesh(THETA1,THETA2,obj_t{1}+obj_t{2}+obj_t{3})
xlabel('DOA(deg)');ylabel('DOD(deg)');zlabel('Cost Function')
