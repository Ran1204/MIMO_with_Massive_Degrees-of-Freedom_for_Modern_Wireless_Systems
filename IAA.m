Nt=4;
Nr=4;
Ki=2;
T_c = 1*10^-4;
c = 3*10^8;
code=randi([1,4],1,L);
SNR=50;

%steering vector of target 
theta1=[150,30]*pi/180;
theta2=[30,60]*pi/180;
theta3=[90,120]*pi/180;
r=0:Nr-1;r=2*r';
t=0:Nt-1;t=2*t';

Sik=exp(-1j*r*cos(theta1(1)));
Sik_hat=exp(-1j*t*cos(theta1(2)));
y1=kron(conj(Sik_hat),Sik);%
y1=awgn(y1,SNR,'measured');

Sik=exp(-1j*r*cos(theta2(1)));
Sik_hat=exp(-1j*t*cos(theta2(2)));
y2=kron(conj(Sik_hat),Sik);%
y2=awgn(y2,SNR,'measured');

Sik=exp(-1j*r*cos(theta3(1)));
Sik_hat=exp(-1j*t*cos(theta3(2)));
y3=kron(conj(Sik_hat),Sik);%
y3=awgn(y3,SNR,'measured');
yt=y1+y2+y3;


M=Nr*Nt;
m=0:M-1;
seta_sweep=0:1:180;n=1;
seta_sweep=seta_sweep*pi./180;
R_inv=diag(ones(1,M));%Initialise the array covariance matrix
K=length(seta_sweep)^2;

% Iteratively solve the array covariance matrix
A_seta=zeros(M,K);P=[];
n=1;
for i=seta_sweep
    for j=seta_sweep
        Sik=exp(-1j*r*cos(i));
        Sik_hat=exp(-1j*t*cos(j));
        hik=kron(conj(Sik_hat),Sik);
        A_seta(:,n)=hik;
        P(n)=A_seta(:,n)'*yt/(A_seta(:,n)'*A_seta(:,n));%Signal power initialisation 
        n=n+1;
    end
end



%Update parameters
for KK=1:5%¸üÐÂ
    R=A_seta*diag(abs(P).^2)*A_seta'+diag(ones(1,M));%Update the covariance matrix 
    R_inv=inv(R);
    for i=1:K   
        P(i)=A_seta(:,i)'*R_inv*yt/(A_seta(:,i)'*R_inv*A_seta(:,i));%update sn
    end

    
end

P_A=[];n=1;
for i=1:length(seta_sweep)
    for j=1:length(seta_sweep)
       P_A(i,j)=P(n);
        n=n+1;
    end
end
 mesh(seta_sweep*180./pi,seta_sweep*180./pi,abs(P_A))
 xlabel('DOA(deg)');ylabel('DOD(deg)');zlabel('Cost Function')