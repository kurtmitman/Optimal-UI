clear all
close all

rho = .9895;      % AR(1) paramter
sigma = .0034;    % standard deviation of shocks
lambda = 2;     % Tauchen coverage parameter
z_mu = 0;       % Mean shock value
znum = 41;      % Number of discrete tauchen values
[Z Pi] = tauchen(z_mu,sigma,rho,lambda,znum);
%[Z Pi] = tauchenhussey(znum,z_mu,rho,sigma,sigma/sqrt(1-rho^2));
P2 = (Pi'-eye(znum));
P2(znum,:)=ones(1,znum);
Pi0 = P2^(-1) * [zeros(znum-1,1); 1];
Z=exp(Z);
ynum=znum;
%Yinv=[zeros(1,(ynum-1)/2) 1 zeros(1,(ynum-1)/2)]';
%yEmis=eye(ynum+1);
%tPi=[Yinv'; Pi];
%yPi=[zeros(ynum+1,1), tPi];
%yseq=hmmgenerate(520000,yPi,yEmis)-1;

%save RandomShocksL2_Nov9 yseq
load RandomShocksL2_Nov9
options=optimset('Display','iter');

% 
% Dyn=0;
% if(Dyn==1)
%     load DynareZseq
%     simlen=length(z);
%     for t=1:simlen
%         index=1;
%         while exp(z(t))>Z(index) && index<znum
%             index=index+1;
%         end
%         yseq(t)=max(1,index-1);
%     
%     end
% end

%%
%x0=[0.2; .55; log(.5643); -5];
%[X val]=fsolve(@(x) calibrate_with_s(x,yseq,Pi,znum,Z),x0,options)

%x0=[0.5;0.1;0.55];
%x0=[0.438;0.079;0.55;1.2641573904;8.772922798];
x0=[0;0;0;0;0];
%x0=[-0.000971291317700;   1.135427590580511;  -0.220249013689227;  -0.071543710063760;   0.413569976903554];
%x0=[-0.002747543695163;   0.173171193282661;  -0.072157023846497;  -0.209910637839363;   1.586708582155921];
%x0=[ -0.004622416914275;  -0.696698192490107;   0.070179823465752;  -0.271176567065061;   3.103387437201057];
% x0=[
%   %-0.000410397667299
%    3.218980404394816
%   -0.660577508885336
%    0.405953172315601
%   -0.949005566151108];

x0=[
  %-0.000410397667299
   3.218980404394816
  -0.660577508885336
   0.405953172315601
  -0.949005566151108];

x0=[  -0.1073
   -0.3770
    0.0
   -1.5618];

x0 =[2.277657641398124
   1.379920912555895
  -0.210019471619202
   3.397420323966712
   0.0];

% x0 =[2.277657641398124
%    1.379920912555895
%   -0.210019471619202
%    3.397420323966712
%    0.0];
x0 =[   5.553163679747865
  -1.626140800294812
   0.967971786499170
   0.647774813939586
   0.279940270013144];




%[X val]=fsolve(@(x) calibrate_3step(x,yseq,Pi,znum,Z),x0,options)
calibrate_3step(x0,yseq,Pi,znum,Z)