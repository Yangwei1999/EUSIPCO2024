clear ;
clc;
% 初始settings
coeff =10;
N = 40 * coeff;
T = 100 * coeff;
% theta_true = [0,5*2*pi/N];
theta_true = [0,pi/4];
k = length(theta_true);
P = [1 0.4; 0.4 1];
% P =[1.9953,0.7981;0.7981,1.9953];
SNRList = 3 ;

ScanArea = [-pi/2 pi/2];
ScanPrec = 4000;

% 变量
VariableList = SNRList;
VariableLabel = 'SNR';

% legend 
ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};

%% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N,T,theta_true,P,VariableList(ii))];
end

% ArrayObject(1).sigma2 = 1;
% ArrayObject.P =    [1.9953,0.7981;0.7981,1.9953] ;
% ArrayObject.GenerateGuass();
nbLoop =100;
% 跟Loop 有关的变量  
ReceivedNum1 = 2;
DoA_Nb = zeros(ReceivedNum1,nbLoop,k);
MSE_Nb =  zeros(ReceivedNum1,nbLoop);
EiValue_Nb= zeros(ReceivedNum1,nbLoop,k);

% 跟自变量VariableList有关的变量 
ReceivedNum2 = ReceivedNum1;
MSE_VList = zeros(ReceivedNum2,length(VariableList));
Var_VList = zeros(ReceivedNum2,length(VariableList));
Bias_VList = zeros(ReceivedNum2,length(VariableList));
EiValue_VList = zeros(2,length(VariableList),k);  % Only ESPRIT-Type methods
CRB_Res = zeros(1,length(VariableList));

%%代码部分
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [DoA_Nb(1,Loop_i,:),MSE_Nb(1,Loop_i),EiValue_Nb(1,Loop_i,:)]  = ObjectNow.GetESPRIT();                
        [DoA_Nb(2,Loop_i,:),MSE_Nb(2,Loop_i),EiValue_Nb(2,Loop_i,:)]  = ObjectNow.GetGESPRIT('Empirical-2');   
    
    end

    for kk = 1:ReceivedNum2
        [MSE_VList(kk,object_i),Var_VList(kk,object_i),Bias_VList(kk,object_i)] = ObjectNow.GetStatNum(squeeze(DoA_Nb(kk,:,:)),MSE_Nb(kk,:));
    end
    
    for kk = 1:2
        EiValue_VList(kk,object_i,:) = mean(squeeze(EiValue_Nb(kk,:,:)),1);
    end
    CRB_Res(1,object_i) = trace(ObjectNow.GetCRB())/ObjectNow.k;
end

ObjectNow = ArrayObject(1);
P./ObjectNow.sigma2
%% 实验部分
% 实验DOA
% MSE
[MSE_VList(1,1),Var_VList(1,1),Bias_VList(1,1)];
[MSE_VList(2,1),Var_VList(2,1),Bias_VList(2,1)];
% 实验特征值
ESPRITEigenValue_E = squeeze(EiValue_VList(1,:,:)).';
GESPRITEigenValue_E = squeeze(EiValue_VList(2,:,:)).';
% ESPRITEigenValue_E = (ESPRITEiValue_Nb())
% GESPRITEigenValue_E= (mean(GESPRITEiValue_Nb,1))

%% 理论部分
% 获取对象信息
U_APA = ObjectNow.UsTrue;
g = (1- ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-2))./...
    (1 + ObjectNow.c .* (ObjectNow.EigsTrue./ObjectNow.sigma2).^(-1));
J_tmp = eye(ObjectNow.N);
n = ObjectNow.N-1;
J1 = J_tmp(1:n,:);
J2 = J_tmp(2:end,:);

% ESPRIT算法理论特征值
u1 = U_APA(:,1);
u2 = U_APA(:,2);
Alpha1 = g(1)  *  u1'*J1'*J2*u1 + g(2) * u2'*J1'*J2*u2;
Alpha2 = g(1)  *  g(2) *(n/N).^2 * exp(1i * theta_true(1)) * exp(1i * theta_true(2));
Alpha2 = g(1)  *  g(2) *(u1'*J1'*J2*u1 *u2'*J1'*J2*u2 - u1'*J1'*J2*u2 *u2'*J1'*J2*u1);
Delta = Alpha1^2 - 4 * Alpha2;
Lambda_Lit_ESPRIT = [(Alpha1 + sqrt(Delta))/2*(N/n)     (Alpha1 - sqrt(Delta))/2*(N/n)];
% GESPRIT算法理论特征值
Lambda_Lit_GESPRIT = (exp(1i*ObjectNow.ThetaTrue));


Angle_Lit_ESPRIT = angle(Lambda_Lit_ESPRIT);
Angle_Emp_ESPRIT = angle(ESPRITEigenValue_E);
Angle_Lit_GESPRIT = angle(Lambda_Lit_GESPRIT);
Angle_Emp_GESPRIT = angle(GESPRITEigenValue_E);

%% 绘图部分
figure;
hold on ;
subplot(2,1,1)
hold on;
%thoery
xline(Angle_Lit_ESPRIT(1),'LineWidth',1,'Color','#D95319','LineStyle','--')
xline(Angle_Lit_ESPRIT(2),'LineWidth',1,'Color','#D95319','LineStyle','--')
%Emperical
xline(Angle_Emp_ESPRIT(1),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
xline(Angle_Emp_ESPRIT(2),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
%true
xline(Angle_Lit_GESPRIT(1),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
xline(Angle_Lit_GESPRIT(2),'LineWidth',2,'Color','#A2142F','LineStyle','-')
legend('theory-1','theory-2','Emp-1','Emp-2','True-1','True-2')
% annotation('doublearrow',Angle_Lit_GESPRIT,[0.5,0.5])
title('Tradition ESPRIT')
axis([-0.2 1.2 0 2])

subplot(2,1,2)
hold on;
%thoery
xline(Angle_Lit_GESPRIT(1),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
xline(Angle_Lit_GESPRIT(2),'LineWidth',1.5,'Color','#A2142F','LineStyle','-')
%Emperical
xline(Angle_Emp_GESPRIT(1),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
xline(Angle_Emp_GESPRIT(2),'LineWidth',1,'Color','#0072BD','LineStyle','-.')
axis([-0.2 1.2 0 2])
legend('theory(True)-1','theory(True)-2','Emp-1','Emp-2')
title('Improved ESPRIT(GESPRIT)')

% miss out

