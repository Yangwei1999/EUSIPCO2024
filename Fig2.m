clear ;
clc   ;

% 初始settings
coeff =[1 2 4 8 16];
N = 40 * 5;
T = 100 * 5;

% theta_true = [0,5*2*pi/N];
theta_true = [0,pi/4];
k = length(theta_true);
P = [1 0.4; 0.4 1];
SNR = 3;

ScanArea = [-pi/2 pi/2];
ScanPrec = 4000;

% Variable
VariableList = coeff;
VariableLabel = 'N';

% legend 
ShowLegend ={'Threshold','ESPRIT','GESPRIT','MUSIC','GMUSIC','CRB'};

% 实例化所有变量对象
ArrayObject = [];
for ii = 1:length(VariableList)
    ArrayObject = [ArrayObject ArraySignalModel(N*VariableList(ii),T*VariableList(ii),theta_true,P,SNR)];
end

nbLoop = 100;
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

% code part
for object_i = 1:length(ArrayObject)
    ObjectNow = ArrayObject(object_i);
    for Loop_i = 1: nbLoop  
        disp([num2str(object_i) '--' num2str(Loop_i)])
        ObjectNow.GenerateGuass();
        [DoA_Nb(1,Loop_i,:),MSE_Nb(1,Loop_i),EiValue_Nb(1,Loop_i,:)]  = ObjectNow.GetESPRIT();                
        [DoA_Nb(2,Loop_i,:),MSE_Nb(2,Loop_i),EiValue_Nb(2,Loop_i,:)]  = ObjectNow.GetGESPRIT('Empirical-2');   
%         [DoA_Nb(3,Loop_i,:),MSE_Nb(3,Loop_i)]                         = ObjectNow.GetMusic(ScanArea,ScanPrec);
%         [DoA_Nb(4,Loop_i,:),MSE_Nb(4,Loop_i)]                         = ObjectNow.GetGMusic(ScanArea,ScanPrec); 
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

figure;
subplot(1,2,1)
hold on ;
plot(VariableList*N,log10(MSE_VList(1,:)),'LineStyle','-','Color',	'#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList*N,log10(Var_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList*N,log10(Bias_VList(1,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
legend('MSE','Var','Bias')
title('ESPRIT')
xlabel(VariableLabel)
axis([min(VariableList*N) max(VariableList*N) -7 0])
%     
subplot(1,2,2)
hold on ;
plot(VariableList*N,log10(MSE_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','x','LineWidth',1.5)
plot(VariableList*N,log10(Var_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','o','LineWidth',1.5)
plot(VariableList*N,log10(Bias_VList(2,:)),'LineStyle','-','Color','#77AC30','Marker','*','LineWidth',1.5)
plot(VariableList*N,log10(CRB_Res(1,:)),'LineStyle','--','Color','#77AC30','Marker','o','LineWidth',1.5)

legend('MSE','Var','Bias','CRB')
title('GESPRIT')
xlabel(VariableLabel)
% axis([min(VariableList*N) max(VariableList*N) -7 0])





