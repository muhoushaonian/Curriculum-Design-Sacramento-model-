function sacOut = SAC_4(canshu)
% canshu: KC UZK PFREE AUZTW(1)
data=xlsread('data.xlsx');
data_P_real = data(:,1);
data_E0_real = data(:,2);
% ��زο�����ֵ����
F = 3893; % �о��������Ϊ3893ƽ������
% KC = 0.9;�ڱ������У�KC���ʶ�  canshu(1)
KC = canshu(1);
PCTIM = 0.01; %
ADIMP = 0.01; %
UZTWM = 20; %
UZFWM = 30; %
LZTWM = 90; %
LZFSM = 20; %
LZFPM = 100; %
RIVA = 0.01; %
LZSK = 0.2; %
LZPK = 0.02; %
% UZK = 0.4;�ڱ�������,UZK���ʶ� canshu(2)
UZK = canshu(2);
% CI = 0.5; �ڱ������У�CI���ʶ� canshu(4)
CI = canshu(4);
% CGS = 0.8; �ڱ������У�CGS���ʶ� canshu(5)
CGS = canshu(5);
% CGP = 0.95; �ڱ������У�CGP���ʶ� canshu(6)
CGP = canshu(6);
% PAREA = 0.98; �ڱ������У�PAREA���ʶ� canshu(7)
PAREA = canshu(7);
% CR = 0.4; �ڱ������У�CR���ʶ� canshu(8)
CR = canshu(8);
% ZPERC =15; �ڱ������У�ZPERC���ʶ� canshu(9)
ZPERC = canshu(9);
% PFREE = 0.3;�ڱ�������,PFREE���ʶ� canshu(3)
PFREE = canshu(3);
% RSERV = 0.3; �ڱ������У�RSERV���ʶ� canshu(10)
RSERV = canshu(10);
SAVED = RSERV*(LZFSM+LZFPM);

%% ģ������
% ��������������
EP = KC*data_E0_real;
% ��͸ˮ��������
ROIMP = PCTIM*data_P_real;

%% �䶯��͸ˮ�����

% �ϲ�����������ˮ����
AUZTW = linspace(0,0,365);
AUZTW(1) = 20; % �ϲ�����ˮ��ʼ�����趨
AE1 = linspace(0,0,365);
AE1(1) = min(20,EP(1)*AUZTW(1)/20); 
for i = 2:365
    AUZTW(i)= min(20,AUZTW(i-1)-AE1(i-1)+data_P_real(i-1));
    AE1(i) = min(AUZTW(1),EP(i)*AUZTW(i)/UZTWM); 
end

% ��Ч�������
PAV = max(0,data_P_real-UZTWM+AUZTW-AE1);

% �²�����������ˮ����
ALZTW = linspace(0,0,365);
ALZTW(1) = 90; % �²�����ˮ��ʼ�����趨
AE3 = linspace(0,0,365);
ADSUR = linspace(0,0,365);
AE3(1) = (EP(1)-AE1(1))*ALZTW(1)/(AUZTW(1)+ALZTW(1));
ADSUR(1) = PAV(1)*(ALZTW(1)-AE3(1))/LZTWM;
for i = 2:365
    ALZTW(i)= min(LZTWM,PAV(i-1)-ADSUR(i-1)+ALZTW(i-1)-AE3(i-1));
    AE3(i) = (EP(i)-AE1(i))*ALZTW(i)/(UZTWM+LZTWM);
    ADSUR(i) = PAV(i)*(ALZTW(i)-AE3(i))/LZTWM;
end

% ���澶������
ARS = max(0,PAV-ADSUR+ALZTW-AE3-LZTWM);

%% ��������
% ��ϼ�����һ

UZTW=linspace(0,0,365);
UZTW(1)=20;
UZFW=linspace(0,0,365);
UZFW(1)=20;
E1=linspace(0,0,365);
E1(1)=min(UZTW(1),EP(1)*(UZTW(1)/UZTWM));
E2=linspace(0,0,365);
E2(1)=0;
UT=linspace(0,0,365);
UT(1)=20;
UF=linspace(0,0,365);
UF(1)=30;
UF1=linspace(0,0,365);
UF1(1)=18;
for i=2:365
    %%% �ϲ�����ˮ��������
    UT(i-1) = min(UZTWM,UZTW(i-1)-E1(i-1)+data_P_real(i-1));
    %%% �ϲ�����ˮ����
    UF(i-1) = min(UZFWM,data_P_real(i-1)+UZTW(i-1)+UZFW(i-1)-E1(i-1)-E2(i-1)-UT(i-1));
    %%% ���º���²�����ˮ��������
    UF1(i-1)=(1-UZK)*UF(i-1);
    %%% �ϲ�����ˮ���ϲ�����ˮ����
    if(UT(i-1)/UZTWM<UF1(i-1)/UZFWM)
        UZTW(i)=UZTWM*(UT(i-1)+UF1(i-1))/(UZTWM+UZFWM);
    else
        UZTW(i)=UT(i-1);
    end

    if(UT(i-1)/UZTWM<UF1(i-1)/UZFWM)
        UZFW(I)=UZFWM*(UT(i-1)+UF1(i-1));
    else
        UZFW(i)=UF1(i-1);
    end
    % ����������
    E1(i) = min(UZTW(i),EP(i)*UZTW(i)/UZTWM);
    % ����������
    E2(i) = min(UZFW(i),EP(i)-E1(i));
end
UT(365) = min(UZTWM,UZTW(365)-E1(365)+data_P_real(365));

% ���������2
RS= linspace(0,0,365);
RI= linspace(0,0,365);
for i=1:365
    %%% ���澶������
    RS(i)=max(0,data_P_real(i)+UZTW(i)+UZFW(i)-E1(i)-E2(i)-UZTWM-UZFWM);
    %%% ����������
    RI(i)=UF(i)*UZK;
end

% ���������3
LZTW=linspace(0,0,366);
LZTW(1)=90;
LZFS=linspace(0,0,366);
LZFS(1)=10;
LZFP=linspace(0,0,365);
LZFP(1)=80;
LT=linspace(0,0,365);
LT(1)=90;
LS=linspace(0,0,365);
LP=linspace(0,0,365);
LS1=linspace(0,0,365);
LT1=linspace(0,0,365);
LP1=linspace(0,0,365);
DEFR=linspace(0,0,365);
PERC=linspace(0,0,365);
RATE=linspace(0,0,365);
PERCT=linspace(0,0,365);
FX=linspace(0,0,365);
E3=linspace(0,0,365);
RATIO=linspace(0,0,365);
PERCS=linspace(0,0,365);
PERCP=linspace(0,0,365);
COEF=linspace(0,0,365);
RGS=linspace(0,0,365);
RGP=linspace(0,0,365);
DEL = linspace(0,0,365);
PBASE=LZFSM*LZSK+LZFPM*LZPK;
for i=1:365    
    %%% ��������������
    E3(i) = (EP(i)-E1(i)-E2(i))*(LZTW(i)/(UZTWM+LZTWM));
    %%% �²�����ˮ����
    LT(i)=LZTW(i)-E3(i);
    %%% �м�ֵ����
    DEFR(i) = 1-(LZFS(i)+LZFP(i)+LT(i))/(LZFSM+LZFPM+LZTWM);
    PERC(i) = PBASE*(1+ZPERC*DEFR(i)^2)*UF1(i)/UZFWM;
    %%% �������ļ���
    RATE(i) = min(PERC(i),LZFSM+LZFPM+LZTWM-LZFS(i)-LZFP(i)-LT(i));
    %%% ������²�����ˮ����
    FX(i) = min(LZFSM+LZFPM-LZFS(i)-LZFP(i),max(RATE(i)-LZTWM+LT(i),RATE(i)*PFREE));
    %%% �������������ˮ
    COEF(i) = (LZFPM/(LZFSM+LZFPM))*(2*(1-(LZFP(i)/LZFPM))/((1-LZFP(i)/LZFPM)+1-LZFS(i)/LZFSM));
    PERCP(i) = min(LZFPM-LZFP(i),max(FX(i)-LZFSM+LZFS(i),COEF(i)*FX(i)));
    %%% �²����ˮ����
    PERCS(i)=FX(i)-PERCP(i);
    LS(i)=PERCS(i)+LZFS(i);
    %%% �²�����ˮ�ļ���
    LP(i)= PERCP(i)+ LZFP(i);
    %%% ���ٵ���ˮ�ļ���
    RGS(i)=LS(i)*LZSK;
    %%% ���ٵ���ˮ�ļ���
    RGP(i)=LP(i)*LZPK;
    %%% ���º���²����ˮ
    LS1(i) = LS(i)-RGS(i);
    %%% ���º���²�����ˮ
    LP1(i) = LP(i)-RGP(i);
    %%% ������²�����ˮ����
    PERCT(i) = RATE(i)-FX(i);
    %%% �����²�����ˮ����
    LT1(i) = LT(i)+PERCT(i);
    %%% �м��������
    RATIO(i) = (LS1(i)+LP1(i)-36+LT1(i))/(LZFSM+LZFPM-SAVED+LZTWM);   
    %%% �²�����ˮ����
    if(LT1(i)/LZTWM<RATIO(i))
        LZTW(i+1)=90*RATIO(i);
    elseif(LT1(i)/LZTWM>=RATIO(i))
        LZTW(i+1)=LT1(i);
    end
    %%% ����������
    DEL(i)=LZFS(i)-LT1(i);
    %%% �²������ˮ����
    if(LT1(i)/LZTWM<RATIO(i))
        LZFS(i+1)=LS1(i)-DEL(i);
    else
        LZFS(i+1) = LS1(i);
    end
    
    %%% �²�������ˮ����
    if(LT1(i)/LZTWM<RATIO(i))
        LZFP(i+1)=LP1(i)-max(0,DEL(i)-LS1(i));
    else
        LZFP(i+1)=LP1(i);
    end
end

%%% ������������
E4=linspace(0,0,365);
for i=1:365
    E4(i)=RIVA*EP(i);
end


%% ��������
QI=linspace(0,0,365);
QI(1)=0;
QS=linspace(0,0,365);
QS(1)=0;
QP=linspace(0,0,365);
QP(1)=0;
QT=linspace(0,0,365);
QT(1)=0;
Q=linspace(0,0,365);
Q(1)=0;
U0=(1-0.01-0.01)*F/(3.6*24); % U0Ϊ��λת��ϵ��
% ���г�������
for i=2:365
    QI(i)=CI*QI(i-1)+(1-CI)*RI(i-1)*U0;
end

% ���ٳ�������
for i=2:365
    
    % ���ٳ�������
    QS(i)=CGS*QS(i-1)+0.2*RGS(i-1)*U0; 
    % ���ٳ�������
    QP(i)=CGP*QP(i-1)+(1-CGP)*RGP(i-1)*U0; 
    % ����������
    QT(i)=QI(i)+QS(i)+QP(i)+((ROIMP(i-1)+ADSUR(i-1)+ARS(i-1))*ADIMP+RS(i-1)*PAREA)*100/(3.6*6);
    % ������������
    Q(i)=CR*Q(i-1)+(1-CR)*0.5*(QT(i-1)+(QT(i)));
end
sacOut = Q;
end
