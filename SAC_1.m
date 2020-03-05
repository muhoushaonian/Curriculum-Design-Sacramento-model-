%% ģ�ͼ���
% SACģ��������˹̹��ģ�ͣ�SWM��,��һ������ģ��ģ��
% SACģ���Ǽ��в����͵����������ȷ��������ˮ��ģ�ͣ�ģ�Ͳ������Ǹ����������ԡ�
% ��������������������
% ģ�ͻ����ṹ����ˮ��Ԥ����������棩P170
% ���γ���Ƹ����������ԡ�����������������������ϣ�ʱ��Ϊ6hr
clc;
clear;
digits(6);
%% �������ϵ�������ز���ֵ�趨
% ���������룬��λΪ����
data_P = linspace(0,0,20);
data_P(1:4)=[40 80 100 30];
% ˮ�����������룬��λΪ����
data_E0=[0.5 0.8 1.2 0.6 0.9 1.5 2.4 1.2 1.2 2.0 3.2 1.6 0.9 1.5 2.4 1.2 1.5 2.0 3.5 1.5];
% ��ز���ֵ����
F = 100; % �������Ϊ100ƽ������
KC = 0.9; % 
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
UZK = 0.4; %
CI = 0.5; %
CGS = 0.8; %
CGP = 0.95; %
PAREA = 0.98; %
CR = 0.4; %
ZPERC =15; %
PFREE = 0.3; %

%% ģ������
% ��������������
EP = KC*data_E0;
% ��͸ˮ��������
ROIMP = PCTIM*data_P;

%% �䶯��͸ˮ�����

% �ϲ�����������ˮ����
AUZTW = linspace(0,0,20);
AUZTW(1) = 20; % �ϲ�����ˮ��ʼ�����趨
AE1 = linspace(0,0,20);
AE1(1) = min(20,EP(1)*AUZTW(1)/20); 
for i = 2:20
    AUZTW(i)= min(20,AUZTW(i-1)-AE1(i-1)+data_P(i-1));
    AE1(i) = min(AUZTW(1),EP(i)*AUZTW(i)/20); 
end

% ��Ч�������
PAV = max(0,data_P-20+AUZTW-AE1);

% �²�����������ˮ����
ALZTW = linspace(0,0,20);
ALZTW(1) = 90; % �²�����ˮ��ʼ�����趨
AE3 = linspace(0,0,20);
ADSUR = linspace(0,0,20);
AE3(1) = (EP(1)-AE1(1))*ALZTW(1)/(AUZTW(1)+ALZTW(1));
ADSUR(1) = PAV(1)*(ALZTW(1)-AE3(1))/90;
for i = 2:20
    ALZTW(i)= min(90,PAV(i-1)-ADSUR(i-1)+ALZTW(i-1)-AE3(i-1));
    AE3(i) = (EP(i)-AE1(i))*ALZTW(i)/(AUZTW(1)+ALZTW(1));
    ADSUR(i) = PAV(i)*(ALZTW(i)-AE3(i))/90;
end

% ���澶������
ARS = max(0,PAV-ADSUR+ALZTW-AE3-90);

%% ��������
% ��ϼ�����һ
PBASE = linspace(6,6,20);
UZTW=linspace(0,0,20);
UZTW(1)=20;
UZFW=linspace(0,0,20);
UZFW(1)=20;
E1=linspace(0,0,20);
E1(1)=0.45;
E2=linspace(0,0,20);
E2(1)=0;
UT=linspace(0,0,20);
UT(1)=20;
UF=linspace(0,0,20);
UF(1)=30;
UF1=linspace(0,0,20);
UF(1)=18;
for i=2:20
    %%% �ϲ�����ˮ��������
    UT(i-1) = min(20,UZTW(i-1)-E1(i-1)+data_P(i-1));
    %%% �ϲ�����ˮ����
    UF(i-1) = min(30,data_P(i-1)+UZTW(i-1)+UZFW(i-1)-E1(i-1)-E2(i-1)-UT(i-1));
    %%% ���º���²�����ˮ��������
    UF1(i-1)=0.6*UF(i-1);
    %%% �ϲ�����ˮ���ϲ�����ˮ����
    if(UT(i-1)/20<UF1(i-1)/30)
        UZTW(i)=20*(UT(i-1)+UF1(i-1))/(20+30);
    else
        UZTW(i)=UT(i-1);
    end

    if(UT(i-1)/20<UF1(i-1)/30)
        UZFW(I)=30*(UT(i-1)+UF1(i-1));
    else
        UZFW(i)=UF1(i-1);
    end
    % ����������
    E1(i) = min(UZTW(i),EP(i)*UZTW(i)/20);
    % ����������
    E2(i) = min(UZFW(i),EP(i)-E1(i));
end
UT(20) = min(20,UZTW(20)-E1(20)+data_P(20));

% ���������2
RS= linspace(0,0,20);
RI= linspace(0,0,20);
for i=1:20
    %%% ���澶������
    RS(i)=max(0,data_P(i)+UZTW(i)+UZFW(i)-E1(i)-E2(i)-UZTWM-UZFWM);
    %%% ����������
    RI(i)=UF(i)*UZK;
end

% ���������3
LZTW=linspace(0,0,21);
LZTW(1)=90;
LZFS=linspace(0,0,21);
LZFS(1)=10;
LZFP=linspace(0,0,20);
LZFP(1)=80;
LT=linspace(0,0,20);
LT(1)=90;
LS=linspace(0,0,20);
LP=linspace(0,0,20);
LS1=linspace(0,0,20);
LT1=linspace(0,0,20);
% LT1(1)=90;
LP1=linspace(0,0,20);
DEFR=linspace(0,0,20);
PERC=linspace(0,0,20);
RATE=linspace(0,0,20);
PERCT=linspace(0,0,20);
FX=linspace(0,0,20);
E3=linspace(0,0,20);
RATIO=linspace(0,0,20);
PERCS=linspace(0,0,20);
PERCP=linspace(0,0,20);
COEF=linspace(0,0,20);
RGS=linspace(0,0,20);
RGP=linspace(0,0,20);
DEL = linspace(0,0,20);
for i=1:20    
    %%% ��������������
    E3(i) = (EP(i)-E1(i)-E2(i))*(LZTW(i)/(UZTWM+LZTWM));
    %%% �²�����ˮ����
    LT(i)=LZTW(i)-E3(i);
    %%% �м�ֵ����
    DEFR(i) = 1-(LZFS(i)+LZFP(i)+LT(i))/(LZFSM+LZFPM+LZTWM);
    PERC(i) = PBASE(i)*(1+ZPERC*DEFR(i)^2)*UF1(i)/30;
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
    RATIO(i) = (LS1(i)+LP1(i)-36+LT1(i))/(20+100-36+90);   
    %%% �²�����ˮ����
    if(LT1(i)/90<RATIO(i))
        LZTW(i+1)=90*RATIO(i);
    elseif(LT1(i)/90>=RATIO(i))
        LZTW(i+1)=LT1(i);
    end
    %%% ����������
    DEL(i)=LZFS(i)-LT1(i);
    %%% �²������ˮ����
    if(LT1(i)/90<RATIO(i))
        LZFS(i+1)=LS1(i)-DEL(i);
    else
        LZFS(i+1) = LS1(i);
    end
    
    %%% �²�������ˮ����
    if(LT1(i)/90<RATIO(i))
        LZFP(i+1)=LP1(i)-max(0,DEL(i)-LS1(i));
    else
        LZFP(i+1)=LP1(i);
    end
end

%%% ������������
E4=linspace(0,0,20);
for i=1:20
    E4(i)=RIVA*EP(i);
end


%% ��������
QI=linspace(0,0,20);
QI(1)=0;
QS=linspace(0,0,20);
QS(1)=0;
QP=linspace(0,0,20);
QP(1)=0;
QT=linspace(0,0,20);
QT(1)=0;
Q=linspace(0,0,20);
Q(1)=0;
U0=(1-0.01-0.01)*100/(3.6*6); % U0Ϊ��λת��ϵ��
% ���г�������
for i=2:20
    QI(i)=CI*QI(i-1)+(1-CI)*RI(i-1)*U0;
end

% ���ٳ�������
for i=2:20
    
    % ���ٳ�������
    QS(i)=CGS*QS(i-1)+0.2*RGS(i-1)*U0; 
    % ���ٳ�������
    QP(i)=CGP*QP(i-1)+(1-CGP)*RGP(i-1)*U0; 
    % ����������
    QT(i)=QI(i)+QS(i)+QP(i)+((ROIMP(i-1)+ADSUR(i-1)+ARS(i-1))*ADIMP+RS(i-1)*PAREA)*100/(3.6*6);
    % ������������
    Q(i)=CR*Q(i-1)+(1-CR)*0.5*(QT(i-1)+(QT(i)));
end

% ģ�������excel���
xlswrite('basic table.xlsx',Q');

