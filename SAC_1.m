%% 模型简述
% SAC模型类似于斯坦福模型（SWM）,是一个连续模拟模型
% SAC模型是集中参数型的连续运算的确定性流域水文模型，模型参数则是根据流域特性、
% 降雨量和流量资料推求。
% 模型基本结构见《水文预报》（第五版）P170
% 本课程设计给出流域特性、降雨量和蒸发量等相关资料，时段为6hr
clc;
clear;
digits(6);
%% 流域资料导入与相关参数值设定
% 降雨量导入，单位为毫米
data_P = linspace(0,0,20);
data_P(1:4)=[40 80 100 30];
% 水面蒸发量导入，单位为毫米
data_E0=[0.5 0.8 1.2 0.6 0.9 1.5 2.4 1.2 1.2 2.0 3.2 1.6 0.9 1.5 2.4 1.2 1.5 2.0 3.5 1.5];
% 相关参数值设置
F = 100; % 流域面积为100平方公里
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

%% 模型主体
% 流域蒸发量计算
EP = KC*data_E0;
% 不透水产流计算
ROIMP = PCTIM*data_P;

%% 变动不透水面积上

% 上层蒸发和张力水计算
AUZTW = linspace(0,0,20);
AUZTW(1) = 20; % 上层张力水初始条件设定
AE1 = linspace(0,0,20);
AE1(1) = min(20,EP(1)*AUZTW(1)/20); 
for i = 2:20
    AUZTW(i)= min(20,AUZTW(i-1)-AE1(i-1)+data_P(i-1));
    AE1(i) = min(AUZTW(1),EP(i)*AUZTW(i)/20); 
end

% 有效降雨计算
PAV = max(0,data_P-20+AUZTW-AE1);

% 下层蒸发和张力水计算
ALZTW = linspace(0,0,20);
ALZTW(1) = 90; % 下层张力水初始条件设定
AE3 = linspace(0,0,20);
ADSUR = linspace(0,0,20);
AE3(1) = (EP(1)-AE1(1))*ALZTW(1)/(AUZTW(1)+ALZTW(1));
ADSUR(1) = PAV(1)*(ALZTW(1)-AE3(1))/90;
for i = 2:20
    ALZTW(i)= min(90,PAV(i-1)-ADSUR(i-1)+ALZTW(i-1)-AE3(i-1));
    AE3(i) = (EP(i)-AE1(i))*ALZTW(i)/(AUZTW(1)+ALZTW(1));
    ADSUR(i) = PAV(i)*(ALZTW(i)-AE3(i))/90;
end

% 地面径流计算
ARS = max(0,PAV-ADSUR+ALZTW-AE3-90);

%% 蒸发计算
% 耦合计算板块一
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
    %%% 上层张力水蓄量计算
    UT(i-1) = min(20,UZTW(i-1)-E1(i-1)+data_P(i-1));
    %%% 上层自由水蓄量
    UF(i-1) = min(30,data_P(i-1)+UZTW(i-1)+UZFW(i-1)-E1(i-1)-E2(i-1)-UT(i-1));
    %%% 更新后的下层张力水蓄量计算
    UF1(i-1)=0.6*UF(i-1);
    %%% 上层张力水和上层自由水计算
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
    % 上张力蒸发
    E1(i) = min(UZTW(i),EP(i)*UZTW(i)/20);
    % 上自由蒸发
    E2(i) = min(UZFW(i),EP(i)-E1(i));
end
UT(20) = min(20,UZTW(20)-E1(20)+data_P(20));

% 耦合运算板块2
RS= linspace(0,0,20);
RI= linspace(0,0,20);
for i=1:20
    %%% 地面径流计算
    RS(i)=max(0,data_P(i)+UZTW(i)+UZFW(i)-E1(i)-E2(i)-UZTWM-UZFWM);
    %%% 壤中流计算
    RI(i)=UF(i)*UZK;
end

% 耦合运算板块3
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
    %%% 下张力蒸发计算
    E3(i) = (EP(i)-E1(i)-E2(i))*(LZTW(i)/(UZTWM+LZTWM));
    %%% 下层张力水蓄量
    LT(i)=LZTW(i)-E3(i);
    %%% 中间值计算
    DEFR(i) = 1-(LZFS(i)+LZFP(i)+LT(i))/(LZFSM+LZFPM+LZTWM);
    PERC(i) = PBASE(i)*(1+ZPERC*DEFR(i)^2)*UF1(i)/30;
    %%% 下渗量的计算
    RATE(i) = min(PERC(i),LZFSM+LZFPM+LZTWM-LZFS(i)-LZFP(i)-LT(i));
    %%% 分配给下层自由水计算
    FX(i) = min(LZFSM+LZFPM-LZFS(i)-LZFP(i),max(RATE(i)-LZTWM+LT(i),RATE(i)*PFREE));
    %%% 分配给慢速自由水
    COEF(i) = (LZFPM/(LZFSM+LZFPM))*(2*(1-(LZFP(i)/LZFPM))/((1-LZFP(i)/LZFPM)+1-LZFS(i)/LZFSM));
    PERCP(i) = min(LZFPM-LZFP(i),max(FX(i)-LZFSM+LZFS(i),COEF(i)*FX(i)));
    %%% 下层快速水计算
    PERCS(i)=FX(i)-PERCP(i);
    LS(i)=PERCS(i)+LZFS(i);
    %%% 下层慢速水的计算
    LP(i)= PERCP(i)+ LZFP(i);
    %%% 快速地下水的计算
    RGS(i)=LS(i)*LZSK;
    %%% 慢速地下水的计算
    RGP(i)=LP(i)*LZPK;
    %%% 更新后的下层快速水
    LS1(i) = LS(i)-RGS(i);
    %%% 更新后的下层慢速水
    LP1(i) = LP(i)-RGP(i);
    %%% 分配给下层张力水计算
    PERCT(i) = RATE(i)-FX(i);
    %%% 更新下层张力水蓄量
    LT1(i) = LT(i)+PERCT(i);
    %%% 中间参数计算
    RATIO(i) = (LS1(i)+LP1(i)-36+LT1(i))/(20+100-36+90);   
    %%% 下层张力水计算
    if(LT1(i)/90<RATIO(i))
        LZTW(i+1)=90*RATIO(i);
    elseif(LT1(i)/90>=RATIO(i))
        LZTW(i+1)=LT1(i);
    end
    %%% 调整量计算
    DEL(i)=LZFS(i)-LT1(i);
    %%% 下层快自由水计算
    if(LT1(i)/90<RATIO(i))
        LZFS(i+1)=LS1(i)-DEL(i);
    else
        LZFS(i+1) = LS1(i);
    end
    
    %%% 下层慢自由水计算
    if(LT1(i)/90<RATIO(i))
        LZFP(i+1)=LP1(i)-max(0,DEL(i)-LS1(i));
    else
        LZFP(i+1)=LP1(i);
    end
end

%%% 河网蒸发计算
E4=linspace(0,0,20);
for i=1:20
    E4(i)=RIVA*EP(i);
end


%% 汇流计算
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
U0=(1-0.01-0.01)*100/(3.6*6); % U0为单位转换系数
% 壤中出流计算
for i=2:20
    QI(i)=CI*QI(i-1)+(1-CI)*RI(i-1)*U0;
end

% 快速出流计算
for i=2:20
    
    % 快速出流计算
    QS(i)=CGS*QS(i-1)+0.2*RGS(i-1)*U0; 
    % 慢速出流计算
    QP(i)=CGP*QP(i-1)+(1-CGP)*RGP(i-1)*U0; 
    % 河网总入流
    QT(i)=QI(i)+QS(i)+QP(i)+((ROIMP(i-1)+ADSUR(i-1)+ARS(i-1))*ADIMP+RS(i-1)*PAREA)*100/(3.6*6);
    % 出口流量计算
    Q(i)=CR*Q(i-1)+(1-CR)*0.5*(QT(i-1)+(QT(i)));
end

% 模型输出到excel表格
xlswrite('basic table.xlsx',Q');

