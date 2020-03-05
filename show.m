canshu = [0.898,0.497,0.433,0.589,0.708,0.994,0.989,0.332,12.596,0.28];
data=xlsread('data.xlsx');
QZ = data(:,3);
QB = data(:,1);
Q = (SAC_4(bestchrom))';
startDate = datenum('1997-01-01');
endDate = datenum('1997-12-31');
xDate=linspace(startDate,endDate,365);
figure;
[AX,H1,H2]=plotyy(xDate,[Q QZ],xDate,QB,'plot','bar');
%画双轴，AX(1)左轴，AX（2）右轴，H为曲线本身
hold on;
set(AX(2),'YDir','reverse','Ylim',[0,4*max(QB)],'YTick',[0:50:max(QB)],'FontSize',12); %设置右边轴为倒立
set(AX(1),'YLim',[0 2400],'YTick',[0:400:2400],'Fontsize',12,'YColor','k');
datetick('x','yyyymmdd','keeplimits');
set(gca,'XLim',[startDate endDate]);
set(get(AX(1),'Xlabel'),'String','时间(d)');
set(get(AX(1),'Ylabel'),'String','流量（m^3/s)');
set(get(AX(2),'Ylabel'),'string','降雨量/mm');
