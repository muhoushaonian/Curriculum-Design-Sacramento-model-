%function DC=calcudc(canshu)
function DC = calcudc(canshu)
c=0;dc=0;dz=0;
SZ=365;
sacOut = SAC_4(canshu);
data=xlsread('data.xlsx');
data_Q=data(:,3);
for i=265:365
    dz=dz+data_Q(i);
    c=c+(data_Q(i)-sacOut(i))^2;
    dc=dc+(data_Q(i)-mean(data_Q))^2;
    % c = c+abs(data_Q(i)-sacOut(i));
end
DC =1/(1-(c/dc));
% DC = sum(sacOut(i)-data_Q(i));
% DC1 = sum(sqrt(mean((sacOut-data_Q).^2)))/365;

end