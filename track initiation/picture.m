    subplot(2,1,1)
%     XX1=[98,90;97.2,94.8;92.4,95.2;97.6,96.8;96.4,98];
 XX1=[94.8,91.2;94.6,94.2;95.2,94.6;96.4,96.2;92.6,97.2];
    Y=[1,2,3,4,5];
MarkerType = 'o';
MarkerSize = 10;
LineWidth = 1.5;
LineStyle = '-';
% 绘制
st = stem(Y-0.2,XX1(:,1),...
    'MarkerEdgeColor','[0.54 0.27 0.07]',...      % 符号轮廓颜色
    'MarkerFaceColor','[1 0.89 0.7]',...      % 符号填充颜色
    'Marker',MarkerType,...       % 符号类型
    'MarkerSize',MarkerSize,...   % 符号尺寸
    'LineWidth',LineWidth,...     % 线宽
    'LineStyle',LineStyle,...     % 线型
    'Color','[0.54 0.27 0.07]');                  % 线的颜色
hold on
st = stem(Y+0.2,XX1(:,2),...
    'MarkerEdgeColor','[0.27 0.5 0.7]',...      % 符号轮廓颜色
    'MarkerFaceColor','[0.52 0.8 0.98]',...      % 符号填充颜色
    'Marker',MarkerType,...       % 符号类型
    'MarkerSize',MarkerSize,...   % 符号尺寸
    'LineWidth',LineWidth,...     % 线宽
    'LineStyle',LineStyle,...     % 线型
    'Color','[0.27 0.5 0.7]');                  % 线的颜色
   set(gca,'XTickLabel',{'','30','60','90','120','240'});
    for i = 1:5
    text(Y(i)-0.2,XX1(i,1)+1,num2str(XX1(i,1),'%g%%'),...
    'HorizontalAlignment','center',...
   'VerticalAlignment','bottom')
    text(Y(i)+0.2,XX1(i,2)+1,num2str(XX1(i,2),'%g%%'),...
        'HorizontalAlignment','center',...
   'VerticalAlignment','bottom')
  end
   legend('the proposed algorithm','the modified Hough transform');
hXLabel = xlabel('Clutter parameter');
hYLabel = ylabel('Track initiation success rate(%)'); 
title('The results of the track initiation success rate');
    axis([0 6 60 100])
    axis on
    grid on
    hold on
        
    subplot(2,1,2)
%     XX2=[1.6,2.59;5.07,9.54;7.23,16.78;9.96,35.64;48.94,83.02];
    XX2=[0.63,2.77;2.84,9.73;7.57,15.84;16.75,37.37;55.9,81.57];
st1 = stem(Y-0.2,XX2(:,1),...
    'MarkerEdgeColor','[0.54 0.27 0.07]',...      % 符号轮廓颜色
    'MarkerFaceColor','[1 0.89 0.7]',...      % 符号填充颜色
    'Marker',MarkerType,...       % 符号类型
    'MarkerSize',MarkerSize,...   % 符号尺寸
    'LineWidth',LineWidth,...     % 线宽
    'LineStyle',LineStyle,...     % 线型
    'Color','[0.54 0.27 0.07]');                  % 线的颜色
hold on
st1 = stem(Y+0.2,XX2(:,2),...
    'MarkerEdgeColor','[0.27 0.5 0.7]',...      % 符号轮廓颜色
    'MarkerFaceColor','[0.52 0.8 0.98]',...      % 符号填充颜色
    'Marker',MarkerType,...       % 符号类型
    'MarkerSize',MarkerSize,...   % 符号尺寸
    'LineWidth',LineWidth,...     % 线宽
    'LineStyle',LineStyle,...     % 线型
    'Color','[0.27 0.5 0.7]');                  % 线的颜色
   set(gca,'XTickLabel',{' ','30','60','90','120','240'});
    for i = 1:5
    text(Y(i)-0.2,XX2(i,1)+2.5,num2str(XX2(i,1),'%g%%'),...
    'HorizontalAlignment','center',...
   'VerticalAlignment','bottom')
    text(Y(i)+0.2,XX2(i,2)+2.5,num2str(XX2(i,2),'%g%%'),...
        'HorizontalAlignment','center',...
   'VerticalAlignment','bottom')
  end
   legend('the proposed algorithm','the modified Hough transform');
hXLabel = xlabel('Clutter parameter');
hYLabel = ylabel('False track occupancy rate(%)'); 
title('The results of the false track occupancy rate');
    axis([0 6 0 100])
    axis on
     grid on
    hold on
    
     pie3([9.413,48.22,5.50,16.29,1.34,19.22]);
legend('中立','偏向正面','完全正面','偏向负面','完全负面','与作品内容无关');

y=[98.2 96.4;97.8 96.2; 97.1 95.7; 97.6 96.9;96.4 93.8];
b=bar(y);
grid on;
ch = get(b,'children');
set(gca,'XTickLabel',{'30','60','90','120','240'});
b(2).FaceColor=[0.96 0.96 0.86];
b(1).FaceColor=[0.53 0.8 0.92];
% xtips1 = b.XEndPoints;
% ytips1 = b.YEndPoints; %获取Bar对象的XEndPoints和YEndPoints属性
% labels1 = string(b.YData); %获取条形末端坐标
% text(xtips1,ytips1,y,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
for i = 1:5
    text(i-0.2,y(i,1)+0.5,num2str(y(i,1),'%g%%'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
    text(i+0.3,y(i,2)+0.5,num2str(y(i,2),'%g%%'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
ylim([70 100]);
legend('the proposed algorithm','the modified logic algorithm');
ylabel('Trace detection rate(%)');
xlabel('Clutter parameters');