function truth= gen_truth_line(model,time,xstart, tbirth,tdeath)

%variables
truth.K= time;                   %length of data/number of scans数据长度/帧数
truth.X= cell(truth.K,1);             %ground truth for states of targets 目标状态真值   
truth.N= zeros(truth.K,1);            %ground truth for number of targets目标数量真值
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)目标标签真值
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)目标标识的绝对索引
truth.total_tracks= 0;          %total number of appearing tracks出现的航迹总数

%target initial states and birth/death times目标初始状态和目标新生/死亡时间
nbirths= length(tbirth);%目标航迹数量

% xstart(:,1)  = [ 0; 0; 0; -10 ];            tbirth(1)  = 1;     tdeath(1)  = 70;
% xstart(:,2)  = [ 400; -10; -600; 5 ];       tbirth(2)  = 1;     tdeath(2)  = truth.K+1;
% xstart(:,3)  = [ -800; 20; -200; -5 ];      tbirth(3)  = 1;     tdeath(3)  = 70;
% 
% xstart(:,4)  = [ 400; -7; -600; -4 ];       tbirth(4)  = 20;    tdeath(4)  = truth.K+1;
% xstart(:,5)  = [ 400; -2.5; -600; 10 ];     tbirth(5)  = 20;    tdeath(5)  = truth.K+1;
% xstart(:,6)  = [ 0; 7.5; 0; -5 ];           tbirth(6)  = 20;    tdeath(6)  = truth.K+1;
% 
% xstart(:,7)  = [ -800; 12; -200; 7 ];       tbirth(7)  = 40;    tdeath(7)  = truth.K+1;
% xstart(:,8)  = [ -200; 15; 800; -10 ];      tbirth(8)  = 40;    tdeath(8)  = truth.K+1;
% 
% xstart(:,9)  = [ -800; 3; -200; 15 ];       tbirth(9)   = 60;   tdeath(9)  = truth.K+1;
% xstart(:,10)  = [ -200; -3; 800; -15 ];     tbirth(10)  = 60;   tdeath(10) = truth.K+1;
% 
% xstart(:,11)  = [ 0; -20; 0; -15 ];         tbirth(11)  = 80;   tdeath(11) = truth.K+1;
% xstart(:,12)  = [ -200; 15; 800; -5 ];      tbirth(12)  = 80;   tdeath(12) = truth.K+1;

%generate the tracks生成航迹
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);%%当前目标状态
    for k=tbirth(targetnum):min(tdeath(targetnum),truth.K)%从新生时间到死亡时间循环
        targetstate = gen_newstate_fn(model,targetstate,'noiseless');
        truth.X{k}= [truth.X{k} targetstate];%产生的连续的目标状态真值存储在 truth.X{k}这个元胞数组中
        truth.track_list{k} = [truth.track_list{k} targetnum];%保存目标的唯一标识，即编号
        truth.N(k) = truth.N(k) + 1;%数量加一，即该时间有目标存在
     end
end
truth.total_tracks= nbirths;
