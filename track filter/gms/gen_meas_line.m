function meas= gen_meas_line(model,truth)

%variables
meas.K= truth.K;%数据长度
meas.Z= cell(truth.K,1);%构建量测数组
meas.clutter = zeros(truth.K,1) ; %构建杂波的量测数组
meas.P_D = zeros(truth.K,1) ; %检测率数组
%generate measurements生成量测
for k=1:truth.K
    if truth.N(k) > 0%有目标真值存在
        idx= rand(truth.N(k),1) <= model.P_D ;                                            %detected target indices寻找该时刻真实目标的标号
        meas.Z{k}= gen_observation_fn(model,truth.X{k}(:,idx),'noise');                          %single target observations if detected被探测到的单目标观测数据 
    end
    N_c= poissrnd(model.lambda_c);                                                               %number of clutter points杂波的数量
    C= repmat(model.range_c(:,1),[1 N_c])+ diag(model.range_c*[ -1; 1 ])*rand(model.z_dim,N_c);  %clutter generation生成杂波
    meas.Z{k}= [ meas.Z{k} C ];                                                                  %measurement is union of detections and clutter
    meas.clutter(k) = N_c ; %%这里存储了每一时刻的杂波数量
end
    