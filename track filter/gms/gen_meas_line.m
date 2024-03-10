function meas= gen_meas_line(model,truth)

%variables
meas.K= truth.K;%���ݳ���
meas.Z= cell(truth.K,1);%������������
meas.clutter = zeros(truth.K,1) ; %�����Ӳ�����������
meas.P_D = zeros(truth.K,1) ; %���������
%generate measurements��������
for k=1:truth.K
    if truth.N(k) > 0%��Ŀ����ֵ����
        idx= rand(truth.N(k),1) <= model.P_D ;                                            %detected target indicesѰ�Ҹ�ʱ����ʵĿ��ı��
        meas.Z{k}= gen_observation_fn(model,truth.X{k}(:,idx),'noise');                          %single target observations if detected��̽�⵽�ĵ�Ŀ��۲����� 
    end
    N_c= poissrnd(model.lambda_c);                                                               %number of clutter points�Ӳ�������
    C= repmat(model.range_c(:,1),[1 N_c])+ diag(model.range_c*[ -1; 1 ])*rand(model.z_dim,N_c);  %clutter generation�����Ӳ�
    meas.Z{k}= [ meas.Z{k} C ];                                                                  %measurement is union of detections and clutter
    meas.clutter(k) = N_c ; %%����洢��ÿһʱ�̵��Ӳ�����
end
    