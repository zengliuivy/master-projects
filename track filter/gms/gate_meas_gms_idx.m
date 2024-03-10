function valid_idx= gate_meas_gms(z,gamma,model,m,P)

valid_idx = [];
zlength = size(z,2); if zlength==0, z_gate= []; return; end
plength = size(m,2);

for j=1:plength
    Sj= model.R + model.H*P(:,:,j)*model.H'; %S=R+HPH��
    Vs= chol(Sj); %����˹����ʽ�ֽ�
    det_Sj= prod(diag(Vs))^2; %����Ԫ�صĳ˻�
    inv_sqrt_Sj= inv(Vs);%�������
    iSj= inv_sqrt_Sj*inv_sqrt_Sj';
    nu= z- model.H*repmat(m(:,j),[1 zlength]);%������ظ�m()��J��zlength��
    dist= sum((inv_sqrt_Sj'*nu).^2);
    valid_idx= unique_faster([ valid_idx find( dist < gamma )]);
end
valid_idx=valid_idx(:)';%��󷵻ص�ʱ��ת��
