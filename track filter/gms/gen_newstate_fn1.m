function X= gen_newstate_fn1(model,Xd,S,V,u,t)%增加输入S,u,t

%linear state space equation (CV model)

if ~isnumeric(V)
    if strcmp(V,'noise')
        V= model.sigma_V*model.B*randn(size(model.B,2),size(Xd,2));
    elseif strcmp(V,'noiseless')
        V= zeros(size(model.B,1),size(Xd,2));
    end
end

A_old=Xd(1,:);
Xd= Xd(2:5,:); 
if isempty(Xd)
    X= [];
else%此处有修改  
    X= model.F*Xd+ V;
%     A= A_old+S;
    X= [S;X;u;t];
end