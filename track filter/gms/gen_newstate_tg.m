function [X,g1,g2]= gen_newstate_tg(filter,model,X_old,u,v)%增加输入s,t,filter
%--- this new state generation function follows the linear state spc eqn.
% x= Ax_old + Bv

if isempty(X_old),
    X= [];
    g1=[];
    g2=[];
else
    mu= X_old(1,:);%取检测率
    mu=min(mu,0.999); mu=max(mu,0.001);%限定范围
    X_old(1,:)= mu;
%    u= X_old(6,:);v=X_old(7,:);
 %   sig= min(0.9*mu.*(1-mu),model.pdvarfac_tg);
   beta_avg_tmp= u./(u+v); 
   beta_avg_tmp=min(beta_avg_tmp,0.999); beta_avg_tmp=max(beta_avg_tmp,0.001);
%    if beta_avg_tmp~=mu
%        t=t
%    end
  % beta_avg_tmp= mu; 
%     beta_var_tmp=min(0.9*beta_avg_tmp.*(1-beta_avg_tmp),model.pdvarfac_tg);
%     beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./beta_var_tmp- 1;
   beta_var_tmp= (u.*v)./((u+v).^2.*(u+v+1)); 
   beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./min(filter.beta_factor.*beta_var_tmp,model.pdvarfac_tg)- 1;
    u= beta_tht_tmp.*beta_avg_tmp;  v= beta_tht_tmp.*(1-beta_avg_tmp);
 %   X= gen_newstate_fn_tg(model,X_old,betarnd((mu.*(1-mu)./sig -1).*mu,(mu.*(1-mu)./sig -1).*(1-mu))-mu, model.B*randn(size(model.B,2),size(X_old,2)));
  %此处改为CV模型
%  X= gen_newstate_fn(model,X_old, model.B*randn(size(model.B,2),size(X_old,2)));
 % X= gen_newstate_fn1(model,X_old,betarnd((mu.*(1-mu)./sig -1).*mu,(mu.*(1-mu)./sig -1).*(1-mu)), model.B*randn(size(model.B,2),size(X_old,2)));
  [r,g1,g2]=betarnd(u,v);
   X= gen_newstate_fn1(model,X_old,r, model.B*randn(size(model.B,2),size(X_old,2)),g1,g2);
end;