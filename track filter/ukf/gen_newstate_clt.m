function [X,u_new,v_new]= gen_newstate_clt(filter,model,X_old,u,v)%接口改变
%--- this new state generation function follows the linear state spc eqn.
% x= Ax_old + Bv

if isempty(X_old),
    X= [];
    u_new=[];
    v_new=[];
else
    mu= X_old(1,:); mu=min(mu,0.999); mu=max(mu,0.001); X_old(1,:)= mu;
    sig= min(0.9*mu.*(1-mu),model.pdvarfac_clt);
    beta_avg_tmp= u./(u+v); 
    beta_avg_tmp=min(beta_avg_tmp,0.999); beta_avg_tmp=max(beta_avg_tmp,0.001);
    beta_var_tmp=min(0.9*beta_avg_tmp.*(1-beta_avg_tmp),model.pdvarfac_clt);
    beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./beta_var_tmp- 1;
%     beta_var_tmp= (u.*v)./((u+v).^2.*(u+v+1)); 
%     beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./min(filter.beta_factor.*beta_var_tmp,model.pdvarfac_clt)- 1;
    u= beta_tht_tmp.*beta_avg_tmp;  v= beta_tht_tmp.*(1-beta_avg_tmp);
   [d_new,u_new,v_new]=betarnd(u,v);
%    d_new=u/(u+v);u_new=u;v_new=v;
   X= gen_newstate_fn_clt(model,X_old,d_new-mu, mvnrnd(zeros(size(X_old,1)-3,1)',model.Q_clt,size(X_old,2))',u_new- X_old(7,:),v_new-X_old(8,:));
    %X= gen_newstate_fn_clt(model,X_old,betarnd((mu.*(1-mu)/sig -1).*mu,(mu.*(1-mu)/sig -1).*(1-mu))-mu, mvnrnd(zeros(size(X_old,1)-3,1)',model.Q_clt,size(X_old,2))',u_new- X_old(7,:),v_new-X_old(8,:));
end;