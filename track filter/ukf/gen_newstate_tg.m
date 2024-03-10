function [X,u_new,v_new]= gen_newstate_tg(filter,model,X_old,u,v)%Ìí¼Ó½Ó¿Ú
%--- this new state generation function follows the linear state spc eqn.
% x= Ax_old + Bv

if isempty(X_old),
    X= [];
    u_new=[];
    v_new=[];
else
    mu= X_old(1,:); mu=min(mu,0.999); mu=max(mu,0.001); X_old(1,:)= mu;
    sig= min(0.9*mu.*(1-mu),model.pdvarfac_tg);
    beta_avg_tmp= u./(u+v); 
    beta_avg_tmp=min(beta_avg_tmp,0.999); beta_avg_tmp=max(beta_avg_tmp,0.001);
    beta_var_tmp= (u.*v)./((u+v).^2.*(u+v+1)); 
    beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./min(filter.beta_factor.*beta_var_tmp,model.pdvarfac_tg)- 1;
%     beta_var_tmp=min(0.9*beta_avg_tmp.*(1-beta_avg_tmp),model.pdvarfac_tg);
%     beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./beta_var_tmp- 1;
    u= beta_tht_tmp.*beta_avg_tmp;  v= beta_tht_tmp.*(1-beta_avg_tmp);
    [r,u_new,v_new]=betarnd(u,v);
 %   r=u./(u+v);u_new=u;v_new=v;
    X= gen_newstate_fn_tg1(model,X_old,r, model.B*randn(size(model.B,2),size(X_old,2)),u_new,v_new);
   % X= gen_newstate_fn_tg(model,X_old,betarnd((mu.*(1-mu)./sig -1).*mu,(mu.*(1-mu)./sig -1).*(1-mu))-mu, model.B*randn(size(model.B,2),size(X_old,2)),u_new- X_old(7,:),v_new-X_old(8,:));
end;