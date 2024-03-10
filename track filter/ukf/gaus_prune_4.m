function [w_new,x_new,P_new,u_new,v_new]= gaus_prune_4(w,x,P,u,v,elim_threshold)

idx1= find( w > elim_threshold );
w_new= w(idx1);
x_new= x(:,idx1);
P_new= P(:,:,idx1);
u_new= u(idx1);
v_new= v(idx1);
