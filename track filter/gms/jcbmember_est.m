%% This is the MATLAB code for the Robust CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, R. Hoseinnezhad, and R. Mahler "Robust Multi-Bernoulli Filtering," IEEE Journal on Selected Topics in Signal Processing, Vol. 7, No. 3, pp. 399-409, 2013.
% http://ba-ngu.vo-au.com/vo/VVHM_JSSP13.pdf
% ---BibTeX entry
% @ARTICLE{RobustCBMEMBER,
% author={B.-T. Vo and B.-N. Vo and R. Hoseinnezhad and R. Mahler},
% journal={IEEE Transactions on Signal Processing},
% title={Robust Multi-Bernoulli Filtering},
% year={2013},
% month={Jun},
% volume={7},
% number={3},
% pages={399-409}} 
%---

function cbmember_out = jcbmember_est(tt_birth , cbmember_in , model , filter ,meas , k)
%输出参数
cbmember_out=cbmember_in; 

%% 预测%%
  T_birth=length(tt_birth);
  % T_birth=4;
    T_predict= cbmember_in.T_update+T_birth+model.T_birth_clt;                                                                              %total number of tracks/components
    J_predict= zeros(T_predict,1);                                                                                  %total number of particles in each track/component
 %   r_predict= zeros(T_predict,1);                                                                                  %existence probability for tracks/components
    x_predict= cell(T_predict,1);                                                                                   %states of particles in each track/component
    w_predict= cell(T_predict,1);                                                                                   %weights of particles in each track/component
    l_predict= cell(T_predict,1);
    r_predict= zeros(T_birth,1); 
    u_b = cell(T_predict,1);     %s参数，初始化可能有问题
    v_b = cell(T_predict,1);     %t参数
    
    %新生目标预测
    for t=1:T_birth
        J_predict(t)= max(round( tt_birth{t}.r_b*filter.J_max),filter.J_min);                                                                            %append birth particles count
        %x_birth_temp1= gen_gms(tt_birth{t}.w,tt_birth{t}.m,tt_birth{t}.P,J_predict(t));
        x_birth_temp1= repmat( tt_birth{t}.m, [1, J_predict(t)])+tt_birth{t}.B*randn(model.x_dim,J_predict(t));
        [d_birth,u_b{t},v_b{t}]=betarnd(tt_birth{t}.u_b_tg,tt_birth{t}.v_b_tg,1,J_predict(t));
        x_predict{t}= [d_birth; x_birth_temp1;u_b{t};v_b{t}];                                                   %append birth particles
        w_predict{t}= ones(J_predict(t),1)/J_predict(t);                                                                                   %append birth weights
        l_predict{t}= ones(J_predict(t),1);
        r_predict(t)= tt_birth{t}.r_b;
    end     
    
    offset= T_birth;

    %杂波导致的新生目标预测
    r_predict= [r_predict; sum(model.r_birth_clt)/model.T_birth_clt*ones(model.T_birth_clt,1)];
%     r_predict= [r_predict; model.r_birth_clt];
    for t=1:model.T_birth_clt             
        J_predict(offset+t)= max(round(model.r_birth_clt(t)*filter.J_max/2),round(filter.J_min/2));                                                                            %append birth particles count
        x_birth_temp0= repmat(model.m_birth_clt{t}, [1, J_predict(offset+t)])+model.B_birth_clt{t}*randn(model.x_dim,J_predict(offset+t));%此处有修改
        [d_clt,u_b{offset+t},v_b{offset+t}]=betarnd(model.u_b_clt(t),model.v_b_clt(t),1,J_predict(offset+t));
        x_predict{offset+t}= [d_clt; x_birth_temp0;u_b{offset+t};v_b{offset+t}];                                                   %append birth particles
        w_predict{offset+t}= ones(J_predict(offset+t),1)/J_predict(offset+t);
        l_predict{offset+t}= zeros(J_predict(offset+t),1);
    end
    
    offset= T_birth+model.T_birth_clt; 
    
        %--- prediction of beta parameters预测贝塔参数
    for j=1:cbmember_in.T_update
%         beta_avg_tmp= cbmember_in.u_update{j}./(cbmember_in.u_update{j}+cbmember_in.v_update{j}); 
%         beta_var_tmp= (cbmember_in.u_update{j}.*cbmember_in.v_update{j})./((cbmember_in.u_update{j}+cbmember_in.v_update{j}).^2.*(cbmember_in.u_update{j}+cbmember_in.v_update{j}+1)); 
%         beta_tht_tmp= (beta_avg_tmp.*(1-beta_avg_tmp))./min(filter.beta_factor.*beta_var_tmp,beta_avg_tmp.*(1-beta_avg_tmp))- 1;
%         u_b{offset+j}= beta_tht_tmp.*beta_avg_tmp;  v_b{offset+j}= beta_tht_tmp.*(1-beta_avg_tmp);
          u_b{offset+j}= cbmember_in.u_update{j};  v_b{offset+j}= cbmember_in.v_update{j};
    end   
    u_predict =u_b;
    v_predict =v_b;
   
    %prediction for surviving tracks存活航迹预测（粒子实现）
    for t=1:cbmember_in.T_update
        idx1= find(cbmember_in.l_update{t} == 1);
        if ~isempty(idx1)
            w_temp1= cbmember_in.w_update{t}(idx1).*compute_pS_tg(model,cbmember_in.x_update{t}(:,idx1));%计算目标存活概率分量Ps
        else
            w_temp1=[];
        end
        idx0= find(cbmember_in.l_update{t} == 0);
        if ~isempty(idx0)
            w_temp0=cbmember_in. w_update{t}(idx0).*compute_pS_clt(model,cbmember_in.x_update{t}(:,idx0));%计算杂波存活概率分量Pu
        else
            w_temp0= []; 
        end
        
        r_predict(offset+t)=cbmember_in. r_update(t)* (sum(w_temp1)+sum(w_temp0));%r(k|k-1)=r(k-1)*<Pu,Ps>==r(k-1)*w*Pu,Ps
        J_predict(offset+t)= cbmember_in.J_update(t);
        w_predict{offset+t}= [w_temp1/(sum(w_temp1)+sum(w_temp0)); w_temp0/(sum(w_temp1)+sum(w_temp0))];%wp和wb的近似值；存储的预测的分量
        [x_new1,u_predict_idx1,v_predict_idx1]=gen_newstate_tg(filter,model,cbmember_in.x_update{t}(:,idx1), u_predict{offset+t}(:,idx1), v_predict{offset+t}(:,idx1));
        [x_new2,u_predict_idx0,v_predict_idx0]=gen_newstate_clt(filter,model,cbmember_in.x_update{t}(:,idx0),u_predict{offset+t}(:,idx0), v_predict{offset+t}(:,idx0));
        x_predict{offset+t}= [x_new1,x_new2];
        u_predict{offset+t}=[u_predict_idx1,u_predict_idx0];
        v_predict{offset+t}=[v_predict_idx1,v_predict_idx0];
        l_predict{offset+t}= [ones(length(idx1),1); zeros(length(idx0),1) ];
%         l_predict{offset+t}(idx1)=1; l_predict{offset+t}(idx0)=0;
    end

    r_predict= limit_range(r_predict);                                                                         %limit range of 0<r<1 for numerical stability
    
    %---construction of pseudo PHD for update
    J_pseudo= sum(J_predict);                                           %number of particles in pseudo-PHD
    x_pseudo= zeros(model.x_dim+3,J_pseudo);                              %states of particles in pseudo-PHD
    beta_pseudo=zeros(2,J_pseudo);
    l_pseudo= zeros(J_pseudo,1);
    w_pseudo= zeros(J_pseudo,1);                                        %weights of particles in pseudo-PHD
    w_pseudo1= zeros(J_pseudo,1);                                       %alt weight (1) of particles in pseudo-PHD - used in CB-MeMBer update later
    w_pseudo2= zeros(J_pseudo,1);                                       %alt weight (2) of particles in pseudo-PHD - used in CB-MeMBer update later
    
    start_pt= 1;
    for t=1:T_predict                  
        end_pt= start_pt+J_predict(t)-1;
        idx1= find(l_predict{t}==1);
        if ~isempty(idx1)
            w_temp1= w_predict{t}(idx1).*x_predict{t}(1,idx1)';
            x_temp1= x_predict{t}(:,idx1);
            beta_temp1=[u_predict{t}(:,idx1); v_predict{t}(:,idx1)];
        else
            w_temp1= [];
            x_temp1= [];
            beta_temp1=[];
        end
        idx0= find(l_predict{t}==0);
        if ~isempty(idx0)
            w_temp0= w_predict{t}(idx0).*x_predict{t}(1,idx0)';
            x_temp0= x_predict{t}(:,idx0);
            beta_temp0=[u_predict{t}(:,idx0); v_predict{t}(:,idx0)];
        else
            w_temp0= [];
            x_temp0= [];
            beta_temp0=[];
        end
        
        x_pseudo(:,start_pt:end_pt)= cat(2,x_temp1, x_temp0);
        beta_pseudo(:,start_pt:end_pt)=cat(2,beta_temp1, beta_temp0);
        l_pseudo(start_pt:end_pt)= cat(1,ones(length(idx1),1), zeros(length(idx0),1));
        w_pseudo(start_pt:end_pt) = cat(1,r_predict(t)/(1-r_predict(t))*w_predict{t}(idx1), r_predict(t)/(1-r_predict(t))*w_predict{t}(idx0));%计算分量
        w_pseudo1(start_pt:end_pt)= cat(1,r_predict(t)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))*w_predict{t}(idx1), r_predict(t)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))*w_predict{t}(idx0));
        w_pseudo2(start_pt:end_pt)= cat(1,r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))^2)*w_predict{t}(idx1), r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))^2)*w_predict{t}(idx0));

        start_pt= end_pt+1;
    end
    %reorder to adhere to tag convention - tag 1 particles, tag 0 particles
    idx1= find(l_pseudo==1); 
    idx0= find(l_pseudo==0); 
    x_pseudo= [x_pseudo(:,idx1) x_pseudo(:,idx0)];
    beta_pseudo= [beta_pseudo(:,idx1) beta_pseudo(:,idx0)];
    l_pseudo= [l_pseudo(idx1);l_pseudo(idx0)];
    w_pseudo= [w_pseudo(idx1);w_pseudo(idx0)];
    w_pseudo1= [w_pseudo1(idx1);w_pseudo1(idx0)];
    w_pseudo2= [w_pseudo2(idx1);w_pseudo2(idx0)];
    
    
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    
    T_update= T_predict+m;                                                                                        %total number of tracks/components
    J_update= zeros(T_update,1);                                                                                  %total number of particles in each track/component
    r_update= zeros(T_update,1);                                                                                  %existence probability for tracks/components
    x_update= cell(T_update,1);                                                                                   %states of particles in each track/component
    w_update= cell(T_update,1);                                                                                   %weights of particles in each track/component
    l_update= cell(T_update,1);
    u_update= cell(T_update,1);
    v_update= cell(T_update,1);
%       %pre calculation for Upsilon0 and Upsilon1
%     
%     missed_factor= 1- (sum(r_predict.*(u_predict./(u_predict+v_predict))));
%     if isnan(missed_factor), missed_factor=0; end %catch the degernerate zero case
%     terms0 = zeros(filter.N_max+1,1);
%     for n=0:filter.N_max
%         idxn= n+1;
%         if n < m
%             terms0(idxn) = eps(0);
%         else
%             terms0(idxn) =  exp(sum(log(1:n))-sum(log(1:n-m))+(n-m)*log(missed_factor));
%         end   
%     end
%     Upsilon0 = terms0;
%     
%     terms1 = zeros(filter.N_max,1);
%     for n=0:filter.N_max
%         idxn= n+1;
%         if n < m+1
%             terms1(idxn) = eps(0);
%         else
%             terms1(idxn) =  exp(sum(log(1:n))-sum(log(1:n-(m+1)))+(n-(m+1))*log(missed_factor));
%         end   
%     end
%     Upsilon1 = terms1;
%     
%     % misdetection term
%     m_update = m_predict;
%     P_update = P_predict;   
%     w_update= 1/(sum(w_predict)+sum(wc_predict))*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict)*(beta(u_predict,v_predict+1)./beta(u_predict,v_predict)).*w_predict;
%     u_update= u_predict; v_update= v_predict+1;
    
    %legacy tracks漏检
    for t=1:T_predict  
        idx1= find(l_predict{t}==1);
        if ~isempty(idx1)
            w_temp1= w_predict{t}(idx1).*x_predict{t}(1,idx1)';%检测率替换
            m_temp1= w_predict{t}(idx1).*(1-x_predict{t}(1,idx1)');
        else
            w_temp1= [];
            m_temp1= [];
        end   
        
        idx0= find(l_predict{t}==0);
        if ~isempty(idx0)
            w_temp0= w_predict{t}(idx0).*x_predict{t}(1,idx0)';
            m_temp0= w_predict{t}(idx0).*(1-x_predict{t}(1,idx0)');
        else
            w_temp0= [];
            m_temp0= [];
        end
        
        J_update(t)= J_predict(t);
        rtrack1= r_predict(t)*sum(m_temp1)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)));
        rtrack0= r_predict(t)*sum(m_temp0)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)));
        r_update(t)= rtrack1+rtrack0;%漏检概率更新
        r_comp1(t)= rtrack1;
        r_comp0(t)= rtrack0;
        
       w_update{t}= [m_temp1/(sum(m_temp1)+sum(m_temp0)) ; m_temp0/(sum(m_temp1)+sum(m_temp0))]; 
       %w_update{t}= [ w_predict{t}(idx1).*beta(uc_predict(t),vc_predict(t))./beta(uc_predict(t),vc_predict(t))/(sum(m_temp1)+sum(m_temp0)) ; w_predict{t}(idx0).*beta(uc_predict(t),vc_predict(t)+1)./beta(uc_predict(t),vc_predict(t))/(sum(m_temp1)+sum(m_temp0))];   
        x_update{t}= x_predict{t};
        l_update{t}= l_predict{t};
        u_update{t}= u_predict{t}; 
        v_update{t}= v_predict{t};
    end
    %measurement updated tracks
    if m~=0
        offset= T_predict;
        idx1= find(l_pseudo==1);
        idx0= find(l_pseudo==0);

        pD_vals_tg= x_pseudo(1,idx1)';%
        pD_vals_clt= x_pseudo(1,idx0)';

        for ell=1:m
            meas_likelihood_tg= compute_likelihood_tg(model,meas.Z{k}(:,ell),x_pseudo(:,idx1))';
            meas_likelihood_clt= compute_likelihood_clt(model,meas.Z{k}(:,ell),x_pseudo(:,idx0))';
          
            r_temp1= sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo2(idx1))/(sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo1(idx1))+sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo1(idx0)));
            r_temp0= sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo2(idx0))/(sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo1(idx1))+sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo1(idx0)));
            r_track= r_temp1+r_temp0;
            
            w_temp1= pD_vals_tg.*meas_likelihood_tg.*w_pseudo(idx1);
            w_temp0= pD_vals_clt.*meas_likelihood_clt.*w_pseudo(idx0);
                     
            r_update(offset+ell)= r_track;
            w_update{offset+ell}= [w_temp1/(sum(w_temp1)+sum(w_temp0)); w_temp0/(sum(w_temp1)+sum(w_temp0))];
           %  w_update{offset+ell}= [ meas_likelihood_tg.*w_pseudo(idx1).*beta(u_predict(offset+ell)+1,v_predict(offset+ell))./beta(u_predict(offset+ell),v_predict(offset+ell))/(sum(w_temp1)+sum(w_temp0)) ; meas_likelihood_clt.*w_pseudo(idx0).*beta(u_predict(offset+ell)+1,v_predict(offset+ell))./beta(u_predict(offset+ell),v_predict(offset+ell))/(sum(w_temp1)+sum(w_temp0))];   
            
            J_update(offset+ell)= J_pseudo;
            x_update{offset+ell}= x_pseudo;
            l_update{offset+ell}= l_pseudo;
            
            r_comp1(offset+ell)= r_temp1;
            r_comp0(offset+ell)= r_temp0;
            
            u_update{offset+ell}= beta_pseudo(1,:); 
            v_update{offset+ell}= beta_pseudo(2,:);
        end
        
    end
    r_update= limit_range(r_update);
    
    %管理航迹
     T_posterior= T_update;
     
     %重采样
     for t=1:T_update
        J_rsp(t)= max(round(r_update(t)*filter.J_max),filter.J_min);
        idx_rsp= randsample(length(w_update{t}),J_rsp(t),true,w_update{t}); %idx_rsp= resample(w_update{t},J_rsp(t));
        w_update{t}= ones(J_rsp(t),1)/J_rsp(t);
        x_update{t}= x_update{t}(:,idx_rsp);
        l_update{t}= l_update{t}(idx_rsp);
        u_update{t}= u_update{t}(:,idx_rsp);
        v_update{t}= v_update{t}(:,idx_rsp);
    end
    J_update= J_rsp;
    
    %航迹的修剪与连接
    [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,u_update, v_update]= prune_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,u_update, v_update,filter.track_threshold);  T_prune= length(r_update);
    [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,u_update, v_update]= cap_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,u_update, v_update,filter.T_max);              T_cap  = length(r_update);  
    
    T_update= T_cap;
    
 %   状态提取
    pD_tmp= [];
%     pD_tmp = (u_update(:)'./(u_update(:)'+v_update(:)'));
%     pD = sum(r_update*w_update(:)'*pD_tmp')/sum(r_update);
    cbmember_out.soft_N(k)= sum(r_comp1);
    cbmember_out.N(k)= round(cbmember_out.soft_N(k));
    [~,idx]= sort(-r_comp1);
    for t=1:min(cbmember_out.N(k),T_update)
        idx1= find(l_update{idx(t)}==1);
        cbmember_out.X{k} = [cbmember_out.X{k} x_update{idx(t)}(:,idx1)*w_update{idx(t)}(idx1)];
        pD_tmp = [pD_tmp; x_update{idx(t)}(1,idx1)*w_update{idx(t)}(idx1)];
    end
    cbmember_out.pD{k} = mean(pD_tmp);   
    
    %估计杂波率
        avg_pd0= zeros(T_update,1);
    for j=1:T_update
        idx0= find(l_update{j}==0);
        if ~isempty(idx0)
            avg_pd0(j)= (w_update{j}(idx0))'*(x_update{j}(1,idx0))';
        end
    end;
    cbmember_out.L(k)= sum(r_comp0'.*avg_pd0);   
   %  cbmember_out.L(k)=  sum(r_comp1);
   %输出参数
   cbmember_out.T_update= T_update;
   cbmember_out.r_update= r_update;
   cbmember_out.J_update=J_update;
   cbmember_out.w_update=w_update;
   cbmember_out.x_update=x_update;
   cbmember_out.l_update=l_update;
   cbmember_out.u_update=u_update;
   cbmember_out.v_update=v_update;
%    cbmember_out.uc_update=uc_update;
%    cbmember_out.vc_update=vc_update;
   
        %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #avg pD=' num2str( cbmember_out.pD{k}),...         
         ' #avg lambda=' num2str( cbmember_out.L(k)),...
         ' #eap target=' num2str(sum(r_comp1)),...
         ' #est card=' num2str( cbmember_out.N(k),4),...
         ' #trax updt=' num2str(T_posterior,4),...
         ' #trax elim=' num2str(T_prune,4),...
         ' #trax filt=' num2str(T_cap,4),...
         ' #samp filt=',num2str(sum(J_rsp))   ]);
    end
end

function clipped_r= limit_range(r)
r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;
end


function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l,new_u,new_v]= prune_tracks(r,r_com1,r_com0,w,x,J,l,u,v,track_threshold)

idx= find(r>track_threshold);

new_r= zeros(length(idx));
new_r_com1= zeros(length(idx));
new_r_com0= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_u= cell(length(idx),1);
new_v= cell(length(idx),1);
new_J= zeros(length(idx),1);
new_l= cell(length(idx),1);

new_r= r(idx); new_r_com1= r_com1(idx); new_r_com0= r_com0(idx);
for i=1:length(idx)
    new_w{i}= w{idx(i)};
    new_x{i}= x{idx(i)};
    new_l{i}= l{idx(i)};
    new_u{i}= u{idx(i)};
    new_v{i}= v{idx(i)};
end
new_J= J(idx);
end

function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l,new_u,new_v]= cap_tracks(r,r_com1,r_com0,w,x,J,l,u,v,T_max)

if length(r) > T_max

[~,idx]= sort(-r);

new_r= zeros(length(idx));
new_r_com1= zeros(length(idx));
new_r_com0= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_u= cell(length(idx),1);
new_v= cell(length(idx),1);
new_J= zeros(length(idx),1);
new_l= cell(length(idx),1);

new_r= r(idx(1:T_max)); new_r_com1= r_com1(idx(1:T_max)); new_r_com0= r_com0(idx(1:T_max));
for i=1:T_max
    new_w{i}= w{idx(i)};
    new_x{i}= x{idx(i)};
    new_l{i}= l{idx(i)};
    new_u{i}= u{idx(i)};
    new_v{i}= v{idx(i)};
end
new_J= J(idx(1:T_max));

else
    new_r= r;
    new_r_com1= r_com1;
    new_r_com0= r_com0;
    new_w= w;
    new_x= x;
    new_J= J;
    new_l= l;
    new_u= u;
    new_v= v;
end
end