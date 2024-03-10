% This is the MATLAB code for the CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, and A. Cantoni, "The Cardinality Balanced Multi-target Multi-Bernoulli filter and its implementations," IEEE Trans. Signal Processing, Vol. 57, No. 2, pp. 409?23, 2009. 
% http://ba-ngu.vo-au.com/vo/VVCmemberSP09.pdf
% ---BibTeX entry
% @ARTICLE{CBMEMBER,
% author={B.-T. Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={The Cardinality Balanced Multi-Target Multi-Bernoulli Filter and Its Implementations},
% year={2009},
% month={Feb},
% volume={57},
% number={2},
% pages={409-423}} 
%---


function cbmember_out = cbmember_est(tt_birth , cbmember_in , model , filter ,meas , k)
%输出参数
cbmember_out.X= cbmember_in.X;
cbmember_out.N= cbmember_in.N;
cbmember_out.L= cbmember_in.L;
cbmember_out.pD= cbmember_in.pD;

%% 预测%%
  %T_birth=length(tt_birth);
   T_birth=4;
    T_predict= cbmember_in.T_update+T_birth;                                                                            %total number of tracks/components
    L_predict= zeros(T_predict,1);                                                                                  %total number of Gaussians in each track/component
    r_predict= zeros(T_predict,1);                                                                                  %existence probability for tracks/components
    m_predict= cell(T_predict,1);                                                                                   %means of Gaussians in each track/component
    P_predict= cell(T_predict,1);                                                                                   %covs of Gaussians in each track/component
    w_predict= cell(T_predict,1);                                                                                   %weights of Gaussians in each track/component
    
    %预测存活的目标航迹
   for t=1:cbmember_in.T_update
        L_predict(t)= cbmember_in.L_update(t);                                                                                  %surviving number of Gaussians in current track 
        r_predict(t)= model.P_S*cbmember_in.r_update(t);                                                                        %surviving existence probability
        [m_predict{t},P_predict{t}] = kalman_predict_multiple(model,cbmember_in.m_update{t},cbmember_in.P_update{t});                       %surviving components
        w_predict{t}= cbmember_in.w_update{t};                                                                                  %surviving weights
   end  
   
   %插入出生航迹
     offset= cbmember_in.T_update;
    for t=1:T_birth
        L_predict(offset+t)= 1;                                                                      %append birth Gaussian count强制为1 
        r_predict(offset+t)= tt_birth{t}.r_b;                                                                      %append birth probabilities
        m_predict{offset+t}= tt_birth{t}.m; P_predict{offset+t}= tt_birth{t}.P;                               %append birth components
        w_predict{offset+t}= tt_birth{t}.w;                                                                      %append birth weights 
    end                     
    
    r_predict= limit_range(r_predict);                                                                         %limit range of 0<r<1 for numerical stability

%构建用于更新的伪▁PHD
    L_pseudo= sum(L_predict);                                           %number of Gaussians in pseudo-PHD
    m_pseudo= zeros(model.x_dim,L_pseudo);                              %means of Gaussians in pseudo-PHD
    P_pseudo= zeros(model.x_dim,model.x_dim,L_pseudo);                  %covs of Gaussians in pseudo-PHD
    w_pseudo= zeros(L_pseudo,1);                                        %weights of Gaussians in pseudo-PHD
    w_pseudo1= zeros(L_pseudo,1);                                       %alt weight (1) of Gaussians in pseudo-PHD - used in CB-MeMBer update later
    w_pseudo2= zeros(L_pseudo,1);                                       %alt weight (2) of Gaussians in pseudo-PHD - used in CB-MeMBer update later
    
    start_pt= 1;
    for t=1:T_predict                  
        end_pt= start_pt+L_predict(t)-1;
        m_pseudo(:,start_pt:end_pt)= m_predict{t};
        P_pseudo(:,:,start_pt:end_pt)= P_predict{t};
        w_pseudo(start_pt:end_pt) = r_predict(t)/(1-r_predict(t))*w_predict{t};
        w_pseudo1(start_pt:end_pt)= r_predict(t)/(1-r_predict(t)*model.P_D)*w_predict{t};
        w_pseudo2(start_pt:end_pt)= r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*model.P_D)^2)*w_predict{t};
        start_pt= end_pt+1;
    end
    
    
    %gating
     if filter.gate_flag
        meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_pseudo,P_pseudo);        
     end
    
    %% 更新 %%
    %量测的数量
    m= size(meas.Z{k},2);
    
    T_update= T_predict+m;                                                                                        %total number of tracks/components
    L_update= zeros(T_update,1);                                                                                  %total number of Gaussians in each track/component
    r_update= zeros(T_update,1);                                                                                  %existence probability for tracks/components
    m_update= cell(T_update,1);                                                                                   %means of Gaussians in each track/component
    P_update= cell(T_update,1);                                                                                   %covs of Gaussians in each track/component
    w_update= cell(T_update,1);                                                                                   %weights of Gaussians in each track/component
    
    %遗漏的航迹
    L_update= L_predict;
    r_update= r_predict.*((1-model.P_D)./(1-r_predict*model.P_D)); 
    m_update= m_predict;
    P_update= P_predict;
    w_update= w_predict;
        
    
     %用量测来更新航迹
     if m~=0
        offset= T_predict;
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_pseudo,P_pseudo);
        for ell=1:m                                        
            r_temp= (model.P_D*w_pseudo2(:)'*qz_temp(:,ell))/(model.lambda_c*model.pdf_c+model.P_D*w_pseudo1(:)'*qz_temp(:,ell));
            w_temp= model.P_D*w_pseudo(:).*qz_temp(:,ell); w_temp= w_temp/sum(w_temp);
            L_temp= length(w_temp);
            
            L_update(offset+ell)= L_temp;
            r_update(offset+ell)= r_temp;
            m_update{offset+ell}= m_temp(:,:,ell);
            P_update{offset+ell}= P_temp;
            w_update{offset+ell}= w_temp;    
        end
        
    end
    
    T_update= T_predict+m;
    r_update= limit_range(r_update);   
    
    %管理航迹
     T_posterior= T_update;
%      
%      %重采样
%      for t=1:T_update
%         J_rsp(t)= max(round(r_update(t)*filter.J_max),filter.J_min);
%         idx_rsp= randsample(length(w_update{t}),J_rsp(t),true,w_update{t}); %idx_rsp= resample(w_update{t},J_rsp(t));
%         w_update{t}= ones(J_rsp(t),1)/J_rsp(t);
%         x_update{t}= x_update{t}(:,idx_rsp);
%         l_update{t}= l_update{t}(idx_rsp);
%     end
%     J_update= J_rsp;
    
    %航迹的修剪与连接
    [r_update,w_update,m_update,P_update,L_update]= prune_tracks(r_update,w_update,m_update,P_update,L_update,filter.track_threshold);  T_prune= length(r_update);
    [r_update,w_update,m_update,P_update,L_update]= cap_tracks(r_update,w_update,m_update,P_update,L_update,filter.T_max);              T_cap  = length(r_update);  
    
    T_update= T_cap;
    
    %(within tracks) pruning,merging,capping
    for t=1:T_update
        [w_update{t},m_update{t},P_update{t}]= gaus_prune(w_update{t},m_update{t},P_update{t},filter.elim_threshold);    L_prune(t)= length(w_update{t});
        [w_update{t},m_update{t},P_update{t}]= gaus_merge(w_update{t},m_update{t},P_update{t},filter.merge_threshold);   L_merge(t)= length(w_update{t});
        [w_update{t},m_update{t},P_update{t}]= gaus_cap(w_update{t},m_update{t},P_update{t},filter.L_max);               L_cap(t)  = length(w_update{t});
    end
    
    L_update= L_cap;
    
%     [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update]= prune_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,filter.track_threshold);  T_prune= length(r_update);
%     [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update]= cap_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,filter.T_max);              T_cap  = length(r_update);  
%     
%     T_update= T_cap;
    
    %状态提取
%     pD_tmp= [];
%     cbmember_out.soft_N(k)= sum(r_comp1);
%     cbmember_out.N(k)= round(cbmember_out.soft_N(k));
%     [~,idx]= sort(-r_comp1);
%     for t=1:min(cbmember_out.N(k),T_update)
%         idx1= find(l_update{idx(t)}==1);
%         cbmember_out.X{k} = [cbmember_out.X{k} x_update{idx(t)}(:,idx1)*w_update{idx(t)}(idx1)];
%         pD_tmp = [pD_tmp; x_update{idx(t)}(1,idx1)*w_update{idx(t)}(idx1)];
%     end
%     cbmember_out.pD{k} = mean(pD_tmp);   
%     
%     %估计杂波率
%         avg_pd0= zeros(T_update,1);
%     for j=1:T_update
%         idx0= find(l_update{j}==0);
%         if ~isempty(idx0)
%             avg_pd0(j)= (w_update{j}(idx0))'*(x_update{j}(1,idx0))';
%         end
%     end;
%     cbmember_out.L(k)= sum(r_comp0'.*avg_pd0);   
     cdn_update= prod(1-r_update)*esf(r_update./(1-r_update));
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    cbmember_out.N(k) = min(length(r_update),map_cdn);
    [~,idx_max_tgt]= sort(-r_update);
    for t=1:cbmember_out.N(k)
         [~,idx]= max(w_update{idx_max_tgt(t)});
        cbmember_out.X{k} = [cbmember_out.X{k} m_update{idx_max_tgt(t)}(:,idx)];
        cbmember_out.L{k}= [];
    end
    
        %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #eap cdn=' num2str(sum(r_update)),...
         ' #var cdn=' num2str(sum(r_update.*(1-r_update)),4),...
         ' #est card=' num2str( cbmember_out.N(k),4),...
         ' #trax updt=' num2str(T_posterior,4),...
         ' #trax elim=' num2str(T_prune,4),...
         ' #trax filt=' num2str(T_cap,4),...
         ' #gaus filt=',num2str(sum(L_cap))   ]);
    end
  
end

function clipped_r= limit_range(r)
r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;
end

function [new_r,new_w,new_m,new_P,new_L]= prune_tracks(r,w,m,P,L,track_threshold)

idx= find(r>track_threshold);

new_r= zeros(length(idx));
new_m= cell(length(idx),1);
new_P= cell(length(idx),1);
new_w= cell(length(idx),1);
new_L= zeros(length(idx),1);

new_r= r(idx);
for i=1:length(idx)
    new_m{i}= m{idx(i)};
    new_P{i}= P{idx(i)};
    new_w{i}= w{idx(i)};
end
new_L= L(idx);
end

function [new_r,new_w,new_m,new_P,new_L]= cap_tracks(r,w,m,P,L,T_max)

if length(r) > T_max
    
    [~,idx]= sort(-r);
    
    new_r= zeros(length(idx));
    new_m= cell(length(idx),1);
    new_P= cell(length(idx),1);
    new_w= cell(length(idx),1);
    new_L= zeros(length(idx),1);
    
    new_r= r(idx(1:T_max));
    for i=1:T_max
        new_m{i}= m{idx(i)};
        new_P{i}= P{idx(i)};
        new_w{i}= w{idx(i)};
    end
    new_L= L(idx(1:T_max));
    
else
    new_r= r;
    new_m= m;
    new_P= P;
    new_w= w;
    new_L= L;
end
end
% function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l]= prune_tracks(r,r_com1,r_com0,w,x,J,l,track_threshold)
% 
% idx= find(r>track_threshold);
% 
% new_r= zeros(length(idx));
% new_r_com1= zeros(length(idx));
% new_r_com0= zeros(length(idx));
% new_w= cell(length(idx),1);
% new_x= cell(length(idx),1);
% new_J= zeros(length(idx),1);
% new_l= cell(length(idx),1);
% 
% new_r= r(idx); new_r_com1= r_com1(idx); new_r_com0= r_com0(idx);
% for i=1:length(idx)
%     new_w{i}= w{idx(i)};
%     new_x{i}= x{idx(i)};
%     new_l{i}= l{idx(i)};
% end
% new_J= J(idx);
% end
% 
% function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l]= cap_tracks(r,r_com1,r_com0,w,x,J,l,T_max)
% 
% if length(r) > T_max
% 
% [~,idx]= sort(-r);
% 
% new_r= zeros(length(idx));
% new_r_com1= zeros(length(idx));
% new_r_com0= zeros(length(idx));
% new_w= cell(length(idx),1);
% new_x= cell(length(idx),1);
% new_J= zeros(length(idx),1);
% new_l= cell(length(idx),1);
% 
% new_r= r(idx(1:T_max)); new_r_com1= r_com1(idx(1:T_max)); new_r_com0= r_com0(idx(1:T_max));
% for i=1:T_max
%     new_w{i}= w{idx(i)};
%     new_x{i}= x{idx(i)};
%     new_l{i}= l{idx(i)};
% end
% new_J= J(idx(1:T_max));
% 
% else
%     new_r= r;
%     new_r_com1= r_com1;
%     new_r_com0= r_com0;
%     new_w= w;
%     new_x= x;
%     new_J= J;
%     new_l= l;
% end
% end