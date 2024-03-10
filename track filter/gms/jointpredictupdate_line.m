function [meas_ru_w , glmb_nextupdate] = jointpredictupdate_line(glmb_update,tt_birth,model,filter,meas,k,N_c,P_D)
%---generate next update

%create surviving tracks - via time prediction (single target CK)
tt_survive= cell(length(glmb_update.tt),1);                                                                                 %initialize cell array
for tabsidx=1:length(glmb_update.tt)
    [mtemp_predict,Ptemp_predict]= kalman_predict_multiple(model,glmb_update.tt{tabsidx}.m,glmb_update.tt{tabsidx}.P);      %kalman prediction for GM
    tt_survive{tabsidx}.m= mtemp_predict;                                                                                   %means of Gaussians for surviving track
    tt_survive{tabsidx}.P= Ptemp_predict;                                                                                   %covs of Gaussians for surviving track
    tt_survive{tabsidx}.w= glmb_update.tt{tabsidx}.w;                                                                       %weights of Gaussians for surviving track
    tt_survive{tabsidx}.l= glmb_update.tt{tabsidx}.l;                                                                       %track label
    tt_survive{tabsidx}.ah= glmb_update.tt{tabsidx}.ah;                                                                     %track association history (no change at prediction)
end

%create predicted tracks - concatenation of birth and survival
glmb_predict.tt= cat(1,tt_birth,tt_survive);                                                                                %copy track table back to GLMB struct

%gating by tracks
if filter.gate_flag
    for tabidx=1:length(glmb_predict.tt)
        glmb_predict.tt{tabidx}.gatemeas= gate_meas_gms_idx(meas.Z{k},filter.gamma,model,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P);
    end
else
    for tabidx=1:length(glmb_predict.tt)
        glmb_predict.tt{tabidx}.gatemeas= 1:size(meas.Z{k},2);
    end
end

%precalculation loop for average survival/death probabilities
NB = length(tt_birth) ; 
r_birth = zeros(NB , 1) ; 
for bidx = 1 : NB
    r_birth(bidx) = tt_birth{bidx}.r_b ; 
end
avps= [r_birth; zeros(length(glmb_update.tt),1)];
for tabidx=1:length(glmb_update.tt)
    avps(NB+tabidx)= model.P_S;
end
avqs= 1-avps;

%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);                                   %number of measurements
avpd= zeros(length(glmb_predict.tt),1);
tt_update= cell((1+m)*length(glmb_predict.tt),1);       %initialize cell array
%missed detection tracks (legacy tracks)
for tabidx= 1:length(glmb_predict.tt)
    tt_update{tabidx}= glmb_predict.tt{tabidx};         %same track table
    tt_update{tabidx}.ah= [tt_update{tabidx}.ah; 0];    %track association history (updated for missed detection)
    avpd(tabidx) = P_D ; 
end
avqd= 1-avpd;
%measurement updated tracks (all pairs)
allcostm= zeros(length(glmb_predict.tt),m);
for tabidx= 1:length(glmb_predict.tt)
    for emm= glmb_predict.tt{tabidx}.gatemeas
        stoidx= length(glmb_predict.tt)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted_tracks*j + i)
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k}(:,emm),model,glmb_predict.tt{tabidx}.m,glmb_predict.tt{tabidx}.P);   %kalman update for this track and this measurement
        w_temp= qz_temp.*glmb_predict.tt{tabidx}.w+eps;                                                                                 %unnormalized updated weights
        tt_update{stoidx}.m= m_temp;                                                                                                    %means of Gaussians for updated track
        tt_update{stoidx}.P= P_temp;                                                                                                    %covs of Gaussians for updated track
        tt_update{stoidx}.w= w_temp/sum(w_temp);                                                                                        %weights of Gaussians for updated track
        tt_update{stoidx}.l = glmb_predict.tt{tabidx}.l;                                                                                %track label
        tt_update{stoidx}.ah= [glmb_predict.tt{tabidx}.ah; emm];                                                                        %track association history (updated with new measurement)
        allcostm(tabidx,emm)= sum(w_temp);                                                                                              %predictive likelihood
    end
end
glmb_nextupdate.tt= tt_update;  

%copy track table back to GLMB struct
%joint cost matrix
jointcostm= [diag(avqs) ...
             diag(avps.*avqd) ...
             repmat(avps.*avpd,[1 m]).*allcostm/(N_c*model.pdf_c)];
%gated measurement index matrix
gatemeasidxs= zeros(length(glmb_predict.tt),m);
for tabidx= 1:length(glmb_predict.tt)
    gatemeasidxs(tabidx,1:length(glmb_predict.tt{tabidx}.gatemeas))= glmb_predict.tt{tabidx}.gatemeas;
end
gatemeasindc= gatemeasidxs>0;
         

%component updates
runidx= 1;
meas_ru = [] ;
for pidx=1:length(glmb_update.w)
    %calculate best updated hypotheses/components
    cpreds= length(glmb_predict.tt);
    nbirths= NB;
    nexists= length(glmb_update.I{pidx});
    ntracks= nbirths + nexists;
    tindices= [1:nbirths nbirths+glmb_update.I{pidx}'];                                                                                 %indices of all births and existing tracks  for current component
    lselmask= false(length(glmb_predict.tt),m); lselmask(tindices,:)= gatemeasindc(tindices,:);                                         %logical selection mask to index gating matrices
    mindices= unique_faster(gatemeasidxs(lselmask));                                                                                    %union indices of gated measurements for corresponding tracks
    costm= jointcostm(tindices,[tindices cpreds+tindices 2*cpreds+mindices]);                                                           %cost matrix - [no_birth/is_death | born/survived+missed | born/survived+detected]
    neglogcostm= -log(costm);                                                                                                           %negative log cost
    [uasses,nlcost]= gibbswrap_jointpredupdt_custom(neglogcostm,round(filter.H_upd*sqrt(glmb_update.w(pidx))/sum(sqrt(glmb_update.w))));%murty's algo/gibbs sampling to calculate m-best assignment hypotheses/components
    uasses(uasses<=ntracks)= -inf;                                                                                                      %set not born/track deaths to -inf assignment
    uasses(uasses>ntracks & uasses<= 2*ntracks)= 0;                                                                                     %set survived+missed to 0 assignment
    uasses(uasses>2*ntracks)= uasses(uasses>2*ntracks)-2*ntracks;                                                                       %set survived+detected to assignment of measurement index from 1:|Z|    
    uasses(uasses>0)= mindices(uasses(uasses>0));                                                                                       %restore original indices of gated measurements
    
    %generate corrresponding jointly predicted/updated hypotheses/components
    for hidx=1:length(nlcost)
        meas_ru_temp = false (size(meas.Z{k},2) , 1) ;
        update_hypcmp_tmp= uasses(hidx,:)'; 
        meas_ru_temp(update_hypcmp_tmp(update_hypcmp_tmp>0)) = true ; 
        meas_ru = [meas_ru , meas_ru_temp] ; 
        update_hypcmp_idx= cpreds.*update_hypcmp_tmp+[(1:nbirths)'; nbirths+glmb_update.I{pidx}];
        glmb_nextupdate.w(runidx)= -N_c+m*log(N_c*model.pdf_c)+log(glmb_update.w(pidx))-nlcost(hidx);                                             %hypothesis/component weight
        glmb_nextupdate.I{runidx}= update_hypcmp_idx(update_hypcmp_idx>0);                                                                                              %hypothesis/component tracks (via indices to track table)
        glmb_nextupdate.n(runidx)= sum(update_hypcmp_idx>0);                                                                                                            %hypothesis/component cardinality
        glmb_nextupdate.clutter(runidx) = m-sum(update_hypcmp_tmp>0) ; 
        runidx= runidx+1;
    end
end

glmb_nextupdate.w= exp(glmb_nextupdate.w-logsumexp(glmb_nextupdate.w));                                                                                                                 %normalize weights
meas_ru_w = meas_ru * glmb_nextupdate.w' ; 
%extract cardinality distribution
for card=0:max(glmb_nextupdate.n)
    glmb_nextupdate.cdn(card+1)= sum(glmb_nextupdate.w(glmb_nextupdate.n==card));                                                                                                       %extract probability of n targets
end

%remove duplicate entries and clean track table
glmb_nextupdate= clean_update(clean_predict(glmb_nextupdate));
end


function glmb_temp= clean_predict(glmb_raw)
%hash label sets, find unique ones, merge all duplicates
for hidx= 1:length(glmb_raw.w)
    glmb_raw.hash{hidx}= sprintf('%i*',sort(glmb_raw.I{hidx}(:)'));
end

[cu,~,ic]= unique(glmb_raw.hash);

glmb_temp.tt= glmb_raw.tt;
glmb_temp.w= zeros(length(cu),1);
glmb_temp.I= cell(length(cu),1);
glmb_temp.n= zeros(length(cu),1);
glmb_temp.clutter = zeros(length(cu),1);
for hidx= 1:length(ic)
        glmb_temp.w(ic(hidx))= glmb_temp.w(ic(hidx))+glmb_raw.w(hidx);
        glmb_temp.I{ic(hidx)}= glmb_raw.I{hidx};
        glmb_temp.n(ic(hidx))= glmb_raw.n(hidx);
        glmb_temp.clutter(ic(hidx))= glmb_raw.clutter(hidx);
end
glmb_temp.cdn= glmb_raw.cdn;
end

function glmb_clean= clean_update(glmb_temp)
%flag used tracks
usedindicator= zeros(length(glmb_temp.tt),1);
for hidx= 1:length(glmb_temp.w)
    usedindicator(glmb_temp.I{hidx})= usedindicator(glmb_temp.I{hidx})+1;
end
trackcount= sum(usedindicator>0);

%remove unused tracks and reindex existing hypotheses/components
newindices= zeros(length(glmb_temp.tt),1); newindices(usedindicator>0)= 1:trackcount;
glmb_clean.tt= glmb_temp.tt(usedindicator>0);
glmb_clean.w= glmb_temp.w;
glmb_clean.clutter= glmb_temp.clutter ; 
for hidx= 1:length(glmb_temp.w)
    glmb_clean.I{hidx}= newindices(glmb_temp.I{hidx});
end
glmb_clean.n= glmb_temp.n;
glmb_clean.cdn= glmb_temp.cdn;
end

