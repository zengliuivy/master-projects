function tt_birth = gen_meas_driven_birth_unline(model , meas_ru , Z , k)
meas_rb = (1 - meas_ru) / sum(1-meas_ru) ;
meas_rb = min(0.001,meas_rb) ; 
meas_rb = max(1e-10,meas_rb); 
NB = size(Z,2) ; 
tt_birth = cell(NB , 1) ; 
for bidx = 1 : NB
    tt_birth{bidx}.m = meas2state_unline(Z(:,bidx));
    tt_birth{bidx}.P = model.P_birth ; 
    tt_birth{bidx}.B =model.B_birth;
    tt_birth{bidx}.w = 1 ; 
    tt_birth{bidx}.r_b = meas_rb(bidx) ; 
    tt_birth{bidx}.l = [k+1;bidx] ; 
    tt_birth{bidx}.ah = [] ; 
    tt_birth{bidx}.st = model.init_st ; 
    tt_birth{bidx}.u_b_tg=90;
    tt_birth{bidx}.v_b_tg=10;
end
    