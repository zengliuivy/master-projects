function Z= gen_observation(model,X,W)

%linear observation equation (position components only)

if ~isnumeric(W)
    if strcmp(W,'noise')
        W= model.D*randn(size(model.D,2),size(X,2));
    elseif strcmp(W,'noiseless')
        W= zeros(size(model.D,1),size(X,2));
    end
end

if isempty(X)
    Z= [];
else
   Z= model.H*X+ W;
%modify below here for user specified measurement model
%     P= X([2 4],:);
%     Z(1,:)= atan2(P(1,:),P(2,:));   
%     Z(2,:)= sqrt(sum(P.^2));
%     Z= Z+ W;
end