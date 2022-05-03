function [paramFinal,fval] = build_optimizer(x,H,J,init_param,K,Q)
%BUILD_W_OPTIM Summary of this function goes here
%   Detailed explanation goes here

reflection_orders = sum(K,2);
%% build
list_index = [1 : 1 : size(H,2)];
H_J = H;
X_J = repmat(x,size(K,1),1);
%% initialization
paramInit = zeros(1,7);
for i = 1 : 1 : size(K,1)
    
    H_J(i,~ismember(list_index,J(i).inliers)) = 0;
    X_J(i,~ismember(list_index,J(i).inliers)) = 0;
    
    if reflection_orders(i)<= 1
        
        if init_param
            
            if reflection_orders(i)==0
                paramInit(1) = sum((H_J(i,:).*X_J(i,:)))/sum(H_J(i,:).^2);
            else
                
                paramInit((find(K(i,:)==1))+1) = paramInit((find(K(i,:)==1))+1) + sum(((paramInit(1)*H_J(i,:)).*X_J(i,:)))/sum((paramInit(1)*H_J(i,:)).^2);
            end
        end
 
    end
end
%% Q
mask = find(reflection_orders<=Q);
X_J(~mask,:) = [];
H_J(~mask,:) = [];
K(~mask,:) = [];

%% constraints : Cx<b
Aoptimiser=@(param)loss(param,X_J,H_J,K);
C = [eye(size(paramInit,2)) ; [zeros(size(paramInit,2)-1,1) (-1)*eye(size(paramInit,2)-1)]] ;
b = [zeros(size(paramInit,2),1) ; -ones(size(paramInit,2)-1,1)];
options = optimoptions('fmincon','Display','none','Algorithm','sqp');
Aeq = []; beq = []; lb = []; ub = []; nonlcon = [];

paramInit(isnan(paramInit)) = 0;
[paramFinal,fval] = fmincon(Aoptimiser,paramInit',(-1)*C,-b,Aeq , beq , lb , ub , nonlcon,options);
