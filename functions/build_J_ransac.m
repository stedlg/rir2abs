function [J_ransac_f,best_S] = build_J_ransac(f,N_iter,Pi,X,H,threshold_ransac,threshold_Pi,natural_score,K)
% Probablistic RIR window selection procedure inspired by RANSAC
%% Algorithm

% Initialization

best_S = 0;

reflection_orders = sum(K,2);

% indices of image sources with orders lower than 1
indices_v_0_1 = find(reflection_orders<=1);

for n = 1:1:N_iter
    
    model = (-1)*ones(1,size(K,2)+1); % [d, w_1, ... w_6]
    S = 0; 
    
    % a tentative room acoustic model is sampled

    for i = 1 : 1 : size(indices_v_0_1,1)
        k = indices_v_0_1(i);
        % windows with probabilities lower than 0.1% are discarded
        
        if n == 1
            struct_Pi(k).cols = find(Pi(k,:)> threshold_Pi);
            struct_Pi(k).probas = Pi(k,struct_Pi(k).cols)/sum(Pi(k,struct_Pi(k).cols));
        end
        
        if size(struct_Pi(k).cols,2) ~= 0
            
            %  windows are sampled from probabilities proportional to Pi_k
            model_k = datasample(struct_Pi(k).cols,1,'Weights',struct_Pi(k).probas);
            
            if reflection_orders(k) == 0
                %model = [model X(f,model_k)./H(k,model_k)];
                model(1) = X(f,model_k)./H(k,model_k);
            else
                %model = [model X(f,model_k)./(H(k,model_k)*model(1))];
                model(find(K(k,:)==1)+1) = X(f,model_k)./(H(k,model_k)*model(1));
            end
            
        else
            continue % equal to -1 when it is impossible to sample a relevant window for k
        end
        
    end
    
    if n == 1
        
        indices_v_qt_1 = find(reflection_orders>1);

        for i = 1 : 1 : size(indices_v_qt_1,1)
            k = indices_v_qt_1(i);
            struct_Pi(k).cols = find(Pi(k,:)> threshold_Pi);
        end
    end
    
    % selection of windows matching the room acoustic model within a relative error of threshold_ransac (default : 10%).
    for k = 1 : 1 : size(H,1)
        
        if ismember(-1,model(find(K(k,:)~=0)+1)) 
            continue
        end

        v_k = model(1)*prod(model(2:end).^K(k,:));
        i_inliers_k = find(abs((X(f,struct_Pi(k).cols)-v_k*H(k,struct_Pi(k).cols))./(v_k*H(k,struct_Pi(k).cols)+eps)) < threshold_ransac);
        struct_Pi(k).inliers = struct_Pi(k).cols(i_inliers_k);
        
        if natural_score % count inliers
            S = S + size(i_inliers_k,2);
        else % score used in the article
            S = S + sum(unique(H(k,struct_Pi(k).inliers)));
        end
    end
    
    if S > best_S
        
        best_S = S ;
        J_ransac_f = rmfield(struct_Pi,{'cols','probas'});
        
    end
    
end
end

