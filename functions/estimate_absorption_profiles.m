function [targets_abs_coeffs_i,estimates_abs_coeffs_i,not_estimated_abs_coeffs_i] = estimate_absorption_profiles(data,fs,sound_speed,Q_max,Q,n_sources,n_receivers,size_rir,Nm,noise,psnr_,sigma_geo,L_v,window_size,overlap,fft_length,r,fqs,f_min_idx,f_max_idx,N_iter,threshold_ransac,threshold_Pi,natural_score,graphic_evaluation,wall_labels)
%ESTIMATE_ABSORPTION_PROFILES Summary of this function goes here
%   Detailed explanation goes here

targets_abs_coeffs_i = zeros(6,size(fqs,2));
estimates_abs_coeffs_i = zeros(6,size(fqs,2));
not_estimated_abs_coeffs_i = zeros(6,size(fqs,2));

room_abs = reshape(data.absorption,6,6).';
    
    %% model geometrical noise
    
    ism.dim = data.dim + normrnd(0,sigma_geo,size(data.dim));
    ism.receivers = reshape([data.receiver.location],3,[])' + normrnd(0,sigma_geo,size(reshape([data.receiver.location],3,[])'));
    ism.sources = reshape([data.source.location],3,[])' + normrnd(0,sigma_geo,size(reshape([data.source.location],3,[])'));
    %% compute toas with image-source model
    %disp('calculation of toas with image-source model')
    
    [toas,K] = compute_ism_toas_and_K(ism.receivers,ism.sources,Nm,ism.dim,sound_speed,n_receivers,n_sources,Q_max);
    reflection_orders = sum(K,2);
    n_features = size(K,1);
    %% split spectrograms beyond highest toa
    
    idx_highest_toa= ceil(max(toas)*fs)+1;
    i_stft = round(window_size/2) + (0 : (size_rir-window_size)/(window_size-overlap)).*(window_size-overlap);
    idx_highest_toa_X = min(find(abs(i_stft-idx_highest_toa) == min(abs(i_stft-idx_highest_toa)))) ;
    idx_split_X = idx_highest_toa_X + round(window_size/2) + 1 ;
    clear i_stft idx_highest_toa idx_highest_toa_X   
    %% H : compensate geometrical spreading
    [H] = build_H(n_sources,n_receivers,size_rir,n_features,sound_speed,window_size,overlap,r,toas,idx_split_X,fs);
    %% Pi : RIR windows validity probabilities
    [Pi] = build_Pi(reflection_orders,n_features,toas,sigma_geo,window_size,overlap,fs,size_rir,sound_speed,idx_split_X,n_receivers,n_sources,L_v);
    %% X : RIR spectrograms
    [X] = build_X(data.RIR,fs,window_size,overlap,fft_length,idx_split_X,n_receivers,n_sources,r,f_min_idx,f_max_idx,noise,psnr_);
    %% ransac
    
    paramFinal = [];
    %disp('ransac algorithm and optimization')
    for f =  1 : 1 : size(X,1)
        
        [J_ransac_f,~] = build_J_ransac(f,N_iter,Pi,X,H,threshold_ransac,threshold_Pi,natural_score,K);
       
        [paramFinal_f,fval] = build_optimizer(X(f,:),H,J_ransac_f,true,K,Q);
        paramFinal = [paramFinal;paramFinal_f'];
        
        for w = 1 : 1 : size(room_abs,1)
            if size(J_ransac_f(w+1).inliers,2) == 0
                not_estimated_abs_coeffs_i(w,f) = not_estimated_abs_coeffs_i(w,f)+1;
            end
        end
    end
    %% evaluation
    
    f = fqs(f_min_idx:f_max_idx);
    room_abs_profiles = interpolate_roomsim_abs_profiles(room_abs,fqs,f_min_idx,f_max_idx);

    for w = 1 : 1 : size(room_abs,1)
        targets_abs_coeffs_i(w,:) = room_abs_profiles(:,w);
        estimates_abs_coeffs_i(w,:) = 1-(paramFinal(:,w+1));
    end
    
    if graphic_evaluation
        fig = figure('WindowState','maximized');
        
        for w = 1 : 1 : size(room_abs,1)

            subplot(2,3,w);
            
            plot((fqs(f_min_idx:f_max_idx)),reshape(targets_abs_coeffs_i(w,:),1,size(fqs,2)),'o-','linewidth',2)
            hold on
            
            plot((fqs(f_min_idx:f_max_idx)),reshape(estimates_abs_coeffs_i(w,:),1,size(fqs,2)),'o-','linewidth',2)
            ylim([0 1])
            
            xlabel('f','fontsize',40)
            ylabel("\alpha(f)",'fontsize',40)
            
            title(strcat(wall_labels(w)),'fontsize',35);
            if w == size(room_abs,1)
                legend(["target", "estimate"],'location','northeast','FontSize',40,'Location','northeast')
            end
            
        end
        pause()
        close all
        
    end
    %clc
end

