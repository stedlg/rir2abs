function [H] = build_H(n_sources,n_receivers,size_rir,n_features,sound_speed,window_size,overlap,r,t_activations,index_split_X,fs)
%BUILD_H_ANALYTIQUE Summary of this function goes here
%   Detailed explanation goes here

H = [];

for m = 1 : n_receivers
    for s = 1: n_sources    
        h = []; 
        [activations,~,i_stft,~] = toas_ioas_temp_and_stft(t_activations,n_sources,m,s,window_size,overlap,fs,n_features,size_rir);
        %% ordonner les doublons
        
        for k = 1 : 1 : n_features
            if r == 0 % rectangular window
               
                h_k = ones(size(i_stft));
                h_k = h_k * (1/(sound_speed*(activations(k))))^2;
                h_k = h_k(1:index_split_X);
                
            else
                disp("please use rectangle windows for stft")
                break
            end
            
            h = [h ; h_k];
            
        end
     
        H = [H h];      
    end   
end
end

