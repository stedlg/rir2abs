function [X,stds] = build_X(D,fs,window_size,overlap,fft_length,index_split_X,n_receivers,n_sources,r,f_split,f_split_max,noise,psnr_)
%BUILD_V_M Summary of this function goes here
%   Detailed explanation goes here

X = [];
%peak_rirs = 0;

% for m = 1:1:n_receivers
%     
%     for s = 1:1:n_sources
%         
%         rir = D{s,m};
%         
%         if peak_rirs <= max(abs(rir))
%             peak_rirs = max(abs(rir));
%         end
%         
%     end
%     
% end


for m = 1:1:n_receivers
    
    for s = 1:1:n_sources
        
        rir = D{s,m};
        peak_rirs = max(abs(rir));
        if noise
            stds = peak_rirs .* sqrt(10.^(-psnr_.'/10));
            rir = rir + normrnd(0,ones(size(rir)).*stds);
        end
        
        rir_stft = stft(rir,fs,'Window',tukeywin(window_size,r),'OverlapLength',overlap,'FFTLength',fft_length,'FrequencyRange',"onesided"); %hann(window_size,'periodic')
        
        X_I = abs(rir_stft(f_split:f_split_max,1:index_split_X)).^2;
        
        X = [X X_I];
        
    end
end
end

