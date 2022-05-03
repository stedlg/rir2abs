function [toas_m_s,idx_toas_m_s_0,i_stft,idx_toas_m_s,n_img_sources_in_window] = toas_ioas_temp_and_stft(toas,n_sources,m,s,window_size,overlap,fs,n_features,size_rir)
%TIME_AXIS Summary of this function goes here
%   Detailed explanation goes here

previous = (m-1)*n_features*n_sources+(s-1)*n_features;
toas_m_s = toas(1,previous+1:previous+n_features);


idx_toas_m_s = ceil(toas_m_s*fs)+1;
idx_toas_m_s_0 = idx_toas_m_s;

% stft windows timings correspond to the center
i_stft = round(window_size/2) + (0 : (size_rir-window_size)/(window_size-overlap)).*(window_size-overlap);

left_border = i_stft - round(window_size/2);
right_border = i_stft + round(window_size/2);

n_img_sources_in_window = zeros(size(i_stft));

for q = 1 : 1 :size(toas_m_s,2)
    
    %idx of arrival (ioa) after stft : idx of the win whose timing (center) is
    %the closer to the toa
    
    idx_toas_m_s(q) = min(find(abs(i_stft-idx_toas_m_s(q)) == min(abs(i_stft-idx_toas_m_s(q))))) ;
    n_img_sources_in_window = n_img_sources_in_window + ((idx_toas_m_s_0(q)>=left_border)&(idx_toas_m_s_0(q)<=right_border));
end

end


