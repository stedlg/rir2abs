function [Pi] = build_Pi(reflection_orders,n_features,toas,sigma_geo,window_size,overlap,fs,size_rir,sound_speed,index_split_X,n_receivers,n_sources,L_v)
%P Summary of this function goes here
%   Detailed explanation goes here

%duration the sound takes to travel 2cm

[~,~,i_stft] = toas_ioas_temp_and_stft(toas,n_sources,1,1,window_size,overlap,fs,n_features,size_rir);

geometrical_std = sigma_geo/sound_speed;
l_k = reflection_orders'+2;
J = index_split_X;
a_j = (i_stft - round(window_size/2))/fs;
b_j = (i_stft + round(window_size/2))/fs;

a_j = a_j(1:J) ; b_j = b_j(1:J);
a_j = repmat(a_j,n_features,1) ; b_j = repmat(b_j,n_features,1);

b_j_ = b_j - L_v; 
a_j_ = a_j - L_v;


Pi = [];

for m = 1 : n_receivers
    for s = 1: n_sources
        
        [activations,~] = toas_ioas_temp_and_stft(toas,n_sources,m,s,window_size,overlap,fs,n_features,size_rir);
        Erf_b_j = erf((b_j-activations')./(sqrt(2)*l_k'*geometrical_std));
        Erf_a_j = erf((a_j-activations')./(sqrt(2)*l_k'*geometrical_std));
        
        Erf_b_j_ = erf((b_j_-activations')./(sqrt(2)*l_k'*geometrical_std));
        Erf_a_j_ = erf((a_j_-activations')./(sqrt(2)*l_k'*geometrical_std));
        
        pi_k_j_ = (1/2)*(Erf_b_j_-Erf_a_j);
        pi_k_j_plus = (1/2)*(Erf_b_j-Erf_a_j_);
        
        Pi_temp = [];
        
        for k = 1 : 1 : n_features
            pi_k_j_comp = pi_k_j_plus;
            pi_k_j_comp(k,:) = [];
            Pi_temp = [ Pi_temp ; prod([pi_k_j_(k,:) ; 1-pi_k_j_comp ],1)];
        end
        
        Pi = [Pi Pi_temp];
    end
end

end

