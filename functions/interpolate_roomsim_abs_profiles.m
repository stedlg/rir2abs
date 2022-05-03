function [room_abs_profiles] = interpolate_roomsim_abs_profiles(roomsim_abs_coeffs,fqs,f_min_idx,f_max_idx)
%ALPHA_TO_F_ALPHA Summary of this function goes here
%   interpolation linéaire des log(beta)

B = [125 250 500 1000 2000 4000];
log_reflectivity_coeffs = log(1-roomsim_abs_coeffs);
f = fqs(f_min_idx:f_max_idx);
room_abs_profiles = zeros(size(f,2),size(roomsim_abs_coeffs,1));

for w = 1 : 1 : size(roomsim_abs_coeffs,1)
    
    for b = 1 : 1 : size(B,2)
       
        if b == 1
            
            f_b = f(f <= B(b));
            room_abs_profiles(find(ismember(f,f_b)),w) = roomsim_abs_coeffs(w,b);
            
        elseif b == size(B,2)
            
            f_b = f(f >= B(b));
            room_abs_profiles(find(ismember(f,f_b)) ,w) = roomsim_abs_coeffs(w,b);
            
            f_left = B(b-1);
            f_right = B(b);
            f_b = f((f >= f_left) & (f <= f_right));
            p = polyfit([f_left f_right],[log_reflectivity_coeffs(w,b-1) log_reflectivity_coeffs(w,b)],1); 
            room_abs_profiles(find(ismember(f,f_b)),w) = 1-exp(polyval(p,f_b));

        else
            
            f_left = B(b-1);
            f_right = B(b);
            
            f_b = f((f >= f_left) & (f <= f_right));
            p = polyfit([f_left f_right],[log_reflectivity_coeffs(w,b-1) log_reflectivity_coeffs(w,b)],1);
            
            room_abs_profiles(find(ismember(f,f_b)),w) = 1-exp(polyval(p,f_b));
            
        end

    end
end

end

