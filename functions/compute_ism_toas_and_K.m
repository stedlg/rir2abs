function [all_toas,K] = compute_ism_toas_and_K(receivers_pos,sources_pos,Nm,dim,sound_speed,n_receivers,n_sources,Q_max)
%UNTITLED Summary of this function goes here
% ﻿Echoes labelling with Image Source Method
% Source : Rigid sphere room impulse response simulation: Algorithm and applications
% D. P. Jarrett, E. A. P. Habets, M. R. P. Thomas, and P. A. Naylor

[x,y,z]=meshgrid(0:1);
P=unique(sort([x(:),y(:),z(:)],3),'rows');

[x,y,z]=meshgrid(-Nm:1:Nm);
M=unique(sort([x(:),y(:),z(:)],3),'rows');

all_toas = [];

for m_ = 1 : 1 : n_receivers
    
    rm = receivers_pos(m_,:);
    for s = 1 : 1 : n_sources
        rs = sources_pos(s,:);
        
        toas = [];
        labels = [];
        
        for i_p = 1 : 1 : size(P,1)
            p = P(i_p,:);
            for i_m = 1 : 1 : size(M,1)
                m = M(i_m,:);
                
                Rp = rs-rm + 2*(p.*rm) ;
                Rm = 2*(m.*dim);
                
                toas = [toas norm((Rp+Rm)/sound_speed)];
                labels = [labels ; [abs(m(1)+p(1)),abs(m(1)),abs(m(2)+p(2)),abs(m(2)),abs(m(3)+p(3)),abs(m(3))]];
            end
        end

        if (m_ == 1) & (s == 1)
            
            reflection_orders = sum(labels,2);
            mask = (reflection_orders <= Q_max);
        end
        
        toas = toas(mask);
        labels = labels(mask,:);
        all_toas = [all_toas toas];
    end
end

K = labels;

end

