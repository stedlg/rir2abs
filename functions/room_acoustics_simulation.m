function [Roomsim,para] = room_acoustics_simulation(n_rooms,fs,T,max_order,n_m_s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rng('shuffle')

% simulator parameters
para = readsetup('sampleroomsetup');
noctave = numel(para.room.surface.frequency);
para.options.fs = fs;
para.options.responseduration = T;
para.options.reflectionorder = [max_order max_order max_order];        % maximum specular reflection order (x,y,z)
para.options.numberofrays = 50000;                 % number of rays in simulation (20*K^2)
para.options.diffusetimestep  = 0.01;           % time resolution in diffuse energy histogram (seconds)
para.options.airabsorption = false;
para.room.surface.diffusion = zeros([6 6]);

% absorption coefficients ranges 
abs_coeffs_lower_bound = [[0.01,0.01,0.01,0.01,0.01,0.01];
    [0.01,0.01,0.01,0.01,0.01,0.01];
    [0.01,0.01,0.01,0.01,0.01,0.01];
    [0.01,0.01,0.01,0.01,0.01,0.01];
    [0.01,0.01,0.05,0.15,0.25,0.30];
    [0.01,0.15,0.40,0.40,0.40,0.30]];

abs_coeffs_upper_bound = [[0.50,0.50,0.30,0.12,0.12,0.12];
    [0.50,0.50,0.30,0.12,0.12,0.12];
    [0.50,0.50,0.30,0.12,0.12,0.12];
    [0.50,0.50,0.30,0.12,0.12,0.12];
    [0.20,0.30,0.50,0.60,0.75,0.80];
    [0.70,1.00,1.00,1.00,1.00,1.00]];

for i = 1 : 1 : n_rooms
    
    % sample surface absorption profiles (see paper)
    j = randi(6,1,6);
    coin_flip = rand(6,1);
    P = (coin_flip >= 0.5).*(abs_coeffs_lower_bound(j,:)+(abs_coeffs_upper_bound(j,:)-abs_coeffs_lower_bound(j,:)).*rand(6,noctave))+(coin_flip < 0.5).*ones(6,noctave).*((0.12-0.01)*rand(6,1)+ 0.01);
    para.room.surface.absorption = P;
    % sample room dimensions
    para.room.dimension(1) = 3 + rand()*(10-3) ;
    para.room.dimension(2) = 3 + rand()*(10-3) ;
    para.room.dimension(3) = 2 + rand()*(5-2) ;
    
    % sample sources and receivers locations 
    
    for n = 1:1:n_m_s
        
        invalid_positions = true ;
        
        while invalid_positions
            r = rand() < 0.5 ;
            if r
                
                para.receiver(n).location = 1 + ([para.room.dimension(1) para.room.dimension(2) para.room.dimension(3)] - 2).*rand(1,3);
                para.source(n).location = 1 + ([para.room.dimension(1) para.room.dimension(2) para.room.dimension(3)] - 2).*rand(1,3);
                
            else
                
                para.source(n).location = 1 + ([para.room.dimension(1) para.room.dimension(2) para.room.dimension(3)] - 2).*rand(1,3);
                para.receiver(n).location = 1 + ([para.room.dimension(1) para.room.dimension(2) para.room.dimension(3)] - 2).*rand(1,3);
                
            end
            
            if norm(para.source(n).location - para.receiver(n).location) > 1
                invalid_positions = false ;
            end
        end
        
        para.source(n).description = 'omnidirectional';
        para.source(n).orientation = [ 0 0 0 ];
        
        para.receiver(n).description = 'omnidirectional';
        para.receiver(n).orientation = [ 0 0 0 ];
    end
   
    % save roomsim structure 
    
    Roomsim(i).inputvar.dim.L = para.room.dimension(1);
    Roomsim(i).inputvar.dim.l = para.room.dimension(2);
    Roomsim(i).inputvar.dim.h = para.room.dimension(3);
    
    Roomsim(i).inputvar.absorption = para.room.surface.absorption;
    Roomsim(i).inputvar.diffusion = para.room.surface.diffusion;
    Roomsim(i).inputvar.source = para.source ;
    Roomsim(i).inputvar.receiver = para.receiver;
    
    % generate room impulse responses (RIRs)
    Roomsim(i).output.RIR = roomsim(para);
    
end
end


