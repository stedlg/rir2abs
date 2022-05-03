%{This script perform the ransac approach described in the paper to estimate the surface absorption profiles %}
%% clean up workspace and set folder path
clear;
close all ;
clc;
path_package = '~/Documents/these/geometry-informed-estimation/';
path_roomsim = '~/Documents/these/matlab/roomsim-1/';
cd(path_package)
addpath(path_package)
addpath(strcat(path_package,'functions'))
addpath(strcat(path_package,'data'))
addpath(strcat(path_package,'results'))
addpath(path_roomsim)
%% simulate data (with ROOMSIM shoebox room acoustics simulator)
simulate_data = false;
n_batchs = 1; n_rooms = 5; fs = 16000; T= 0.25; max_order = 20 ; n_m_s = 5 ; % n_m_s : number of microphones and sources 
if simulate_data
    for i = 1:1: n_batchs
        roomsim = room_acoustics_simulation(n_rooms,fs,T,max_order,n_m_s);
        save(strcat(path_package,'data/roomsim_',num2str(i),'.mat'),"roomsim");
    end
end
clearvars -except simulate_data path_package 
%% concatenate data
if ~simulate_data
    if isfile(strcat(path_package,'room_acoustics_dataset.mat'))
        disp('loading data');
        load(strcat(path_package,'room_acoustics_dataset.mat'));
    else
        disp('please set simulate_data = true;');
        return
    end
else
    disp('concatenate data');
    path_data = strcat(path_package,'data');
    listing = dir(path_data);
    n_files = size(listing);
    
    for i = 1:1:n_files
        if listing(i).isdir
            continue
        end
        file = strcat(listing(i).folder,"/",listing(i).name);
        load(file);
        
        if exist('data') ~= 1
            data = roomsim;
        else
            data = [data roomsim];
        end
    end
    
    data = extraction(data);
    save(strcat(path_package,'room_acoustics_dataset.mat'),'data');
end
dataset = data;
clear data
%% set main parameters
fs = 16000;
sound_speed = 3.424396008126316e+02;

Q_max = 2;
Q = 2;

n_sources = 3;
n_receivers = 3;

size_rir = size(dataset(1).RIR{1,1},1); % dataset(i).RIR{s, m}

% ism
Nm = 1;

% distortions
noise = true; % awgn
psnr_ = 50;
sigma_geo = 0.02;

% spectrograms
w_ms = 2;
L_v = 0.5*10^(-3);
window_size = w_ms *(fs/1000);
fft_length = window_size;
overlap = window_size-1;
r = 0; % 0 : rectangle / 1 : hanning
f_0 = 0;
f_max = fs/2 ;
fqs = round((fs/2)/((fft_length/2)+1).*[0:((fft_length/2)+1)-1]);
f_min_idx = min(find(fqs>=f_0));
f_max_idx = min([find(fqs>f_max),size(fqs,2)]);

% ransac
N_iter = 1000;
threshold_ransac = 0.1;
threshold_Pi = 0.001;
natural_score = false;

% visualisation / metrics
graphic_evaluation = true;
wall_labels = ["west wall", "east wall","south wall","north wall","floor","ceiling"];
threshold_error = 0.1;

warning('off','all')
%% loop over rooms
clc;

size_data = size(dataset,2);

disp(strcat("estimation of surface absorption profiles over the " ,num2str(size_data)," rooms "));

targets_abs_coeffs = [];
estimates_abs_coeffs = [];
not_estimated_abs_coeffs = [];

for i_room = 1 : 1 : size_data
    
    data = dataset(i_room);
    
    [targets_abs_coeffs_i,estimates_abs_coeffs_i,not_estimated_abs_coeffs_i] = estimate_absorption_profiles(data,fs,sound_speed,Q_max,Q,n_sources,n_receivers,size_rir,Nm,noise,psnr_,sigma_geo,L_v,window_size,overlap,fft_length,r,fqs,f_min_idx,f_max_idx,N_iter,threshold_ransac,threshold_Pi,natural_score,graphic_evaluation,wall_labels);
    
    targets_abs_coeffs = cat(3,targets_abs_coeffs,targets_abs_coeffs_i);
    estimates_abs_coeffs = cat(3,estimates_abs_coeffs,estimates_abs_coeffs_i);
    not_estimated_abs_coeffs = cat(3,not_estimated_abs_coeffs,not_estimated_abs_coeffs_i);
    
end

clc
disp('finished')
%% compute metrics
% Mean absolute error (MAE), percentage of correct estimates (CE)
% /!\ first frequency bin excluded (see paper)
abs_error = abs(estimates_abs_coeffs(:,2:end,:) - targets_abs_coeffs(:,2:end,:));
abs_error(not_estimated_abs_coeffs(:,2:end,:)==1) = 0.5; % see paper
mae = mean(abs_error,[1 2 3]);
ce = mean(abs_error < threshold_error,[1 2 3]);
save(strcat(path_package,"results/",'results.mat'),'targets_abs_coeffs','estimates_abs_coeffs','not_estimated_abs_coeffs','psnr_','sigma_geo','n_sources','n_receivers');