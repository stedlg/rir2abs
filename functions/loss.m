function [val] = loss(model,X_J,H_J,K) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% [d,w_1,...,w_6] estimation
val =  sum((X_J - H_J.*(model(1)*prod(model(2:end)'.^K,2))).^2,'all');
end
