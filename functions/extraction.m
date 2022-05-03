function [dataStructure] = extraction(data)
%EXTRACTION Summary of this function goes here

n_samples = size(data);
n_samples = n_samples(1,2);

for n = [1 : 1 : n_samples]
    dataStructure(n).RIR = data(n).output.RIR ;
    abs = data(n).inputvar.absorption.';
    diff = data(n).inputvar.diffusion.';
    dim = struct2array(data(n).inputvar.dim);
    dataStructure(n).receiver = data(n).inputvar.receiver;
    dataStructure(n).source = data(n).inputvar.source;
    dataStructure(n).receiver = data(n).inputvar.receiver;
    dataStructure(n).dim = data(n).inputvar.dim;
    dataStructure(n).hpfilter = data(n).inputvar.hpfilter;
    dataStructure(n).stdev = data(n).inputvar.stdev;
    dataStructure(n).absorption = abs(:);
    dataStructure(n).diffusion = diff(:,1);
    dataStructure(n).dim = dim; % L l h 
    dataStructure(n).G_noise = data(n).inputvar.G_noise;
    %dataStructure(n).exact = data(n).inputvar.exact;
end 


end

