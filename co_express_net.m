function [G adjC]=co_express_net(gene_data,P_th)

Corr_data = corr(log(gene_data')); 
Corr_data(abs(Corr_data)<P_th) = 0; 
Corr_data = Corr_data-triu(tril(Corr_data));
adjC = logical(Corr_data);
G = graph(adjC);
figure;h = plot(G,'Layout','force','UseGravity',true);
title({'Differentially expressed genes';'Co-expression network in Normal cells'});

end