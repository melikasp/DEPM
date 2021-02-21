load('FinalData.mat')

%% parameters
P_th = 0.6; % Pearson Correlation Coefficient (P<0.6) threshold
fct=1.2;
fdrt=0.05;

%% 2-volcano plot

FC = mean(log2(dataC'./ dataN'))'; % Average Flood Change
[h,p] = ttest2(dataC',dataN'); % p-value
pvalue=-(log10(p));

%call volcano plot function
volcano_plot(FC,pvalue,fct,fdrt);
clear h p pvalue fct fdrt

%% 3-a)compute the gene co-expression networks related to the 2 conditions (cancer, normal) 
% considering: Pearson’s correlation and Binary adjacency matrix

% The Normal gene co-expression networks
[G_N,adjN]=co_express_net(dataN,P_th);

% the cancer gene co-expression networks
[G_C adjC]=co_express_net(dataC,P_th);

%% 3-b) Analysis: Compute the degree index of Normal and Cancer gene

%Normal
degree_N = centrality(G_N,'degree');
degree_N(degree_N==0) = [];
figure;histogram(degree_N) 
title('Degree of Normal distribution')
xlabel('k');ylabel('Occurrences')

% Cancer
degree_CC = centrality(G_C,'degree');
degree_CC(degree_CC==0) = [];
figure;histogram(degree_CC) ;
title('Degree of Cancer distribution');
xlabel('k');ylabel('Occurrences');
clear degree_N degree_C

%% Centrality Measures: normal data

HubsN_Degree=CenterM(adjN,'degree',geneID);

HubsN_closeness=CenterM(adjN,'closeness',geneID);

HubsN_betweenness=CenterM(adjN,'betweenness',geneID);

HubsN_eigenvector=CenterM(adjN,'eigenvector',geneID);

%% Centrality Measures: cancer data

HubsC_Degree=CenterM(adjC,'degree',geneID);

HubsC_closeness=CenterM(adjC,'closeness',geneID);

HubsC_betweenness=CenterM(adjC,'betweenness',geneID);

HubsC_eigenvector=CenterM(adjC,'eigenvector',geneID);
%% compare Hubs of Normal and Cancer

mutual_gene_NC = intersect(cell2mat(HubsC_Degree),cell2mat(HubsN_Degree));

%% 4-Differential Co-expressed Network: Computation

CC = corr(log(dataC')); 
CC = CC-triu(tril(CC));
Z1=atanh(CC);
clear CC

CN = corr(log(dataN')); 
CN = CN-triu(tril(CN));
Z2=atanh(CN);
clear CN

Z = (Z1-Z2)/sqrt(1/(length(dataC)-3)+1/(length(dataC)-3));
clear Z2 Z1
Z(abs(Z)<3) = 0;
adjZ = logical(Z);
G_Z = graph(adjZ);

degree_Z = centrality(G_Z,'degree');
degree_Z(degree_Z==0) = [];
figure;histogram(degree_Z) ;
title('Degree of distribution')

%% 4-Differential Co-expressed Network:Analysis:

HubsZ_Degree=CenterM(adjZ,'degree',geneID);


%% find upregulated and downregulated

CentrMeasures = 'degree';
index = centrality(G_Z,CentrMeasures);
Y = prctile(index,99);
Hubs_ind = find(index>Y);
HubsFC = FC(Hubs_ind);
figure; h3 = plot(G_Z,'Layout','force','UseGravity',true);
highlight(h3,Hubs_ind,'NodeColor','r','Marker','h','MarkerSize',4)
if find(HubsFC<0)
    highlight(h3,Hubs_ind((HubsFC<0)),'NodeColor','cyan','Marker','h','MarkerSize',4)
end

title(CentrMeasures)

downregulated=geneID(HubsFC<0);
upregulated=geneID(HubsFC>0);

%% compare part3 and part4
mutual_gene_ZN = intersect(cell2mat(HubsN_Degree),cell2mat(HubsZ_Degree));
mutual_gene_ZC = intersect(cell2mat(HubsC_Degree),cell2mat(HubsZ_Degree));









