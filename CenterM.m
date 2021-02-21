function Hubs=CenterM(adjm,cm,geneID)
K = sum(adjm,2);
adjm((K==0),:) = []; adjm(:,(K==0)) = [];
G = graph(adjm);
geneID((K==0)) = [];


index = centrality(G,cm);
Y = prctile(index,95);
figure; h2 = plot(G,'Layout','force','UseGravity',true);
highlight(h2,(index>Y),'NodeColor','r','Marker','h','MarkerSize',4)
title(cm)
Hubs= geneID(index>Y);

end