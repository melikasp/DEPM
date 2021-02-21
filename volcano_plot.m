function volcano_plot(FC,pvalue,fct,fdrt)

FC_uperT=[];FC_midT=[];FC_downT=[]; 
a=1;b=1;c=1;
for i=1:length(FC)
    n=FC(i);
    if n>fct
        FC_uperT(a,1)=FC(i);
        FC_uperT(a,2)=pvalue(i);
        a=a+1;
    end
    if n<-fct
        FC_downT(b,1)=FC(i);
        FC_downT(b,2)=pvalue(i);
        b=b+1;
    else
        FC_midT(c,1)=FC(i);
        FC_midT(c,2)=pvalue(i);
        c=c+1;
    end
end

plot(FC_midT(:,1),FC_midT(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3);
hold on 
plot(FC_uperT(:,1),FC_uperT(:,2),'o','MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',3);
hold on 
plot(FC_downT(:,1),FC_downT(:,2),'o','MarkerEdgeColor','g','MarkerFaceColor','w','MarkerSize',3);
hold on
grid on;title('volcano plot');xlabel('Folld Change(log2)');ylabel('-log10(pvalue)');
hold on
plot([-fct,-fct],[0,max(pvalue)+1],'r');
hold on
plot([fct,fct],[0,max(pvalue)+1],'r');
hold on
plot([min(FC),max(FC)],[fdrt,fdrt],'r');
end