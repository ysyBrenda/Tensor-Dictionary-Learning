function [AUC,Zauc] = computeAUC(xi,yi)
%compute AUC£¬ xi= positive sample£¨is ore£©£¬yi=negative sample£¨not ore£©
q=length(yi);  %the number of TF£¬TP
p=length(xi);
sum=0;
for i=1:length(xi)
    sum_a=0;
    for j=1:length(yi)
        if (xi(i)-yi(j))>10e-6
            a=1;
        elseif abs(xi(i)-yi(j))<10e-6
            a=0.5;
        else
            a=0;
        end
        sum_a=sum_a+a;
    end
    sum=sum_a+sum;
end
AUC=sum/(q*p);
SEauc=sqrt((AUC*(1-AUC)+(p-1)*(AUC/(2-AUC)-AUC*AUC)+(q-1)*(2*AUC*AUC/(1+AUC)-AUC*AUC))/(p*q));
Zauc=(AUC-0.5)/SEauc;
end

