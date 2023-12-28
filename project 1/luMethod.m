function x=luMethod(A,b)

[L,U,P]=lu(A);


z=P.'\b;




aug1=[L z];
[mL,nL]=size(A);
[mB,nB]=size(b);
for i=1:mL
        pivot=aug1(i,:);
        for j=i+1:mL
            ratio=aug1(j,i)/pivot(i);
            aug1(j,:)=aug1(j,:)-ratio*pivot;
        end
end

[mA,nA]=size(U);

aug2=[U aug1(:,end)];
aug2(1,:)=aug2(1,:)/aug2(1,1);
for i=mA:-1:2
        pivot=aug2(i,:);
        for j=i-1:-1:1
            ratio=aug2(j,i)/pivot(i);
            aug2(j,:)=aug2(j,:)-ratio*pivot;
        end
        aug2(i,:)=aug2(i,:)/aug2(i,i);
end

x=aug2(:,end);


end