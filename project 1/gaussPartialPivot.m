function x=gaussPartialPivot(A,b)

    aug=[A,b];
    [mA,nA]=size(A);
    [mB,nB]=size(b);
    x=zeros(length(b),1);
    if nA~=mB
        error('Wrong matrix dimensions');
    end
    
    for i=1:nA
        [value,id] =max(abs(aug(i:end,i)));
        id=id+i-1;
        
        temp=aug(i,:);
        aug(i,:)=aug(id,:);
        aug(id,:)=temp;
        pivot=aug(i,:);
        for j=i+1:nA
           ratio=aug(j,i)/pivot(i);
           aug(j,:)=aug(j,:)-ratio*pivot;
        end
        %rref(aug)
    end
    for i=mA:-1:2
        pivot=aug(i,:);
        for j=i-1:-1:1
            ratio=aug(j,i)/pivot(i);
            aug(j,:)=aug(j,:)-ratio*pivot;
        end
        aug(i,:)=aug(i,:)/aug(i,i);
    end
    x=aug(:,end);
end