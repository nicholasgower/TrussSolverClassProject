function x=naiveGauss(A,b)

    aug=[A,b];
    [mA,nA]=size(A);
    [mB,nB]=size(b);
    x=zeros(length(b),1);
    if nA~=mB
        error('Wrong matrix dimensions');
    end

    for i=1:mA
        pivot=aug(i,:);
        for j=i+1:mA
            ratio=aug(j,i)/pivot(i);
            aug(j,:)=aug(j,:)-ratio*pivot;
        end
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
