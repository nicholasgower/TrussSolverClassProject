function x = cramersRule(A,b)

    [mA,nA]=size(A);
    [mB,nB]=size(b);
    x=zeros(length(b),1);
    if nA~=mB
        error('Wrong matrix dimensions');
    end
    for n=1:nA
        aMod=A;
        aMod(:,n)=b;
        
        x(n)=det(aMod)/det(A);
    end
    

end