function x = doolittleMethod(A,b)
[m,n]=size(A)
if m ~= n
    error("must be square matrix")
end
L=eye(m);
U=A;

for i=1:m
        pivot=U(i,:);
        for j=i+1:m
            ratio=U(j,i)/pivot(i);
            if ratio ~=0
                U(j,:)=U(j,:)-ratio*pivot;
                L(j,i)=ratio;
            end
        end
    end

L
U
if L*U ~= A
    error("Function broken. L*U != A, but it should."); 
end

y=gaussPartialPivot(L,b)
x=gaussPartialPivot(U,y)
end