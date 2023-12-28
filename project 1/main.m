% Numerical Truss Solver
% By Nicholas Gower
% 02/12/21

% Based on problem 6-4 in Engineering Mechanics: Statics 14th Ed.
% "Determine the force in each member of the truss
% and state if the members are in tension or compression."

% The steps for solving the problem:
% 1. Declare origin and coordinate system. Right=+x, Up=+y, A=(0,0)
% 2. Describe the truss, writing the position and connections of each pin.
% 3. Describe all non-member reaction and applied forces
% 4. For each pin, write two equations that would be made from a FBD
% Note: Distances are measured in feet, and forces are measured in kip.
% 5. With this information, generate augmented matrix that describes the
% entire system.
% 6. With various methods, row reduce the matrix.

% For every member, there is a tension force parallel to the member
saving=0; %Change to 1 to enable writing to disk
pinsX=[0,10,20,30,40, 30,20,10,10 ]; %Positions of each pin
pinsY=[0,0,0,0,0,     4  ,8   ,12 ,4];

% Pin coordinate data can be tested with scatter(pinsX,pinsY).

% Pin A=1
% Pin B=2
% Pin C=3
% Pin D=4
% Pin E=5
% Pin F=6
% Pin G=7
% Pin H=8
% Pin I=9
connection=zeros(9,5);

connection(1,:)=[2,9,0,0,0]; % Pin A is connected with pins B and I
connection(2,:)=[1,3,9,0,0]; % Pin B is connected with pins A, C, and I
connection(3,:)=[2,9,7,6,4]; % etc.
connection(4,:)=[3,6,5,0,0];
connection(5,:)=[4,6,0,0,0];
connection(6,:)=[7,3,4,5,0];
connection(7,:)=[8,9,3,6,0];
connection(8,:)=[9,7,0,0,0];
connection(9,:)=[1,2,3,7,8];
arraySize=size(connection);
numMembers=0;
for i=1:arraySize(1)
    for j=1:arraySize(2)
        if connection(i,j) ~=0
            numMembers=numMembers+0.5; %For every two pin-pin connections, there is one member.
        end
    end
end

memberForceIDsystem=zeros(numMembers); %Enter two pins, get column of force in forceMatrix
forceNames=["A","Ex","Ey"];
numNonMemberForces=length(forceNames); 



forceMatrix=zeros(2*length(pinsX),numMembers+numNonMemberForces+1); %Augmented matrix
%Will be 18 x 19 matrix, so problem is solvable.
numForces=numMembers+numNonMemberForces;

%sum(forces)
% This loop generates the forceMatrix. 
for i=1:length(pinsX) %For every pin in the truss
    pin=[pinsX(i),pinsY(i)];
    forces=zeros(2,numForces+1);
    connections=connection(i,:);
    for n=1:length(connections) %Adds tension forces
        if connections(n)~=0
            otherPin=[pinsX(connections(n)),pinsY(connections(n))];
            displacement=otherPin-pin;
            magnitude=getMagnitude(displacement);
            
            xAxis=[1,0]; 
            angle=acosd(dot(displacement,xAxis)/(magnitude)); %Finds angle of tension vector using dot product
            if displacement(2)<0 %Flips sign of angle if y component is negative
                angle=-angle;
            end
            
            pinIDs=[i,connections(n)];
            %forceID=(length(pinsX)+1)*min(pinIDs)+max(pinIDs);
            firstPin=min(pinIDs);
            secondPin=max(pinIDs);
            forceID=memberForceIDsystem(firstPin,secondPin); % Gets force ID from array
            if forceID==0 %If force ID not set yet
               forceID=length(forceNames)+1; %Generate new ID
               memberForceIDsystem(firstPin,secondPin)=forceID; %Send new ID to array, so 3rd law pair can use the same ID
               forceName=strcat('T_',number2letter(firstPin),number2letter(secondPin)); %Generate new force name
               forceNames(forceID)=forceName; %Send name of force to name list
            end
            
            forces(1,forceID)=cosd(angle); %Add tension force with angle
            forces(2,forceID)=sind(angle);
            
        end
    end
    
    %Applied forces and non-tension reaction forces
    
    if i==1 %Pin A
        forces(2,1)=1; % Rocker A reaction force
    elseif i==5 %Pin E
        forces(1,2)=1; % Fixed pin E reaction force
        forces(2,3)=1;
    elseif i==6 
        forces(2,length(forces))=3; % -3 kip applied force down
    elseif i==7
        forces(2,length(forces))=3; % -3 kip applied force down
    elseif i==8
        forces(1,length(forces))=-2; % 2 kip applied force right
    elseif i==9
        forces(1,length(forces))=-1.5; % 1.5 kip applied force right
    end
    
    
    forceMatrix(i*2-1,:)=forces(1,:); %Adds force equations to main matrix
    forceMatrix(i*2,:)=forces(2,:);
end

%save data/forceNames forceNames 
if saving
    for i=1:length(forceNames)
        forceNames(i)=strcat(num2str(i)," : ",forceNames(i));
    end
    writematrix(forceNames.',"data/forceNames")
    save data/forceMatrix.dat forceMatrix -ASCII
end


A=forceMatrix(:,1:end-1);
b=forceMatrix(:,end);
%A*x=b

% custom Cramer's Rule
% Time: 470 mu-s
cramer=cramersRule(A,b);
if saving
    save data/customCramer.dat cramer -ASCII
end
% custom naive Gaussian elimination
% Time: 174 mu-s
% ---This method does not work with this problem.
gNaive=naiveGauss(A,b);
if saving
    save data/customNaive.dat gNaive -ASCII
end

% custom Gaussian elimination with partial pivoting
% Time: 329 mu-s
gPartial=gaussPartialPivot(A,b);
if saving
    save data/customPartial.dat gPartial -ASCII
end

% custom LU method

xLu=luMethod(A,b);
save data/luMethod.dat xLu -ASCII

% built-in inverse matrix multiplication (A\b)
% Time: 97 mu-s
inverse=inverseWrapper(A,b); 
if saving
    save data/builtInInverse.dat inverse -ASCII
end

% built-in inverse matrix multiplication (A^-1*b)
inverse2=inverseWrapper2(A,b); 
if saving
    save data/builtInInverse2.dat inverse2 -ASCII
end
% built-in inverse matrix multiplication (A^-1*b)
inverse3=inverseWrapper3(A,b); 
if saving
    save data/builtInInverse3.dat inverse3 -ASCII
end
% built-in rref function
% Time 3502 mu-s
builtInGauss=rref(forceMatrix);
values=builtInGauss(:,end);
if saving
    save data/builtInGauss.dat values -ASCII
end

