function Data_o = CHOA(Data_i,Data_o)
rand_num=[];                        
ti=clock;                              
Positions = Data_i.X;                        
Cost=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);  
Max_iter=Data_i.maxIter;        ub=Data_i.ub;           lb=Data_i.lb;       fobj=Data_i.fobj;
dim=Data_i.dim;
% initialize Attacker, Barrier, Chaser, and Driver
Attacker_pos=zeros(1,dim);
Attacker_score=inf; %change this to -inf for maximization problems

Barrier_pos=zeros(1,dim);
Barrier_score=inf; %change this to -inf for maximization problems

Chaser_pos=zeros(1,dim);
Chaser_score=inf; %change this to -inf for maximization problems

Driver_pos=zeros(1,dim);
Driver_score=inf; %change this to -inf for maximization problems
l=1;
while l<=Max_iter
    for i=1:size(Positions,1)  
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness=fobj(Positions(i,:));

        if fitness<Attacker_score 
            Attacker_score=fitness; % Update Attacker
            Attacker_pos=Positions(i,:);
        end

        if fitness>Attacker_score && fitness<Barrier_score 
            Barrier_score=fitness; % Update Barrier
            Barrier_pos=Positions(i,:);
        end
        
        if fitness>Attacker_score && fitness>Barrier_score && fitness<Chaser_score 
            Chaser_score=fitness; % Update Chaser
            Chaser_pos=Positions(i,:);
        end

        if fitness>Attacker_score && fitness>Barrier_score && fitness>Chaser_score && fitness>Driver_score 
            Driver_score=fitness; % Update Driver
            Driver_pos=Positions(i,:);
        end

        %% 更新全局最优值和最优位置
        if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Attacker_score
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Attacker_score;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    f=2-l*((2)/Max_iter);

    %Group 1
    C1G1=1.95-((2*l^(1/3))/(Max_iter^(1/3)));
    C2G1=(2*l^(1/3))/(Max_iter^(1/3))+0.5;
        
    %Group 2
    C1G2= 1.95-((2*l^(1/3))/(Max_iter^(1/3)));
    C2G2=(2*(l^3)/(Max_iter^3))+0.5;
    
    %Group 3
    C1G3=(-2*(l^3)/(Max_iter^3))+2.5;
    C2G3=(2*l^(1/3))/(Max_iter^(1/3))+0.5;
    
    %Group 4
    C1G4=(-2*(l^3)/(Max_iter^3))+2.5;
    C2G4=(2*(l^3)/(Max_iter^3))+0.5;

    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
            r11=C1G1*rand(); % r1 is a random number in [0,1]
            r12=C2G1*rand(); % r2 is a random number in [0,1]
            
            r21=C1G2*rand(); % r1 is a random number in [0,1]
            r22=C2G2*rand(); % r2 is a random number in [0,1]
            
            r31=C1G3*rand(); % r1 is a random number in [0,1]
            r32=C2G3*rand(); % r2 is a random number in [0,1]
            
            r41=C1G4*rand(); % r1 is a random number in [0,1]
            r42=C2G4*rand(); % r2 is a random number in [0,1]
            
            A1=2*f*r11-f; % Equation (3)
            C1=2*r12; % Equation (4)

            m=chaos(3,1,1); % Equation (5)
            D_Attacker=abs(C1*Attacker_pos(j)-m*Positions(i,j)); % Equation (6)
            X1=Attacker_pos(j)-A1*D_Attacker; % Equation (7)
                       
            A2=2*f*r21-f; % Equation (3)
            C2=2*r22; % Equation (4)
      
            D_Barrier=abs(C2*Barrier_pos(j)-m*Positions(i,j)); % Equation (6)
            X2=Barrier_pos(j)-A2*D_Barrier; % Equation (7)     

            A3=2*f*r31-f; % Equation (3)
            C3=2*r32; % Equation (4)
            
            D_Driver=abs(C3*Chaser_pos(j)-m*Positions(i,j)); % Equation (6)
            X3=Chaser_pos(j)-A3*D_Driver; % Equation (7)      
            
            A4=2*f*r41-f; % Equation (3)
            C4=2*r42; % Equation (4)
            
            D_Driver=abs(C4*Driver_pos(j)-m*Positions(i,j)); % Equation (6)
            X4=Chaser_pos(j)-A4*D_Driver; % Equation (7)       
            Positions(i,j)=(X1+X2+X3+X4)/4;% Equation (8)
        end
    end
    IterCurve(l)=Attacker_score;

    l=l+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);            
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function O=chaos(index,max_iter,Value)
    O=zeros(1,max_iter);
    x(1)=0.7;
    switch index
    %Chebyshev map
        case 1
    for i=1:max_iter
        x(i+1)=cos(i*acos(x(i)));
        G(i)=((x(i)+1)*Value)/2;
    end
        case 2
    %Circle map
    a=0.5;
    b=0.2;
    for i=1:max_iter
        x(i+1)=mod(x(i)+b-(a/(2*pi))*sin(2*pi*x(i)),1);
        G(i)=x(i)*Value;
    end
        case 3
    %Gauss/mouse map
    for i=1:max_iter
        if x(i)==0
            x(i+1)=0;
        else
            x(i+1)=mod(1/x(i),1);
        end
        G(i)=x(i)*Value;
    end
    
        case 4
    %Iterative map
    a=0.7;
    for i=1:max_iter
        x(i+1)=sin((a*pi)/x(i));
        G(i)=((x(i)+1)*Value)/2;
    end
    
        case 5
    %Logistic map
    a=4;
    for i=1:max_iter
        x(i+1)=a*x(i)*(1-x(i));
        G(i)=x(i)*Value;
    end
        case 6
    %Piecewise map
    P=0.4;
    for i=1:max_iter
        if x(i)>=0 && x(i)<P
            x(i+1)=x(i)/P;
        end
        if x(i)>=P && x(i)<0.5
            x(i+1)=(x(i)-P)/(0.5-P);
        end
        if x(i)>=0.5 && x(i)<1-P
            x(i+1)=(1-P-x(i))/(0.5-P);
        end
        if x(i)>=1-P && x(i)<1
            x(i+1)=(1-x(i))/P;
        end    
        G(i)=x(i)*Value;
    end
    
        case 7
    %Sine map
    for i=1:max_iter
         x(i+1) = sin(pi*x(i));
         G(i)=(x(i))*Value;
     end
        case 8
     %Singer map 
     u=1.07;
     for i=1:max_iter
         x(i+1) = u*(7.86*x(i)-23.31*(x(i)^2)+28.75*(x(i)^3)-13.302875*(x(i)^4));
         G(i)=(x(i))*Value;
     end
        case 9
    %Sinusoidal map
     for i=1:max_iter
         x(i+1) = 2.3*x(i)^2*sin(pi*x(i));
         G(i)=(x(i))*Value;
     end
     
        case 10
     %Tent map
     x(1)=0.6;
     for i=1:max_iter
         if x(i)<0.7
             x(i+1)=x(i)/0.7;
         end
         if x(i)>=0.7
             x(i+1)=(10/3)*(1-x(i));
         end
         G(i)=(x(i))*Value;
     end
    
    end
    O=G;
end