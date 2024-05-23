function Data_o = AOAA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Ffun=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
M_Iter=Data_i.maxIter;  LB=Data_i.ub; Dim=Data_i.dim;
F_obj=Data_i.fobj;   UB=Data_i.ub;   LB=Data_i.lb;

Best_P=zeros(1,Dim);
Best_FF=inf;
IterCurve=zeros(1,Data_i.maxIter);

%Initialize the positions of solution
Xnew=X;

Ffun_new=Ffun;% (fitness values)

MOP_Max=1;
MOP_Min=0.2;
C_Iter=1;
Alpha=5;
Mu=0.499;

for i=1:size(X,1)
    if Ffun(i)<Best_FF
        Best_FF=Ffun(i);
        Best_P=X(i,:);
    end
end
    
while C_Iter<M_Iter+1  %Main loop
    MOP=1-((C_Iter)^(1/Alpha)/(M_Iter)^(1/Alpha));   % Probability Ratio 
    MOA=MOP_Min+C_Iter*((MOP_Max-MOP_Min)/M_Iter); %Accelerated function
   
    %Update the Position of solutions
    for i=1:size(X,1)   % if each of the UB and LB has a just value 
        for j=1:size(X,2)
           r1=rand();
            if (size(LB,2)==1)
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB-LB)*Mu+LB);
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB-LB)*Mu+LB);
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB-LB)*Mu+LB);
                    end
                end               
            end
            
           
            if (size(LB,2)~=1)   % if each of the UB and LB has more than one value 
                r1=rand();
                if r1<MOA
                    r2=rand();
                    if r2>0.5
                        Xnew(i,j)=Best_P(1,j)/(MOP+eps)*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)*MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                else
                    r3=rand();
                    if r3>0.5
                        Xnew(i,j)=Best_P(1,j)-MOP*((UB(j)-LB(j))*Mu+LB(j));
                    else
                        Xnew(i,j)=Best_P(1,j)+MOP*((UB(j)-LB(j))*Mu+LB(j));
                    end
                end               
            end
            
        end
        
        Flag_UB=Xnew(i,:)>UB; % check if they exceed (up) the boundaries
        Flag_LB=Xnew(i,:)<LB; % check if they exceed (down) the boundaries
        Xnew(i,:)=(Xnew(i,:).*(~(Flag_UB+Flag_LB)))+UB.*Flag_UB+LB.*Flag_LB;
 
        Ffun_new(i)=F_obj(Xnew(i,:));  % calculate Fitness function 
        if Ffun_new(i)<Ffun(i)
            X(i,:)=Xnew(i,:);
            Ffun(i)=Ffun_new(i);
        end
        if Ffun(i)<Best_FF
            Best_FF=Ffun(i);
            Best_P=X(i,:);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Best_FF;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter);
        end
    end
    %Update the convergence curve
    IterCurve(C_Iter)=Best_FF;
    
    %Print the best solution details after every 50 iterations     
    C_Iter=C_Iter+1;  % incremental iteration
end   
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=C_Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end