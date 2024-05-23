function Data_o=TSO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 

%% Problem Definition
Max_iter = Data_i.maxIter;      
Particles_no = Data_i.pop;      
Dim = Data_i.dim;              
T = Data_i.X;                   
Tuna1_fit = Data_i.F_value;     
IterCurve=zeros(1,Data_i.maxIter);
[fvalbest,index]=min(Tuna1_fit);
Tuna1=T(index,:);              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Tuna1_fit);

Iter=0;
aa=0.7;
z=0.05;
while Iter<Max_iter
    C=Iter/Max_iter;
    a1=aa+(1-aa)*C;
%     rand_num(1,Iter+1) = a1;
    a2=(1-aa)-(1-aa)*C;
%     rand_num(2,Iter+1) = a2;

    %---------------- Memory saving-------------------
    if Iter==0
        Tuna2=T;
    end

    for i=1:Particles_no
        if Tuna1_fit(i)<fvalbest
            fvalbest=Tuna1_fit(i); 
            Tuna1=T(i,:);
        end
    end

    %-------------------------------------------------
    
    t=(1-Iter/Max_iter)^(Iter/Max_iter);                   
    
    if rand<z
        T(1,:)= (Data_i.ub-Data_i.lb)*rand+Data_i.lb;
    else
        if  0.5<rand
            r1=rand;
            Beta=exp(r1*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*r1));
            if  C>rand
                T(1,:)=a1.*(Tuna1+Beta*abs(Tuna1-T(1,:)))+a2.*T(1,:); %Equation (8.3)
                
            else
                IndivRand=rand(1,Dim).*(Data_i.ub-Data_i.lb)+Data_i.lb;
                T(1,:)=a1.*(IndivRand+Beta*abs(IndivRand-T(1,:)))+a2.*T(1,:);%Equation (8.1)
            end
        else
            TF = (rand>0.5)*2-1;
            if 0.5>rand
                T(1,:)=Tuna1+rand(1,Dim).*(Tuna1-T(1,:))+TF.*t^2.*(Tuna1-T(1,:));%Equation (9.1)
            else
                T(1,:) =TF.* t^2.*T(1,:);%Equation (9.2)
            end
            
        end
        
    end
    
    for i=2:Particles_no
        if rand<z    
            T(i,:)= (Data_i.ub-Data_i.lb)*rand+Data_i.lb;
        else
            if  0.5<rand
                r1=rand;
                Beta=exp(r1*exp(3*cos(pi*((Max_iter-Iter+1)/Max_iter))))*(cos(2*pi*r1));
                if  C>rand
                    T(i,:)=a1.*(Tuna1+Beta*abs(Tuna1-T(i,:)))+a2.*T(i-1,:);%Equation (8.4)
                else
                    
                    IndivRand=rand(1,Dim).*(Data_i.ub-Data_i.lb)+Data_i.lb;
                    T(i,:)=a1.*(IndivRand+Beta*abs(IndivRand-T(i,:)))+a2.*T(i-1,:);%Equation (8.2)
                end
            else
                TF = (rand>0.5)*2-1;
                if 0.5>rand
                    T(i,:)=Tuna1+rand(1,Dim).*(Tuna1-T(i,:))+TF*t^2.*(Tuna1-T(i,:)); %Equation (9.1)
                else
                    T(i,:) = TF*t^2.*T(i,:);%Equation (9.2)
                end
            end
        end
    end
    
    for i=1:size(T,1)
        T=BoundaryCheck(T,Data_i.ub,Data_i.lb,Data_i.dim);
        fitness(i)=Data_i.fobj(T(i,:));
        if fitness(i)<Tuna1_fit(i)
            Tuna1_fit(i)=fitness(i);  
            Tuna2(i,:)=T(i,:);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    
    Iter=Iter+1;
    fval = min(Tuna1_fit);
    IterCurve(Iter) = fval;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = fval;
    
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end



