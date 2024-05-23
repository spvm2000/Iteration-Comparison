function Data_o = ASO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                            
Atom_Pop = Data_i.X;                        
Fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
Atom_Num=Data_i.pop; alpha=50;  beta=0.2;   Low=Data_i.lb;  Up=Data_i.ub;  Dim=Data_i.dim;
Iteration=1;
Atom_V=rand(Atom_Num,Dim).*(Up-Low)+Low;      Max_Iteration=Data_i.maxIter;
Functon_Best=zeros(Max_Iteration,1);
[Max_Fitness,Index]=min(Fitness);
Functon_Best(1)=Fitness(Index);
X_Best=Atom_Pop(Index,:);
IterCurve(1)=Max_Fitness;
Atom_Acc=Acceleration(Atom_Pop,Fitness,Iteration,Max_Iteration,Dim,Atom_Num,X_Best,alpha,beta);
for Iteration=2:Max_Iteration
    Functon_Best(Iteration)=Functon_Best(Iteration-1);
    Atom_V=rand(Atom_Num,Dim).*Atom_V+Atom_Acc;
    Atom_Pop=Atom_Pop+Atom_V;
    for i=1:Atom_Num
        Atom_Pop(i,:)=BoundaryCheck(Atom_Pop(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
        Fitness(i)=Data_i.fobj(Atom_Pop(i,:));
    end    
    [Max_Fitness,Index]=min(Fitness);
    if Max_Fitness<Functon_Best(Iteration)
        Functon_Best(Iteration)=Max_Fitness;
        X_Best=Atom_Pop(Index,:);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Max_Fitness;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Index;
    else
        r=fix(rand*Atom_Num)+1;
        Atom_Pop(r,:)=X_Best;
    end
    
    Atom_Acc=Acceleration(Atom_Pop,Fitness,Iteration,Max_Iteration,Dim,Atom_Num,X_Best,alpha,beta);
    IterCurve(Iteration)=Max_Fitness;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iteration;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function Acc=Acceleration(Atom_Pop,Fitness,Iteration,Max_Iteration,Dim,Atom_Num,X_Best,alpha,beta)
%Calculate mass 
M=exp(-(Fitness-max(Fitness))./(max(Fitness)-min(Fitness)));
M=M./sum(M);  
G=exp(-20*Iteration/Max_Iteration); 
Kbest=Atom_Num-(Atom_Num-2)*(Iteration/Max_Iteration)^0.5;
Kbest=floor(Kbest)+1;
[Des_M Index_M]=sort(M,'descend');
 
for i=1:Atom_Num       
    E(i,:)=zeros(1,Dim);   
    MK(1,:)=sum(Atom_Pop(Index_M(1:Kbest),:),1)/Kbest;
    Distance=norm(Atom_Pop(i,:)-MK(1,:));   
    for k=1:Kbest
        j=Index_M(k);       
        %Calculate LJ-potential
        Potential=LJPotential(Atom_Pop(i,:),Atom_Pop(j,:),Iteration,Max_Iteration,Distance);                   
        E(i,:)=E(i,:)+rand(1,Dim)*Potential.*((Atom_Pop(j,:)-Atom_Pop(i,:))/(norm(Atom_Pop(i,:)-Atom_Pop(j,:))+eps));             
    end
    E(i,:)=alpha*E(i,:)+beta*(X_Best-Atom_Pop(i,:));
    %Calculate acceleration
    a(i,:)=E(i,:)./M(i); 
end
Acc=a.*G;
end

function Potential=LJPotential(Atom1,Atom2,Iteration,Max_Iteration,s)
r=norm(Atom1-Atom2);  
c=(1-(Iteration-1)/Max_Iteration).^3;  
rsmin=1.1+0.1*sin(Iteration/Max_Iteration*pi/2);
rsmax=1.24;

if r/s<rsmin
    rs=rsmin;
else
    if  r/s>rsmax
        rs=rsmax;  
    else
        rs=r/s;
    end
end           
Potential=c*(12*(-rs)^(-13)-6*(-rs)^(-7)); 
end