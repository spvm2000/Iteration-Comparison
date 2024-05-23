function Data_o = HBA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
tmax = Data_i.maxIter;      
N = Data_i.pop;            
dim=Data_i.dim;
X = Data_i.X;               
fitness = Data_i.F_value;         
lb = Data_i.lb;
ub = Data_i.ub;
IterCurve=zeros(1,Data_i.maxIter);
beta    = 6;     % the ability of HB to get the food  Eq.(4)
C       = 2;     %constant in Eq. (3)
vec_flag=[1,-1];
[GYbest, gbest] = min(fitness);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=GYbest;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=gbest;
Xprey = X(gbest,:);

for t = 1:tmax
    alpha=C*exp(-t/tmax);   %density factor in Eq. (3)
    I=Intensity(N,Xprey,X); %intensity in Eq. (2)
    for i=1:N
        r =rand();
        F=vec_flag(floor(2*rand()+1));
        for j=1:1:dim
            di=((Xprey(j)-X(i,j)));
            if r<.5
                r3=rand;                r4=rand;                r5=rand;
                
                Xnew(i,j)=Xprey(j) +F*beta*I(i)* Xprey(j)+F*r3*alpha*(di)*abs(cos(2*pi*r4)*(1-cos(2*pi*r5)));
            else
                r7=rand;
                Xnew(i,j)=Xprey(j)+F*r7*alpha*di;
            end
        end
        FU=Xnew(i,:)>ub;FL=Xnew(i,:)<lb;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        
        tempFitness = Data_i.fobj(Xnew(i,:));
        if tempFitness<fitness(i)
            fitness(i)=tempFitness;
            X(i,:)= Xnew(i,:);
        end
    end
    FU=X>ub;FL=X<lb;X=(X.*(~(FU+FL)))+ub.*FU+lb.*FL;
    [Ybest,index] = min(fitness);
    IterCurve(t)=Ybest;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Ybest;
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    if Ybest<GYbest
%         GYbest=Ybest;
        Xprey = X(index,:);
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);          
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end


function I=Intensity(N,Xprey,X)
    for i=1:N-1
        di(i) =( norm((X(i,:)-Xprey+eps))).^2;
        S(i)=( norm((X(i,:)-X(i+1,:)+eps))).^2;
    end
    di(N)=( norm((X(N,:)-Xprey+eps))).^2;
    S(N)=( norm((X(N,:)-X(1,:)+eps))).^2;
    for i=1:N
        r2=rand;
        I(i)=r2*S(i)/(4*pi*di(i));
    end
end
