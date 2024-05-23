function Data_o = BES(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
MaxIt = Data_i.maxIter;      
nPop = Data_i.pop;            
dim = Data_i.dim;              
low = Data_i.lb;
high = Data_i.ub;
pop.pos = Data_i.X;               
pop.cost = Data_i.F_value;        
Curve=zeros(1,Data_i.maxIter);
% Initialize Best Solution
[BestSol.cost,index]=min(pop.cost);
BestSol.pos=(pop.pos(index,:));
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(pop.cost);
    for t=1:MaxIt
        %%               1- select_space 
        [pop ,BestSol,s1]=select_space(Data_i.fobj,pop,nPop,BestSol,low,high,dim);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=s1;
        %%                2- search in space
        [pop ,BestSol ,s2]=search_space(Data_i.fobj,pop,BestSol,nPop,low,high);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=s2;
        %%                3- swoop
        [pop ,BestSol ,s3]=swoop(Data_i.fobj,pop,BestSol,nPop,low,high);
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=s3;

        Curve(t)=BestSol.cost;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestSol.cost;
        
    end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Curve;
end

function [pop ,BestSol ,s1]=select_space(fobj,pop,npop,BestSol,low,high,dim)
    Mean=mean(pop.pos);
    % Empty Structure for Individuals
    empty_individual.pos = [];
    empty_individual.cost = [];
    lm= 2;
    s1=inf;
    for i=1:npop
        newsol=empty_individual;
        newsol.pos= BestSol.pos+ lm*rand(1,dim).*(Mean - pop.pos(i,:));
        newsol.pos = max(newsol.pos, low);
        newsol.pos = min(newsol.pos, high);
        newsol.cost=fobj(newsol.pos);
        if newsol.cost<pop.cost(i)
           pop.pos(i,:) = newsol.pos;
           pop.cost(i)= newsol.cost;
           
             if pop.cost(i) < BestSol.cost
                  BestSol.pos= pop.pos(i,:);
                  BestSol.cost=pop.cost(i);
                  s1=i;
             end
        end
    end
end

function [pop ,best ,s2]=search_space(fobj,pop,best,npop,low,high)
    Mean=mean(pop.pos);
    a=10;
    R=1.5;
    s2=inf;
    % Empty Structure for Individuals
    empty_individual.pos = [];
    empty_individual.cost = [];
    for i=1:npop-1
        A=randperm(npop);
        pop.pos=pop.pos(A,:);
        pop.cost=pop.cost(A);
        [x ,y]=polr(a,R,npop);
        newsol=empty_individual;
        Step = pop.pos(i,:) - pop.pos(i+1,:);
        Step1=pop.pos(i,:)-Mean;
        newsol.pos = pop.pos(i,:) +y(i)*Step+x(i)*Step1;
        newsol.pos = max(newsol.pos, low);
        newsol.pos = min(newsol.pos, high);
        newsol.cost=fobj(newsol.pos);
        if newsol.cost<pop.cost(i)
           pop.pos(i,:) = newsol.pos;
           pop.cost(i)= newsol.cost;
           if pop.cost(i) < best.cost
               best.pos= pop.pos(i,:);
               best.cost=pop.cost(i);
               s2=i;
           end
        end
    end
end

function [pop ,best ,s3]=swoop(fobj,pop,best,npop,low,high)
    Mean=mean(pop.pos);
    a=10;
    s3=inf;
    % Empty Structure for Individuals
    empty_individual.pos = [];
    empty_individual.cost = [];
    for i=1:npop
        A=randperm(npop);
        pop.pos=pop.pos(A,:);
        pop.cost=pop.cost(A);
        [x ,y]=swoo_p(a,npop);
        newsol=empty_individual;
        Step = pop.pos(i,:) - 2*Mean;
        Step1= pop.pos(i,:)-2*best.pos;
        newsol.pos = rand(1,length(Mean)).*best.pos+x(i)*Step+y(i)*Step1;
        newsol.pos = max(newsol.pos, low);
        newsol.pos = min(newsol.pos, high);
        newsol.cost=fobj(newsol.pos);
        if newsol.cost<pop.cost(i)
           pop.pos(i,:) = newsol.pos;
           pop.cost(i)= newsol.cost;
           if pop.cost(i) < best.cost
                best.pos= pop.pos(i,:);
                best.cost=pop.cost(i);
                s3=i;
           end
        end
    end
end

function [xR ,yR]=swoo_p(a,N)
    th = a*pi*exp(rand(N,1));
    r  =th; %R*rand(N,1);
    xR = r.*sinh(th);
    yR = r.*cosh(th);
    xR=xR/max(abs(xR));
    yR=yR/max(abs(yR));
end
 
function [xR ,yR]=polr(a,R,N)
    %// Set parameters
    th = a*pi*rand(N,1);
    r  =th+R*rand(N,1);
    xR = r.*sin(th);
    yR = r.*cos(th);
    xR=xR/max(abs(xR));
    yR=yR/max(abs(yR));
end