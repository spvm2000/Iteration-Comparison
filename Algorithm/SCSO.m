function Data_o = SCSO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
Max_iter = Data_i.maxIter;     
SearchAgents_no = Data_i.pop;           
Positions = Data_i.X;               
fitness = Data_i.F_value;         
lb = Data_i.lb;
ub = Data_i.ub;

[Best_Score,index]=min(fitness);
BestFit=Positions(index,:);
Convergence_curve=zeros(1,Max_iter);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Best_Score;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
t=0;
p=[1:360];
while t<Max_iter

    S=2;                                    %%% S is maximum Sensitivity range 
    rg=S-((S)*t/(Max_iter));                %%%% guides R
   for i=1:size(Positions,1)
        r=rand*rg;
        R=((2*rg)*rand)-rg;                 %%%%   controls to transtion phases  
        for j=1:size(Positions,2)
        teta=RouletteWheelSelection(p);
           if((-1<=R)&&(R<=1))              %%%% R value is between -1 and 1
                Rand_position=abs(rand*BestFit(j)-Positions(i,j));
                Positions(i,j)=BestFit(j)-r*Rand_position*cos(teta);
           else                 
                cp=floor(SearchAgents_no*rand()+1);
                CandidatePosition =Positions(cp,:);
                Positions(i,j)=r*(CandidatePosition(j)-rand*Positions(i,j));
            end
        end
    end
    t=t+1;
    for i=1:size(Positions,1)
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness=Data_i.fobj(Positions(i,:));
        if fitness<Best_Score
            Best_Score=fitness;
            BestFit=Positions(i,:);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    Convergence_curve(t)=Best_Score;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Best_Score;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end

function j=RouletteWheelSelection(P) 
    r=rand; 
    s=sum(P);
    P=P./s;
    C=cumsum(P); 
    j=find(r<=C,1,'first'); 
end