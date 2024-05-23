function Data_o = COOT(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
Max_iter = Data_i.maxIter;    
N = Data_i.pop;               
dim = Data_i.dim;             
pos = Data_i.X;               
fit = Data_i.F_value;         

NLeader=ceil(0.1*N);  
Ncoot=N-NLeader;
Convergence_curve = zeros(1,Max_iter);
% gBest=zeros(1,dim);
% gBestScore=inf;
%Initialize the positions of Coots
CootPos=pos(1:Ncoot,:);                                %CootPos=rand(Ncoot,dim).*(ub-lb)+lb;
CootFitness=fit(1,1:Ncoot);                              %CootFitness=zeros(1,Ncoot);
%Initialize the locations of Leaders
LeaderPos=pos(Ncoot+1:N,:);                            %LeaderPos=rand(NLeader,dim).*(ub-lb)+lb;
LeaderFit=fit(1,Ncoot+1:N);                              %LeaderFit=zeros(1,NLeader);

% for i=1:size(CootPos,1)
%     CootFitness(1,i)=Data_i.fobj(CootPos(i,:));
%       if(gBestScore>CootFitness(1,i))
%             gBestScore=CootFitness(1,i);
%             gBest=CootPos(i,:);
%       end  
% end
% for i=1:size(LeaderPos,1)  
%     LeaderFit(1,i)=Data_i.fobj(LeaderPos(i,:));
%       if(gBestScore>LeaderFit(1,i))
%             gBestScore=LeaderFit(1,i);
%             gBest=LeaderPos(i,:);
%       end   
% end

[gBestScore,gBestPos]=min(fit);
gBest=pos(gBestPos,:);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
Convergence_curve(1)=gBestScore;

l=2; % Loop counter
while l<Max_iter+1
 B=2-l*(1/Max_iter);
 A=1-l*(1/Max_iter);

    for i=1:size(CootPos,1) 
         if rand<0.5
            R=-1+2*rand;
%             rand_num(l-1,1) = R;
            R1=rand();
%             rand_num(l-1,2) = R1;
         else  
            R=-1+2*rand(1,dim);
%             rand_num(l-1,1) = R;
            R1=rand(1,dim);
%             rand_num(l-1,2) = R1;
         end
        k=1+mod(i,NLeader);  %%除余函数
        if rand<0.5 %|| i==1
            CootPos(i,:)=2*R1.*cos(2*pi*R).*(LeaderPos(k,:)-CootPos(i,:))+LeaderPos(k,:); 
             % Check boundries
            CootPos=BoundaryCheck(CootPos,Data_i.ub,Data_i.lb,Data_i.dim);
        else
            if rand<0.5 && i~=1 %i>2*size(CootPos,1)/3%
                 CootPos(i,:)=(CootPos(i,:)+CootPos(i-1,:))/2;
            else
                Q=rand(1,dim).*(Data_i.ub-Data_i.lb)+Data_i.lb;
                  %R1=0.2+ 0.6*rand;                
                  CootPos(i,:)=CootPos(i,:)+A*R1.*(Q-CootPos(i,:));               
            end
            CootPos=BoundaryCheck(CootPos,Data_i.ub,Data_i.lb,Data_i.dim);  
        end
    end
    % fitness of location of Coots
   for i=1:size(CootPos,1)
    CootFitness(1,i)=Data_i.fobj(CootPos(i,:));
    k=1+mod(i,NLeader); 
    % Update the location of coot
    if CootFitness(1,i)<LeaderFit(1,k)
         Temp=LeaderPos(k,:);
         TemFit= LeaderFit(1,k);
       LeaderFit(1,k)=CootFitness(1,i);
       LeaderPos(k,:)=CootPos(i,:);
         CootFitness(1,i)=TemFit;
         CootPos(i,:)=Temp;      
    end
   end
    % fitness of location of Leaders
    for i=1:size(LeaderPos,1)
         if rand<0.5
            R2=-1+2*rand;
%             rand_num(l-1,3) = R2;
            R3=rand();
%             rand_num(l-1,4) = R3;
         else  
            R2=-1+2*rand(1,dim);
%             rand_num(l-1,3) = R2;
            R3=rand(1,dim);
%             rand_num(l-1,4) = R3;
         end
        if rand<0.5             
            Temp=B*R3.*cos(2*pi*R2).*(gBest-LeaderPos(i,:))+gBest;           
        else           
            Temp=B*R3.*cos(2*pi*R2).*(gBest-LeaderPos(i,:))-gBest;           
        end
        % Check boundries
        Temp=BoundaryCheck(Temp,Data_i.ub,Data_i.lb,Data_i.dim);
        TempFit=Data_i.fobj(Temp);
        % Update the location of Leader
         if(gBestScore>TempFit)
             LeaderFit(1,i)=gBestScore;
             LeaderPos(i,:)=gBest;
             Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
             gBestScore=TempFit;
             gBest=Temp;
         end
    end
    Convergence_curve(l)=gBestScore;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = gBestScore;
    l = l + 1;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;

