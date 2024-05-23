function Data_o = CHIO(Data_i,Data_o)
rand_num=[];                        
ti=clock;                               
swarm = Data_i.X;                        
ObjVal=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
PopSize=Data_i.pop;         Max_iter=Data_i.maxIter;
MaxAge = 100;       C0 = 1;         SpreadingRate = 0.05;   
runs = 1;
ObjVal = zeros(1,PopSize);      Age = zeros(1,PopSize);
BestResults = zeros(runs,1);
lb=Data_i.lb(1);       ub=Data_i.ub(1);       dim=Data_i.dim;     fobj=Data_i.fobj;

Fitness=calculateFitness(ObjVal);
Status=zeros(1,PopSize);
for i=1:C0
    Status(fix(rand*(PopSize))+1)=1;  
end
itr=1;
while itr<=Max_iter
  for i=1:PopSize
      NewSol=swarm(i,:);
      CountCornoa = 0;
      % find the set of confirmed solutions
      confirmed = randperm(size(find(Status==1),2));
      confirmed1 = find(Status==1);
      %find(Status==1);
      % find the set of normal solutions
      normal = randperm(size(find(Status==0),2));
      normal1 = find(Status==0);
      % find the set of recovered solutions
      recovered = find(ObjVal & Status==2);
      [cost,Index3]=min(recovered);
      for j=1: dim 
          r = rand();  % select a number within range 0 to 1.
         if ((r < SpreadingRate/3)&&(size(confirmed1,2)>0))
              % select one of the confirmed solutions
              z=round(1+(size(confirmed1,2)-1)*rand);   
              zc= confirmed1(z);
              % modify the curent value
              NewSol(j) = swarm(i,j)+(swarm(i,j)-swarm(zc,j))*(rand-0.5)*2;
              % manipulate range between lb and ub
              NewSol(j)= min(max(NewSol(j),lb),ub); 
              CountCornoa = CountCornoa + 1;                              
          
          elseif ((r < SpreadingRate/2) &&size(normal1,2)>0)
              % select one of the normal solutions
              z=round(1+(size(normal1,2)-1)*rand);
              zn= normal1(z);
              % modify the curent value
              NewSol(j) = swarm(i,j)+(swarm(i,j)-swarm(zn,j))*(rand-0.5)*2;
              % manipulate range between lb and ub
              NewSol(j)= min(max(NewSol(j),lb),ub); 
          
          elseif (r < SpreadingRate && size(recovered,2)>0)
              % modify the curent value
              NewSol(j) = swarm(i,j)+(swarm(i,j)-swarm(Index3,j))*(rand-0.5)*2;
              % manipulate range between lb and ub
              NewSol(j)= min(max(NewSol(j),lb),ub);  
          end
     end

      %evaluate new solution
      ObjValSol=fobj(NewSol);
      FitnessSol=calculateFitness(ObjValSol);
      
      % Update the curent solution  & Age of the current solution
      if (ObjVal(i)>ObjValSol) 
        swarm(i,:)=NewSol;
        Fitness(i)=FitnessSol;
        ObjVal(i)=ObjValSol;
      else
          if(Status(i)==1)
              Age(i) = Age(i) + 1;
          end
      end            
                   
      % change the solution from normal to confirmed
      if ((Fitness(i) < mean(Fitness))&& Status(i)==0 && CountCornoa>0)
          Status(i) = 1;
          Age(i)=1;
      end
      
      % change the solution from confirmed to recovered
      if ((Fitness(i) >= mean(Fitness))&& Status(i)==1)
          Status(i) = 2; 
          Age(i)=0;
      end
      
      % killed the current soluion and regenerated from scratch
      if(Age(i)>=MaxAge)
          NewSolConst = CHIO_initialization(1,dim,ub,lb);
          swarm(i,:) = NewSolConst(:);
          Status(i) = 0;
      end
  end
       
  IterCurve(itr)=min(ObjVal);
  if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>min(ObjVal)
      [Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(ObjVal);      
  end
  itr=itr+1;    
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=itr;                            
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function fFitness=calculateFitness(fObjV)
    fFitness=zeros(size(fObjV));
    ind=find(fObjV>=0);
    fFitness(ind)=1./(fObjV(ind)+1);
    ind=find(fObjV<0);
    fFitness(ind)=1+abs(fObjV(ind));
end

function Positions=CHIO_initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end