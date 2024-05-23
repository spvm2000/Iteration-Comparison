function Data_o = MPA(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
%% Problem Definition
Max_iter = Data_i.maxIter;    
SearchAgents_no = Data_i.pop;               
dim = Data_i.dim;             
Prey = Data_i.X;               
fit = Data_i.F_value;         

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
[Top_predator_fit,index]=min(fit); 
Top_predator_pos=Prey(index,:);

Convergence_curve=zeros(1,Max_iter);
stepsize=zeros(SearchAgents_no,dim);
fitness=inf(SearchAgents_no,1);
  
Xmin=repmat(ones(1,dim).*Data_i.lb,SearchAgents_no,1);
Xmax=repmat(ones(1,dim).*Data_i.ub,SearchAgents_no,1);
         

Iter=0;
FADs=0.2;
P=0.5;

while Iter<Max_iter    
     
     %------------------- Detecting top predator -----------------    
     for i=1:size(Prey,1)  
               
        Prey=BoundaryCheck(Prey,Data_i.ub,Data_i.lb,Data_i.dim);                    
        fitness(i,1)=Data_i.fobj(Prey(i,:));
                         
         if fitness(i,1)<Top_predator_fit 
           Top_predator_fit=fitness(i,1); 
           Top_predator_pos=Prey(i,:);
           Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
         end          
     end

     %------------------- Marine Memory saving ------------------- 
    
 if Iter==0
   fit_old=fitness;    Prey_old=Prey;
 end
     
  Inx=(fit_old<fitness);
  Indx=repmat(Inx,1,dim);
  Prey=Indx.*Prey_old+~Indx.*Prey;
  fitness=Inx.*fit_old+~Inx.*fitness;
        
  fit_old=fitness;    Prey_old=Prey;

     %------------------------------------------------------------   
     
 Elite=repmat(Top_predator_pos,SearchAgents_no,1);  %(Eq. 10) 
 CF=(1-Iter/Max_iter)^(2*Iter/Max_iter);
                             
 RL=0.05*levy(SearchAgents_no,dim,1.5);   %Levy random number vector
 RB=randn(SearchAgents_no,dim);          %Brownian random number vector
           
  for i=1:size(Prey,1)
     for j=1:size(Prey,2)        
       R=rand();
          %------------------ Phase 1 (Eq.12) ------------------- 
       if Iter<Max_iter/3 
          stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));                    
          Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j); 
             
          %--------------- Phase 2 (Eqs. 13 & 14)----------------
       elseif Iter>Max_iter/3 && Iter<2*Max_iter/3 
          
         if i>size(Prey,1)/2
            stepsize(i,j)=RB(i,j)*(RB(i,j)*Elite(i,j)-Prey(i,j));
            Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j); 
         else
            stepsize(i,j)=RL(i,j)*(Elite(i,j)-RL(i,j)*Prey(i,j));                     
            Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);  
         end  
         
         %----------------- Phase 3 (Eq. 15)-------------------
       else 
           
           stepsize(i,j)=RL(i,j)*(RL(i,j)*Elite(i,j)-Prey(i,j)); 
           Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);  
    
       end  
      end                                         
  end    
        
     %------------------ Detecting top predator ------------------        
  for i=1:size(Prey,1)  
        
      Prey=BoundaryCheck(Prey,Data_i.ub,Data_i.lb,Data_i.dim);
      fitness(i,1)=Data_i.fobj(Prey(i,:));
        
      if fitness(i,1)<Top_predator_fit 
         Top_predator_fit=fitness(i,1);
         Top_predator_pos=Prey(i,:);
         Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
      end     
  end
        
     %---------------------- Marine Memory saving ----------------
    
 if Iter==0
    fit_old=fitness;    Prey_old=Prey;
 end
     
    Inx=(fit_old<fitness);
    Indx=repmat(Inx,1,dim);
    Prey=Indx.*Prey_old+~Indx.*Prey;
    fitness=Inx.*fit_old+~Inx.*fitness;
        
    fit_old=fitness;    Prey_old=Prey;

     %---------- Eddy formation and FADs effect (Eq 16) ----------- 
                             
  if rand()<FADs
     U=rand(SearchAgents_no,dim)<FADs;                                                                                              
     Prey=Prey+CF*((Xmin+rand(SearchAgents_no,dim).*(Xmax-Xmin)).*U);
  else
     r=rand();  Rs=size(Prey,1);
     stepsize=(FADs*(1-r)+r)*(Prey(randperm(Rs),:)-Prey(randperm(Rs),:));
     Prey=Prey+stepsize;
  end
                                                        
  Iter=Iter+1;
  Convergence_curve(Iter) = Top_predator_fit;
  Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = Top_predator_fit;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end

function [z] = levy(n,m,beta)

    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator 
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator

    sigma_u = (num/den)^(1/beta);% Standard deviation

    u = random('Normal',0,sigma_u,n,m); 
    
    v = random('Normal',0,1,n,m);

    z =u./(abs(v).^(1/beta));

  
  end