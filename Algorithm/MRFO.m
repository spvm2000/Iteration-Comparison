function Data_o = MRFO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                  
%% Problem Definition
MaxIt = Data_i.maxIter;             
nPop = Data_i.pop;                 
Dim = Data_i.dim;              
PopPos = Data_i.X;                   
PopFit = Data_i.F_value;     
Up=Data_i.ub;
Low=Data_i.lb;
IterCurve=zeros(1,Data_i.maxIter);
[BestF,index]=min(PopFit);
BestX=PopPos(index,:);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(PopFit);

for It=1:MaxIt  
     Coef=It/MaxIt; 
     
       if rand<0.5
          r1=rand;                         
          Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
          if  Coef>rand                                                      
              newPopPos(1,:)=BestX+rand(1,Dim).*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); %Equation (4)
          else
              IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
              newPopPos(1,:)=IndivRand+rand(1,Dim).*(IndivRand-PopPos(1,:))+Beta*(IndivRand-PopPos(1,:)); %Equation (7)         
          end              
       else 
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(1,:)=PopPos(1,:)+rand(1,Dim).*(BestX-PopPos(1,:))+Alpha.*(BestX-PopPos(1,:)); %Equation (1)
       end
     
    for i=2:nPop
        if rand<0.5
           r1=rand;                         
           Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
             if  Coef>rand                                                      
                 newPopPos(i,:)=BestX+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(BestX-PopPos(i,:)); %Equation (4)
             else
                 IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
                 newPopPos(i,:)=IndivRand+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(IndivRand-PopPos(i,:));  %Equation (7)       
             end              
        else
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(i,:)=PopPos(i,:)+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Alpha.*(BestX-PopPos(i,:)); %Equation (1)
       end         
    end
         
       for i=1:nPop        
           newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
           newPopFit(i)=Data_i.fobj(newPopPos(i,:));    
          if newPopFit(i)<PopFit(i)
             PopFit(i)=newPopFit(i);
             PopPos(i,:)=newPopPos(i,:);
          end
       end
     
            S=2;

        for i=1:nPop           
            newPopPos(i,:)=PopPos(i,:)+S*(rand*BestX-rand*PopPos(i,:)); %Equation (8)
        end
     
     for i=1:nPop        
         newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
         newPopFit(i)=Data_i.fobj(newPopPos(i,:));    
         if newPopFit(i)<PopFit(i)
            PopFit(i)=newPopFit(i);
            PopPos(i,:)=newPopPos(i,:);
         end
     end
     
     for i=1:nPop
        if PopFit(i)<BestF
           BestF=PopFit(i);
           BestX=PopPos(i,:);
           Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
     end
    IterCurve(It)=BestF;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestF;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=It;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function  X=SpaceBound(X,Up,Low)
    Dim=length(X);
    S=(X>Up)+(X<Low);    
    X=(rand(1,Dim).*(Up-Low)+Low).*S+X.*(~S);
end