function Data_o = SDO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
X = Data_i.X;                        
Cost=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
MaxIt=Data_i.maxIter;        MarketSize=Data_i.pop;  
Low=Data_i.lb;  Up=Data_i.ub;  Dim=Data_i.dim;      BenFunctions=Data_i.fobj;
OneMarket.CommPrice=[];         % Commodity price
OneMarket.CommPriceFit=[];      % Fitness of commodity price

OneMarket.CommQuantity=[];      % Commodity quantity
OneMarket.CommQuantityFit=[];   % Fitness of commodity quantity

Market=repmat(OneMarket,MarketSize,1);
Matr=[1 Dim];
BestF=inf;
BestX=[];

for i=1:MarketSize
     Market(i).CommPrice=X(i,:);
     Market(i).CommPriceFit=BenFunctions(Market(i).CommPrice);
     Market(i).CommQuantity=rand(1,Dim).*(Up-Low)+Low; 
     Market(i).CommQuantityFit=BenFunctions(Market(i).CommQuantity);    
     if Market(i).CommQuantityFit<=Market(i).CommPriceFit
        Market(i).CommPriceFit=Market(i).CommQuantityFit;
        Market(i).CommPrice=Market(i).CommQuantity;
     end    
end

for i=1:MarketSize
   if Market(i).CommPriceFit<=BestF
      BestX=Market(i).CommPrice;
      BestF=Market(i).CommPriceFit;
   end
end

for Iter=1:MaxIt   
    a=2*(MaxIt-Iter+1)/MaxIt;    
    F=zeros(MarketSize,1); 
    
    MeanQuantityFit = mean([Market.CommQuantityFit]);     
    for i=1:MarketSize       
        F(i) =(abs(Market(i).CommQuantityFit-MeanQuantityFit)+10^(-15)); % Equation (7)
    end    
    FQ=F/sum(F);  % Equation (8)
    
    MeanPriceFit=mean([Market.CommPriceFit]);  
    for i=1:MarketSize       
        F(i) =(abs(Market(i).CommPriceFit-MeanPriceFit)+10^(-15));%Equation (10)
    end 
    FP=F/sum(F);       % Equation (11)
    MeanPrice=(mean(reshape([Market.CommPrice],Dim,MarketSize),2))';  
          
    for i=1:MarketSize 
        Ind=round(rand)+1;
        k=find(rand<=cumsum(FQ),1,'first'); % Equation (9)
        CommQuantityEqu=Market(k).CommQuantity;
        
        Alpha=a*sin((2*pi)*rand(1,Matr(Ind))); % Equation (16)
        Beta=2*cos((2*pi)*rand(1,Matr(Ind))); % Equation (17)
            
         if rand>0.5               
            CommPriceEqu=rand*MeanPrice;% Equation (12)
          else
            k=find(rand<=cumsum(FP),1,'first');   
            CommPriceEqu=Market(k).CommPrice; % Equation (12)        
         end
          
          % Supply function (supply relation of producers)
         NewCommQuantity=CommQuantityEqu+Alpha.*(Market(i).CommPrice-CommPriceEqu); %Equation (13)       
         NewCommQuantity=BoundaryCheck(NewCommQuantity,Up,Low,Data_i.dim); 
         NewCommQuantityFit= BenFunctions(NewCommQuantity);   
         
         if NewCommQuantityFit<=Market(i).CommQuantityFit
            Market(i).CommQuantityFit=NewCommQuantityFit;
            Market(i).CommQuantity=NewCommQuantity;
         end
         
         % Demand function (demand relation of consumers)
         NewCommPrice=CommPriceEqu-Beta.*(NewCommQuantity-CommQuantityEqu );% Equation (14)
         NewCommPrice=BoundaryCheck(NewCommPrice,Up,Low,Data_i.dim);
         NewCommPriceFit= BenFunctions(NewCommPrice);    
         
         if NewCommPriceFit<=Market(i).CommPriceFit
            Market(i).CommPriceFit=NewCommPriceFit;
            Market(i).CommPrice=NewCommPrice;
         end       
    end    
   
    % Replacement
    for i=1:MarketSize        
        if Market(i).CommQuantityFit<=Market(i).CommPriceFit
           Market(i).CommPriceFit=Market(i).CommQuantityFit;
           Market(i).CommPrice=Market(i).CommQuantity;
        end        
    end
     
    for i=1:MarketSize        
        if Market(i).CommPriceFit<=BestF
           BestX=Market(i).CommPrice;
           BestF=Market(i).CommPriceFit;
           if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>BestF
                Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestF;
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
           end
        end            
    end       
  
    IterCurve(Iter)=BestF;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=Iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end