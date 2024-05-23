function Data_o = ARO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                 
%% Problem Definition
MaxIt = Data_i.maxIter;      
nPop = Data_i.pop;            
Dim = Data_i.dim;              
PopPos = Data_i.X;               
PopFit = Data_i.F_value;         
Low = Data_i.lb;
Up = Data_i.ub;

[BestF,index]=min(PopFit);
Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=BestF;
Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
HisBestF=zeros(1,Data_i.maxIter);
for It=1:MaxIt
    Direct1=zeros(nPop,Dim);
    Direct2=zeros(nPop,Dim);
    theta=2*(1-It/MaxIt);
    for i=1:nPop
        L=(exp(1)-exp(((It-1)/MaxIt)^2))*(sin(2*pi*rand)); %Eq.(3)
        rd=ceil(rand*(Dim));
        Direct1(i,randperm(Dim,rd))=1;
        c=Direct1(i,:); %Eq.(4)
        R=L.*c; %Eq.(2)
        
        A=2*log(1/rand)*theta;%Eq.(15)

        if A>1

            K=[1:i-1 i+1:nPop];
            RandInd=K(randi([1 nPop-1]));
            newPopPos=PopPos(RandInd,:)+R.*( PopPos(i,:)-PopPos(RandInd,:))...
                +round(0.5*(0.05+rand))*randn; %Eq.(1)
        else

            Direct2(i,ceil(rand*Dim))=1;
            gr=Direct2(i,:); %Eq.(12)
            H=((MaxIt-It+1)/MaxIt)*randn; %Eq.(8)
            b=PopPos(i,:)+H*gr.*PopPos(i,:); %Eq.(13)
            newPopPos=PopPos(i,:)+ R.*(rand*b-PopPos(i,:)); %Eq.(11)

        end
        newPopPos=BoundaryCheck(newPopPos,Up,Low,Data_i.dim);
        newPopFit=Data_i.fobj(newPopPos);
        if newPopFit<PopFit(i)
            PopFit(i)=newPopFit;
            PopPos(i,:)=newPopPos;
        end

    end

    for i=1:nPop
        if PopFit(i)<BestF
            BestF=PopFit(i);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    HisBestF(It)=BestF;
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestF;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=It;                               
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=HisBestF;
end
