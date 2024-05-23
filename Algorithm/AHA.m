function Data_o = AHA(Data_i,Data_o)
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

BestF=min(PopFit);
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(PopFit);

% Initialize visit table
Curve=zeros(1,Data_i.maxIter);
VisitTable=zeros(nPop) ;
VisitTable(logical(eye(nPop)))=NaN;    
    
    for It=1:MaxIt
        DirectVector=zeros(nPop,Dim);% Direction vector/matrix

        for i=1:nPop
            r=rand;
            if r<1/3     % Diagonal flight
                RandDim=randperm(Dim);
                if Dim>=3
                    RandNum=ceil(rand*(Dim-2)+1);
                else
                    RandNum=ceil(rand*(Dim-1)+1);
                end
                DirectVector(i,RandDim(1:RandNum))=1;
            else
                if r>2/3  % Omnidirectional flight
                    DirectVector(i,:)=1;
                else  % Axial flight
                    RandNum=ceil(rand*Dim);
                    DirectVector(i,RandNum)=1;
                end
            end

            if rand<0.5   % Guided foraging
                [MaxUnvisitedTime,TargetFoodIndex]=max(VisitTable(i,:));
                MUT_Index=find(VisitTable(i,:)==MaxUnvisitedTime);

                if length(MUT_Index)>1
                    [~,Ind]= min(PopFit(MUT_Index));
                    TargetFoodIndex=MUT_Index(Ind);
                end

                newPopPos=PopPos(TargetFoodIndex,:)+randn*DirectVector(i,:).*...
                    (PopPos(i,:)-PopPos(TargetFoodIndex,:));
                newPopPos=SpaceBound(newPopPos,Up,Low);
                newPopFit=Data_i.fobj(newPopPos);
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(i,TargetFoodIndex)=0;
                end
            else    % Territorial foraging
                newPopPos= PopPos(i,:)+randn*DirectVector(i,:).*PopPos(i,:);
                newPopPos=SpaceBound(newPopPos,Up,Low);
                newPopFit=Data_i.fobj(newPopPos);
                if newPopFit<PopFit(i)
                    PopFit(i)=newPopFit;
                    PopPos(i,:)=newPopPos;
                    VisitTable(i,:)=VisitTable(i,:)+1;
                    VisitTable(:,i)=max(VisitTable,[],2)+1;
                    VisitTable(i,i)=NaN;
                else
                    VisitTable(i,:)=VisitTable(i,:)+1;
                end
            end
        end

        if mod(It,2*nPop)==0 % Migration foraging
            [~, MigrationIndex]=max(PopFit);
            PopPos(MigrationIndex,:) =rand(1,Dim).*(Up-Low)+Low;
            PopFit(MigrationIndex)=Data_i.fobj(PopPos(MigrationIndex,:));
            VisitTable(MigrationIndex,:)=VisitTable(MigrationIndex,:)+1;
            VisitTable(:,MigrationIndex)=max(VisitTable,[],2)+1;
            VisitTable(MigrationIndex,MigrationIndex)=NaN;            
        end

        for i=1:nPop
            if PopFit(i)<BestF
                BestF=PopFit(i);
                Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
            end
        end
        Curve(It) = BestF;
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = BestF;
    end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=It;                                
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Curve;
end

function  X=SpaceBound(X,Up,Low)
    Dim=length(X);
    S=(X>Up)+(X<Low);    
    X=(rand(1,Dim).*(Up-Low)+Low).*S+X.*(~S);
end    