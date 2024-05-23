function Data_o = WCA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                              
X = Data_i.X;                        
fitness=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
LB=Data_i.lb;       UB=Data_i.ub;       
objective_function=Data_i.fobj;   max_it=Data_i.maxIter;
nvars=Data_i.dim;   
Npop=Data_i.pop;
dmax=1e-16;
Nsr=4;
N_stream=Npop-Nsr;

ind.position=[];
ind.cost=[];

pop=repmat(ind,Npop,1);
for i=1:Npop
    pop(i).position=X(i,:);
    pop(i).cost=fitness(i);
end
[~, index]=sort([pop.cost]);
sea=pop(index(1));

river=repmat(ind,Nsr-1,1);
for i=1:Nsr-1
    river(i)=pop(index(1+i));
end
%------------ Forming Streams----------------------------------------------
stream=repmat(ind,N_stream,1);
for i=1:N_stream
    stream(i)=pop(index(Nsr+i));
end
%--------- Designate streams to rivers and sea ------------------------
cs=[sea.cost;[river.cost]';stream(1).cost];
f=0;
if length(unique(cs))~=1
    CN=cs-max(cs);
else
    CN=cs;
    f=1;
end
NS=round(abs(CN/sum(CN))*N_stream);
if f~=1
    NS(end)=[];
end
NS=sort(NS,'descend');
i=Nsr;
while sum(NS)>N_stream
    if NS(i)>1
        NS(i)=NS(i)-1;
    else
        i=i-1;
    end
end

i=1;
while sum(NS)<N_stream
    NS(i)=NS(i)+1;
end

if find(NS==0)
    index=find(NS==0);
    for i=1:size(index,1)
        while NS(index(i))==0
            NS(index(i))=NS(index(i))+round(NS(i)/6);
            NS(i)=NS(i)-round(NS(i)/6);
        end
    end
end

NS=sort(NS,'descend');
NB=NS(2:end);
for i=1:max_it
    %---------- Moving stream to sea---------------------------------------
    for j=1:NS(1)
        stream(j).position=stream(j).position+2.*rand(1).*(sea.position-stream(j).position);
        
        stream(j).position=min(stream(j).position,UB);
        stream(j).position=max(stream(j).position,LB);
        
        stream(j).cost=objective_function(stream(j).position);
        
        if stream(j).cost<sea.cost
            new_sea=stream(j);
            stream(j)=sea;
            sea=new_sea;
        end
    end
    %---------- Moving Streams to rivers-----------------------------------
    for k=1:Nsr-1
        for j=1:NB(k)
            stream(j+sum(NS(1:k))).position=stream(j+sum(NS(1:k))).position+2.*rand(1,nvars).*(river(k).position-stream(j+sum(NS(1:k))).position);
            
            stream(j+sum(NS(1:k))).position=min(stream(j+sum(NS(1:k))).position,UB);
            stream(j+sum(NS(1:k))).position=max(stream(j+sum(NS(1:k))).position,LB);
            
            stream(j+sum(NS(1:k))).cost=objective_function(stream(j+sum(NS(1:k))).position);
            
            if stream(j+sum(NS(1:k))).cost<river(k).cost
                new_river=stream(j+sum(NS(1:k)));
                stream(j+sum(NS(1:k)))=river(k);
                river(k)=new_river;
  
                if river(k).cost<sea.cost
                    new_sea=river(k);
                    river(k)=sea;
                    sea=new_sea;
                end
            end
        end
    end
    %---------- Moving rivers to Sea --------------------------------------
    for j=1:Nsr-1
        river(j).position=river(j).position+2.*rand(1,nvars).*(sea.position-river(j).position);
        
        river(j).position=min(river(j).position,UB);
        river(j).position=max(river(j).position,LB);
        
        river(j).cost=objective_function(river(j).position);
        
        if river(j).cost<sea.cost
            new_sea=river(j);
            river(j)=sea;
            sea=new_sea;
        end
    end
    %-------------- Evaporation condition and raining process--------------
    % Check the evaporation condition for rivers and sea
    for k=1:Nsr-1
        if ((norm(river(k).position-sea.position)<dmax) || rand<0.1)
            for j=1:NB(k)
                stream(j+sum(NS(1:k))).position=LB+rand(1,nvars).*(UB-LB);
            end
        end
    end
    % Check the evaporation condition for streams and sea
    for j=1:NS(1)
        if ((norm(stream(j).position-sea.position)<dmax))
            stream(j).position=LB+rand(1,nvars).*(UB-LB);
        end
    end
    %----------------------------------------------------------------------
    dmax=dmax-(dmax/max_it);
    IterCurve(i)=sea.cost;
end
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);       
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=i;                        
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end