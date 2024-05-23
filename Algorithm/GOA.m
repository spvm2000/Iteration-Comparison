function Data_o = GOA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                            
GrassHopperPositions = Data_i.X;                        
GrassHopperFitness=Data_i.F_value;                 
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
N=Data_i.pop;     Max_iter=Data_i.maxIter;   
lb=Data_i.lb;   ub=Data_i.ub;   dim=Data_i.dim;     fobj=Data_i.fobj;

fitness_history=zeros(N,Max_iter);
position_history=zeros(N,Max_iter,dim);
Trajectories=zeros(N,Max_iter);

cMax=1;
cMin=0.00004;

for i=1:size(GrassHopperPositions,1)
    if flag == 1
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
    else
        GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
    end
    fitness_history(i,1)=GrassHopperFitness(1,i);
    position_history(i,1,:)=GrassHopperPositions(i,:);
    Trajectories(:,1)=GrassHopperPositions(:,1);
end

[sorted_fitness,sorted_indexes]=sort(GrassHopperFitness);
for newindex=1:Data_i.pop
    Sorted_grasshopper(newindex,:)=GrassHopperPositions(sorted_indexes(newindex),:);
end

TargetPosition=Sorted_grasshopper(1,:);
TargetFitness=sorted_fitness(1);
IterCurve(1)=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
l=2;
while l<Max_iter+1
    c=cMax-l*((cMax-cMin)/Max_iter); % Eq. (2.8) in the paper
    for i=1:size(GrassHopperPositions,1)
        temp= GrassHopperPositions';
        for k=1:2:dim
            S_i=zeros(2,1);
            for j=1:N
                if i~=j
                    if k==dim
                        k=dim-1;
                    end    
                    Dist=distance(temp(k:k+1,j), temp(k:k+1,i)); % Calculate the distance between two grasshoppers
                    
                    r_ij_vec=(temp(k:k+1,j)-temp(k:k+1,i))/(Dist+eps); % xj-xi/dij in Eq. (2.7)
                    xj_xi=2+rem(Dist,2); % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij=((ub(k) - lb(k))*c/2)*S_func(xj_xi).*r_ij_vec; % The first part inside the big bracket in Eq. (2.7)
                    S_i=S_i+s_ij;
                end
            end
            S_i_total(k:k+1, :) = S_i;
        end
        X_new = c * S_i_total'+ (TargetPosition);    % Eq. (2.7) in the paper      
        GrassHopperPositions_temp(i,:)=X_new'; 
    end
    % GrassHopperPositions
    GrassHopperPositions=GrassHopperPositions_temp;
    
    for i=1:size(GrassHopperPositions,1)
        % Relocate grasshoppers that go outside the search space 
        Tp=GrassHopperPositions(i,:)>ub;
        Tm=GrassHopperPositions(i,:)<lb;
        GrassHopperPositions(i,:)=(GrassHopperPositions(i,:).*(~(Tp+Tm)))+ub.*Tp+lb.*Tm;
        
        % Calculating the objective values for all grasshoppers
        if flag == 1
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,1:end-1));
        else
            GrassHopperFitness(1,i)=fobj(GrassHopperPositions(i,:));
        end
        fitness_history(i,l)=GrassHopperFitness(1,i);
        position_history(i,l,:)=GrassHopperPositions(i,:);
        
        Trajectories(:,l)=GrassHopperPositions(:,1);
        
        % Update the target
        if GrassHopperFitness(1,i)<TargetFitness
            TargetPosition=GrassHopperPositions(i,:);
            TargetFitness=GrassHopperFitness(1,i);
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=TargetFitness;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
    end
    IterCurve(l)=TargetFitness;
    l = l + 1;
end
%% end 扫尾工作
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);             
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function d = distance(a,b)
    d=sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
end

function o=S_func(r)
    f=0.5;
    l=1.5;
    o=f*exp(-r/l)-exp(-r);  % Eq. (2.3) in the paper
end