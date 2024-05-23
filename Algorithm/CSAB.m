function Data_o = CSAB(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
current_X = Data_i.X;                       
current_X_fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
pop=Data_i.pop;     gbestNum=3;     Max_iteration=Data_i.maxIter;
lb=Data_i.lb;   ub=Data_i.ub;   dim=Data_i.dim;     fobj=Data_i.fobj;
global_Best_position=zeros(gbestNum,dim);% global_Best_position
global_Best_fitness=inf(1,gbestNum);
Destination_position=zeros(1,dim);%Destination_position
Destination_fitness=inf;
Person_Best_position=current_X;%Person_Best_position
Person_Best_fitness=current_X_fitness;

[~,index_sorted]=sort(current_X_fitness);
for i=1:gbestNum
    global_Best_position(i,:)=current_X(index_sorted(i),:);
    global_Best_fitness(i)=current_X_fitness(1,index_sorted(i));
end
iter=1;
iter=1;
while iter<=Max_iteration
    ave_Pbest=mean(Person_Best_position,1);
    ave_gbest=mean(global_Best_position,1);
    for i=1:size(current_X,1)
        for j=1:size(current_X,2)
            alpha = 0.10;%alpha
            beta = 0.15;%beta
            %------------Eq.(4)
            num_AK=log(1.0/Phi(0,1))*(global_Best_position(randi(gbestNum),j)-current_X(i,j));
            %------------Eq.(5)
            num_BK=alpha*Phi(0,1)*(ave_gbest(j)-current_X(i,j));
            %------------Eq.(6)
            num_CK=beta*Phi(0,1)*(ave_Pbest(j)-current_X(i,j));
            %------------Eq.(3)
            next_X(i,j)=current_X(i,j)+num_AK+num_BK+num_CK;         
        end     
        %boundary check
        next_X(i,:)=BoundaryCheck(next_X(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
    end
    reflect_X=zeros(pop,dim);
    for j=1:size(current_X,2)
        for i=1:size(current_X,1)
            %------------Eq.(10)
            num_C= (ub(j) + lb(j)) * 0.5;
            gailv=abs(num_C- next_X(i,j))/(ub(j) - lb(j));
            if  next_X(i,j)>=num_C  %------------Eq.(7) First
                if gailv<Phi(0,1)
                    reflect_X(i,j)= Phi((ub(j) + lb(j))- next_X(i,j),num_C)  ;%------------Eq.(8) First
                else
                    reflect_X(i,j)= Phi(lb(j),(ub(j) + lb(j))- next_X(i,j))  ;%------------Eq.(8) Second
                end
            else%------------Eq.(7) Second
                if gailv<Phi(0,1)
                    reflect_X(i,j)= Phi(num_C,(ub(j) + lb(j))- next_X(i,j))  ;%------------Eq.(9) First
                else
                    reflect_X(i,j)= Phi((ub(j) + lb(j))- next_X(i,j),ub(j))  ;%------------Eq.(9) Second
                end
            end          
        end % end j
    end

    reflect_X(i,:)=BoundaryCheck(reflect_X(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
    for i=1:size(current_X,1)
        reflect_X_fit=fobj( reflect_X(i,:));
        next_X_fit=fobj( next_X(i,:));
        if next_X_fit>reflect_X_fit
            current_X(i,:)=reflect_X(i,:);
        else
            current_X(i,:)=next_X(i,:);
        end
    end

    for i=1:size(current_X,1)
        current_X_fitness(1,i)=fobj(current_X(i,:));
        % Update the pbest
        if current_X_fitness(1,i)<Person_Best_fitness(1,i)
            Person_Best_fitness(1,i)=current_X_fitness(1,i);
            Person_Best_position(i,:)=current_X(i,:);
        end
    end

    tem_X=zeros(gbestNum*2,dim);
    tem_X_fitness=inf(1,gbestNum*2);
    [~,index_sorted]=sort(current_X_fitness);
    for i=1:gbestNum
        tem_X(i,:)=current_X(index_sorted(i),:);
        tem_X_fitness(i)=current_X_fitness(index_sorted(i));
    end

    tem_X=[tem_X;global_Best_position];
    tem_X_fitness=[tem_X_fitness,global_Best_fitness];
    [~,index_sorted]=sort(tem_X_fitness);
    global_Best_position=tem_X(index_sorted([1:gbestNum]),:);
    global_Best_fitness=tem_X_fitness(index_sorted([1:gbestNum]));
    [~,index_sorted]=sort(global_Best_fitness);

    Destination_position=global_Best_position(index_sorted(1),:);
    Destination_fitness=global_Best_fitness(index_sorted(1));
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Destination_fitness
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index_sorted(1);
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Destination_fitness;
    end
    IterCurve(iter)=Destination_fitness;
    iter=iter+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=iter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function o=Phi(num1,num2)
    if num1<num2
        o=num1+rand()*abs(num2-num1);
    else
        o=num2+rand()*abs(num2-num1);
    end
end