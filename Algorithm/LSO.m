function Data_o = LSO(Data_i,Data_o)
rand_num=[];                         
ti=clock;                             
LightRays = Data_i.X;                       
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fitness);
IterCurve=zeros(1,Data_i.maxIter);
GBest_fitness=Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter);
index=Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter);

GBestRayColor=LightRays(index,:);    
NewWLightRays=LightRays;
red=1.3318; violet=1.3435;    % From assumption (2)
it=1; 
Ps=0.05; %% Probaility of first and second scattering stages
Pe=0.6;  %% Contolling parameter to exchange between the first and second scattering 
Ph=0.4;  %% Pobability of Hybridization between two various search boundary methods 
B=0.05;  %% Exploitation probability in the first scattering stage

Max_iter=Data_i.maxIter;                
NoF_LightRays=Data_i.pop;
dim=Data_i.dim;
lb=Data_i.lb;
ub=Data_i.ub;
while it<=Max_iter
    for i=1:NoF_LightRays
        nA=LightRays(randi(NoF_LightRays),:);        
        nB= LightRays(i,:); 
        nC=GBestRayColor; 
        xbar=(sum(LightRays)/NoF_LightRays);
        norm_nA=nA/norm(nA);   % the normal vector of inner refraction  "equation (6)"
        norm_nB=nB/norm(nB);   % the normal vector of inner reflection  "equation (7)"
        norm_nC=nC/norm(nC);   % the normal vector of outer refraction  "equation (8)"

        Incid_norm=xbar/norm(xbar);
        k=red+rand.*(violet-red);
        p=rand; 
        q=rand; 
        L1=(1./k).*(Incid_norm-(norm_nA.*dot(norm_nA,Incid_norm)))-(norm_nA.*(abs((1-(1./k.^2)+((1./k.^2).*dot(norm_nA,Incid_norm).^2)))).^(1/2)); % "equaton (12)"
        L2=L1-((2.*norm_nB).*dot(L1,norm_nB)); % "equation (13)" 
        L3=k.*(L2-(norm_nC.*dot(norm_nC,L2)))+norm_nC.*(abs(1-(k.^2)+(k.^2).*((dot(norm_nC,L2)).^2))).^(1/2); % "equation (14)" 
       
        a=rand*(1-it/Max_iter) ;    % "equation (20)" 
        ginv = gammaincinv(a,1);    % compute the new ginv value
        GI=a*(1/rand)*ginv;         %"equation (19)"
        Epsln=a.*randn(1,dim);      %  "equation (18)" 

        if p<=q
            NewWLightRays(i,:)=LightRays(i,:)+GI.*Epsln.*rand(1,dim).*(L1-L3).*(LightRays(randi(NoF_LightRays),:)-LightRays(randi(NoF_LightRays),:));   % "equation (17)"  
        else
            NewWLightRays(i,:)=(LightRays(i,:))+GI.*Epsln.*rand(1,dim).*(L2-L3).*(LightRays(randi(NoF_LightRays),:)-LightRays(randi(NoF_LightRays),:));   % "equation (18)"
        end

        if rand<Ph
            U=NewWLightRays(i,:)>ub;
            L=NewWLightRays(i,:)<lb;
            NewWLightRays(i,:)=(NewWLightRays(i,:).*(~(U+L)))+ub.*U+lb.*L;
        else
            for j=1:size(LightRays,2)
               if  NewWLightRays(i,j)>ub(j)
                   NewWLightRays(i,j)=lb(j)+rand*(ub(j)-lb(j));
               elseif  NewWLightRays(i,j)<lb(j)
                   NewWLightRays(i,j)=lb(j)+rand*(ub(j)-lb(j));
               end
            end
        end
        % Calculate objective function for each search agent
        Fnew=Data_i.fobj(NewWLightRays(i,:));
        %If fitness improves (better solutions found), update then
        if (Fnew<=fitness(i)),LightRays(i,:)=NewWLightRays(i,:);
             fitness(i)=Fnew;
        end

        %% update Global best solution
        if  Fnew<=GBest_fitness
            GBestRayColor=NewWLightRays(i,:);    % Update Global best solution
            GBest_fitness=Fnew;
            Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Fnew;
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
        end
        [in1]=sort(fitness);
        IterCurve(it)=GBest_fitness;
        it=it+1;
        if (it>Max_iter)
           break;
        end

        F=abs((fitness(i)-GBest_fitness)/(GBest_fitness-(in1(NoF_LightRays)))); %  "equation (25)" 
        if F<rand || rand<Ps
            if rand<Pe
                  NewWLightRays(i,:)=(LightRays(i,:))+(rand)*(LightRays(randi(NoF_LightRays),:)-LightRays(randi(NoF_LightRays),:))+(rand<B).*rand(1,dim).*((GBestRayColor - LightRays(i,:))); %  "equation (21)" 
            else          
                  NewWLightRays(i,:)=(((2*cos(rand*180)).*(GBestRayColor.*LightRays(i,:))));   % "equation (22)"
            end
        else
            U=(rand(1,dim)>rand(1,dim));
            NewWLightRays(i,:)= U.*(LightRays(randi(NoF_LightRays),:)+abs(randn).*(LightRays(randi(NoF_LightRays),:)-LightRays(randi(NoF_LightRays),:)))+(1-U).*LightRays(i,:);  % "equation (24)"
        end

        if rand<Ph
               U=NewWLightRays(i,:)>ub;
               L=NewWLightRays(i,:)<lb;
               NewWLightRays(i,:)=(NewWLightRays(i,:).*(~(U+L)))+ub.*U+lb.*L;
        else
            for j=1:size(LightRays,2)
            if  NewWLightRays(i,j)>ub(j)
                NewWLightRays(i,j)=lb(j)+rand*(ub(j)-lb(j));
            elseif  NewWLightRays(i,j)<lb(j)
                NewWLightRays(i,j)=lb(j)+rand*(ub(j)-lb(j));
            end
            end
        end

        Fnew=Data_i.fobj(NewWLightRays(i,:));
         % If fitness improves (better solutions found), update then
        if (Fnew<=fitness(i)),
            LightRays(i,:)=NewWLightRays(i,:);
            fitness(i)=Fnew;
        end
        %% update Global best solution
          if  Fnew<=GBest_fitness
              GBestRayColor=NewWLightRays(i,:);    % Update Global best solution
              GBest_fitness=Fnew;
              Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Fnew;
              Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=i;
          end
          IterCurve(it)=GBest_fitness;
          it=it+1;
          if (it>Max_iter)
             break;
          end
    end
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=it;                              
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end