function Data_o = ESDA(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
options.x = Data_i.X;                        
options.F=Data_i.F_value;             
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);

options.d=Data_i.dim;                      % dimension
options.lb=Data_i.lb;                      % lower bound
options.ub=Data_i.ub;                      % upper bound
options.ProblemSize=length(options.ub);    % dimension of the problem.
options.ObjectsSize=Data_i.pop;            % m: number of objects
options.MaxIter=Data_i.maxIter;            % MAXITER: maximum number of iterations
options.ObjFunction=Data_i.fobj;           % the name of the objective function
options.Display_Flag=1;                    % Flag for displaying results over iterations
options.run_parallel_index=0;              % 1 for parallel processing
options.run=1;

[Xbest,Fbest,IterCurve,index]=ESDA_v1(options);
if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>Fbest
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=Fbest; 
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index(1);
end
% if options.run_parallel_index
%     stream = RandStream('mrg32k3a');
%     parfor index=1:options.run
%         set(stream,'Substream',index);
%         RandStream.setGlobalStream(stream)
%         
%         bestX_M(index,:)=Xbest;
%         Fbest_M(index)=Fbest;
%         IterCurve=FunctionEvolution_best;
%     end
% else
%     rng('default')
%     for index=1:options.run
%         [Xbest, Fbest,FunctionEvolution_best]=ESDA_v1(options);
%         bestX_M(index,:)=Xbest;
%         Fbest_M(index)=Fbest;
%         IterCurve=FunctionEvolution_best;
%     end
% end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=options.MaxIter;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end

function [Xbest,Fbest,FunctionEvolution_best,index]=ESDA_v1(options)
%--------------------------------------------------------------------------
% ESDA Electrostatic Discharge Algorithm
%--------------------------------------------------------------------------
FunctionEvolution_best=[];
ub=options.ub;
lb=options.lb;
ProblemSize=options.ProblemSize;
ObjFunction=options.ObjFunction;
ObjectsSize = options.ObjectsSize;
MaxIter = options.MaxIter;
archSize=round(ObjectsSize/3);

x=options.x;

x(:,ProblemSize+1)=0;
F=options.F;
[F, ind] = sort(F);
x = x(ind,:);
cc=2;

for Iter=1:MaxIter
    x_arch  = x(1:archSize,:);
    F_arch=F(1:archSize);
    newobj=[];
    
    for i=1:2*archSize
        SP=sort(randperm(archSize,3));
        ii=SP(1);
        jj=SP(2);
        kk=SP(3);
        if rand>0.5
            newobj(i,1:ProblemSize)= x(jj,1:ProblemSize) + cc*normrnd(0.7,0.2) * (x(ii,1:ProblemSize) - x(jj,1:ProblemSize));
            x(jj,ProblemSize+1)=x(jj,ProblemSize+1)+1;
        else
            newobj(i,1:ProblemSize)= x(kk,1:ProblemSize) + cc*normrnd(0.7,0.2) * (x(ii,1:ProblemSize) - x(kk,1:ProblemSize))+...
                cc*normrnd(0.7,0.2) * (x(jj,1:ProblemSize) - x(kk,1:ProblemSize));
            x(jj,ProblemSize+1)=x(jj,ProblemSize+1)+1;
        end        
    end
    
    newobj =BoundaryCheck(newobj(:,1:ProblemSize),ub,lb,options.d);
    newobj(i,ProblemSize+1)=x(jj,ProblemSize+1)+1;
    
    for i1=1:ObjectsSize
        if (x(i1,ProblemSize+1)>=1)&&(x(i1,ProblemSize+1)<=3)
            for j=1:ProblemSize
                if (rand<0.2)
                    pos= randi(archSize(1));
                    newobj(i1,j)=  x_arch(pos,j);
                    x(i1,ProblemSize+1)= 0;
                end
            end
        elseif (x(i1,ProblemSize+1)>3)
            newobj(i1,1:ProblemSize)=lb+rand(1,ProblemSize).*(ub-lb);
            x(i1,ProblemSize+1)= 0;
        end
    end
    
    x_arch  = x(1:archSize,:);
    
    x_all=[x_arch; newobj];
    F_all=[];
    F_all=[F_all F_arch];
    for i=1:size(newobj(:,1:ProblemSize),1)
        F_all(length(F_all)+1)=feval(ObjFunction,newobj(i,1:ProblemSize));
    end        
    [F_all, index ] = sort(F_all);
    x_all = x_all( index, : );
    
    x = x_all (1:ObjectsSize,:);
    F = F_all(1:ObjectsSize);
    
    Xbest=x(1,1:ProblemSize);
    Fbest=F(1);
    %--------------------------------------------------------------------------------
    FunctionEvolution_best(Iter)=Fbest;
    %--------------------------------------------------------------------------------
       
%     if options.Display_Flag==1
%         fprintf('Iteration N° is %g Best Fitness is %g\n',Iter,Fbest)
%     end
end
Fbest=F(1);
end

