function Data_o = WHO(Data_i,Data_o)
t_c=clock;                             
cnt=1;                                
%% Problem Definition
Max_iter = Data_i.maxIter;     
N = Data_i.pop;           
dim=Data_i.dim;
x = Data_i.X;             
fit = Data_i.F_value;         
lb = Data_i.lb;
ub = Data_i.ub;

PS=0.2;    
PC=0.13;    
NStallion=ceil(PS*N); 
Nfoal=N-NStallion;

Convergence_curve = zeros(1,Max_iter);

empty.pos=[];
empty.cost=[];

group=repmat(empty,Nfoal,1);

for i=1:Nfoal 
   group(i).pos=x(i,:);
   group(i).cost=fit(i);
end

Stallion=repmat(empty,NStallion,1);

for i=1:NStallion 
   Stallion(i).pos=x(i+Nfoal,:);
   Stallion(i).cost=fit(i+Nfoal);
 
end
  ngroup=length(group);
  a=randperm(ngroup);
  group=group(a);

i=0;
k=1;
for j=1:ngroup
    i=i+1;    
    Stallion(i).group(k)=group(j);  
    if i==NStallion
        i=0;
        k=k+1;
    end
end
Stallion=exchange(Stallion);
[~,index]=min([Stallion.cost]);
WH=Stallion(index); % global

[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(fit);
Convergence_curve(1)=WH.cost;
l=2; % Loop counter
while l<Max_iter+1
TDR=1-l*((1)/Max_iter);

for i=1:NStallion
    
   ngroup=length(Stallion(i).group);
    [~,index]=sort([Stallion(i).group.cost]);
    Stallion(i).group=Stallion(i).group(index);
   
   for j=1:ngroup
    
    if rand>PC
            z=rand(1,dim)<TDR;
            r1=rand;
            r2=rand(1,dim);
            idx=(z==0);
            r3=r1.*idx+r2.*~idx;
           rr=-2+4*r3;

           Stallion(i).group(j).pos= 2*r3.*cos(2*pi*rr).*(Stallion(i).pos-Stallion(i).group(j).pos)+(Stallion(i).pos);
    else
    A=randperm(NStallion);
    A(A==i)=[];
    a=A(1);
    c=A(2);

    x1=Stallion(c).group(end).pos;
    x2=Stallion(a).group(end).pos;

       y1=(x1+x2)/2;   % Crossover

    Stallion(i).group(j).pos=y1;
    end
   

    Stallion(i).group(j).pos=min(Stallion(i).group(j).pos,ub);
    Stallion(i).group(j).pos=max(Stallion(i).group(j).pos,lb);
    
    Stallion(i).group(j).cost=Data_i.fobj(Stallion(i).group(j).pos);
    
   end
    
% end
% 
% for i=1:NStallion

    R=rand;
%     z=rand(1,dim)<TDR;
%             r1=rand;
%             r2=rand(1,dim);
%             idx=(z==0);
%             r3=r1.*idx+r2.*~idx;
%             rr=-2+4*r3;

       if R<0.5
        k= 2*r3.*cos(2*pi*rr).*(WH.pos-(Stallion(i).pos))+WH.pos;
        else
        k= 2*r3.*cos(2*pi*rr).*(WH.pos-(Stallion(i).pos))-WH.pos;
       end

    k=min(k,ub);
    k=max(k,lb);
    fk=Data_i.fobj(k);
    if fk<Stallion(i).cost
      Stallion(i).pos  =k;
      Stallion(i).cost=fk;
    end
end
    Stallion=exchange(Stallion);
     [value,index]=min([Stallion.cost]);
   if value<WH.cost
       WH=Stallion(index);
   end
   gBestScore=WH.cost;
   Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
   Convergence_curve(l)=gBestScore;
   Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter) = gBestScore;
    l = l + 1;
end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);           
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=l;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=Convergence_curve;
end

function Stallion=exchange(Stallion)

nStallion=length(Stallion);
for i=1:nStallion
   [value,index]=min([Stallion(i).group.cost]);
   if value<Stallion(i).cost
       
       bestgroup=Stallion(i).group(index);
       
       Stallion(i).group(index).pos=Stallion(i).pos;
       Stallion(i).group(index).cost=Stallion(i).cost;
       
       
       Stallion(i).pos=bestgroup.pos;
       Stallion(i).cost=bestgroup.cost;
   end
end
end

