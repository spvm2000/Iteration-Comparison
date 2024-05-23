function [X] = BoundaryCheck(x,ub,lb,dim)
    if size(x,1)>1
        for i=1:size(x,1)
            for j=1:size(x,2)
                if x(i,j)>ub(j)
                   x(i,j) = ub(j); 
                end
                if x(i,j)<lb(j)
                    x(i,j) = lb(j);
                end
                if isnan(x(i,j)) || isinf(x(i,j))
                    x(i,j)=rand*(ub(j)-lb(j));
                end    
            end    
        end
    else
        for i = 1:dim
            if x(i)>ub(i)
               x(i) = ub(i); 
            end
            if x(i)<lb(i)
                x(i) = lb(i);
            end
            if isnan(x(i)) || isinf(x(i))
                x(i)=rand*(ub(i)-lb(i));
            end
        end
    end
    X = x;
    %dimΪ���ݵ�ά�ȴ�С
    %xΪ�������ݣ�ά��Ϊ[1,dim];
    %ubΪ�����ϱ߽磬ά��Ϊ[1,dim]
    %lbΪ�����±߽磬ά��Ϊ[1,dim]
end