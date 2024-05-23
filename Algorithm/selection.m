function [dad,mom] = selection(pop,fitvalue)
select={'轮盘赌','排序选择','精英保存','锦标赛'}; 
switch select{3}
    case '轮盘赌'
        PP = cumsum( fitvalue ./ sum(fitvalue) );       
        [row, ~] = size(pop);
        for i = 1:row
            for j = 1:row
                r = rand;
                if r <= PP(j)
                    dad(i,:) = pop(j,:);
                    break;
                end
            end
            mom(i,:) = pop(randi([1 row]),:);
        end
    case '排序选择'
        [sort_value,sort_index]=sort(fitvalue); 
        pop=pop(sort_index,:);                 
        q=0.3;                                  
        p=q*(1-q).^((1:size(pop,1)-1))/(1-(1-q)^size(pop,1));
        PP=cumsum(p);
        for i = 1:size(pop,1)
            for j = 1:size(pop,1)
                r = rand;
                if r <= PP(j)
                    dad(i,:) = pop(j,:);
                    break;
                end
            end
            
            mom(i,:) = pop(randi([1 size(pop,1)]),:);
        end
    case '精英保存'
       
        PP = cumsum( fitvalue ./ sum(fitvalue) );       
        [row, ~] = size(pop);
        
        for i = 1:row
            
            for j = 1:row
                r = rand;
                if r <= PP(j)
                    dad(i,:) = pop(j,:);
                    break;
                end
            end
            
            mom(i,:) = pop(randi([1 row]),:);
        end
    case '锦标赛'
        [row, ~] = size(pop);
        for i = 1:row
            sn=randi(row);
            r = randi([1 row],sn,1);
            [~,bestindex] = min(fitvalue(r));   
            dad(i,:) = pop(r(bestindex),:);
            mom(i,:) = pop(randi([1 row]),:);
        end
end    