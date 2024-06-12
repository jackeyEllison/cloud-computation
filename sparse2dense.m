% 将稀疏矩阵Rup转换为稠密矩阵Rup_new 


function dense = sparse2dense(matrix)
    M = size(matrix,1); % rows
    dense = zeros(M,1);
    for m=1:M
        % 找到最大值的位置 而不是找到不为0的那个位置
        % 输入的矩阵有可能在m行全为0，如果找不为0的位置有bug
        [~,index] = max(matrix(m,:)); 
        dense(m) = matrix(m,index);
    end   
end