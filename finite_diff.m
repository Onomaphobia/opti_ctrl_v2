% size = 10;
% M = eye(size);
% i = 1;
% while i<=size
%     j = 1;
%     while j<=size
%         if j == i - 1
%             M(i,j) = -1;
%         elseif j == i
%             M(i,j) = 2;
%         elseif j == i + 1
%             M(i,j) = -1;
%         else
%             M(i,j) = 0;
%         end
%         j = j + 1;
%     end
%     i = i + 1;
% end
% M

Diff_rules(1, 1) = 0;
Diff_rules(1, 2) = -1 / 12.0;
Diff_rules(1, 3) = 16 / 12.0;
Diff_rules(1, 4) = -30 / 12.0;
Diff_rules(1, 5) = 16 / 12.0;
Diff_rules(1, 6) = -1 / 12.0;
Diff_rules(1, 7) = 0;
size = 15;
i = 0;
M = zeros(size);
while i <=size
    i = i + 1;
    j = -3;
    while j < 3
        j = j + 1;
        index = i + j;
        if index <= 0
            continue
        end
        if index >size
            continue
        end
        M(i,index) = Diff_rules(1,j + 4);
        
    end

end
           
            






