function [ AngX] = AngX( AB )
%   ����������X��ĽǶ�
%%������������2X1��ʽ
                C = [1;0];
                A = [0;0];
                B = AB;
                if AB(2) == 0
     if AB(1) > 0
        AngX = 0;
    else
        AngX = 180;
     end
    elseif AB(2) > 0
    AngX = acosd((norm(A-B)^2+norm(A-C)^2- norm(B-C)^2)/(2*(norm(A-B)*norm(A-C))));
    elseif AB(2) < 0
    AngX = 180+acosd((norm(A-B)^2+norm(A+C)^2- norm(B+C)^2)/(2*(norm(A-B)*norm(A+C))));
                end
    if AngX>270
     AngX= AngX-360;
    end
%%%%%%%%%%%%%%%%%%%%%����Ϊ�����򣨽����AngBAX--AB��X����нǣ�
end

