function [v b]=LP2FEM_ShapeFunc( p, N, u)

if p<2,
    b = [];
    v = zeros(N+1,1);
    v(1) = -u(2)+u(1);
    v(2) = u(2)+u(1);
    
    for i=1:N-1,
        v(i+2) = v(i+1)+2*u(i*2+2,1);
    end
    
elseif p==2,
    v = zeros(N+1,1);
    b = zeros(N,1);
    
    b(1) = sqrt(6)*u(3);
    v(2) = u(3)+b(2)/sqrt(10)+u(2)+u(1);
    v(1) = u(3)-b(2)/sqrt(10)-u(2)+u(1);
    
    for i=1:N-1,
        v(i+2) = v(i+1)+2*u(i*3+2);
    end
    
else
    v = zeros(N+1,1);
    b = zeros(N*(p-1),1);
    tmpA = [sqrt(2*(2*p-1)), sqrt(2*(2*p-3)), 1/sqrt(6), 1/sqrt(10), 2/sqrt(10)];
    b(p-1) = tmpA(1)*u(p+1);
    b(p-2) = tmpA(2)*u(p);
    tmpB = zeros(1, p-3);
    tmpC = zeros(1, p-3);
    
    for j=2:p-2,
        tmpB(j-1) = sqrt(2*(2*(p-j)-1));
        tmpC(j-1) = tmpB(j-1)/sqrt(2*(2*(p-j)+3));
        b(p-1-j) = tmpB(j-1)*u(p-j+1)+tmpC(j-1)*b(p-j+1);
    end
    v(2) = tmpA(3)*b(1)+tmpA(4)*b(2)+u(2)+u(1);
    v(1) = tmpA(3)*b(1)-tmpA(4)*b(2)-u(2)+u(1);
    
    for i=1:N-1,
        ptmp1 = i*(p-1);
        ptmp2 = i*(p+1);
        b(ptmp1+p-1) = tmpA(1)*u(ptmp2+p+1);
        b(ptmp1+p-2) = tmpA(2)*u(ptmp2+p);
        for j=2:p-2,
            b(ptmp1+p-1-j) = tmpB(j-1)*u(ptmp2+p-j+1)+tmpC(j-1)*b(ptmp1+p-j+1);
        end
        v(i+2) = v(i+1)+tmpA(5)*b(ptmp1+2)+2*u(ptmp2+2,1);
    end
    
end

