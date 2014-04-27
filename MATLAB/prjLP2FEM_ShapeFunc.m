function [Pv Pb] = prjLP2FEM_ShapeFunc( p, N, u)

Pv = zeros( N+1, 1);
Pb = zeros( N*(p-1), 1);

Pv(1) = u(1)-u(2)/3;
for i=1:N-1,
    Pv(i+1) = u((p+1)*(i-1)+1)+u((p+1)*(i-1)+2)/3+...
        +u((p+1)*i+1)-u((p+1)*i+2)/3;
end
Pv(N+1) = u((p+1)*(N-1)+1)+u((p+1)*(N-1)+2)/3;
if p==2,
    tmpA = [ -(2/sqrt(6)), (2/sqrt(6)/5)];
    for i=0:N-1,
        ptmp = 3*i+1;
        Pb(i+1) = tmpA(1)*u(ptmp)+tmpA(2)*u(ptmp+2);
    end
elseif p==3,
    tmpA = [-(2/sqrt(6)), (2/sqrt(6)/5), -(2/sqrt(10)/3), (2/sqrt(10)/7)];
    for i=0:N-1,
        ptmp1 = i*2+1; 
        ptmp2 = i*4+1; 
        Pb(ptmp1) = tmpA(1)*u(ptmp2)+tmpA(2)*u(ptmp2+2);
        Pb(ptmp1+1) = tmpA(3)*u(ptmp2+1)+tmpA(4)*u(ptmp2+3);
    end
elseif p>3,
    for i=0:N-1,
        ptmp1 = i*(p-1);
        ptmp2 = i*(p+1);
        tmpB = zeros( p-1, 1);
        tmpC = zeros( p-1, 1);
        for j=1:p-1,
            tmpB(j) = 2/(sqrt(2*(2*j+1)))/(2*j+3);
            tmpC(j) = -2/(sqrt(2*(2*j+1)))/(2*j-1);
        end
        for j=1:p-1,
            Pb(ptmp1+j) =  tmpB(j)*u(ptmp2+j+2)+tmpC(j)*u(ptmp2+j);
        end
    end
end