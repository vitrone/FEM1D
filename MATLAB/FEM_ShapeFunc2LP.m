function u=FEM_ShapeFunc2LP(p, N, v, b)

u = zeros(N*(p+1),1);
if p<2,
    for i=0:N-1,
        u(i*2+1,1) = 0.5*(v(i+1)+v(i+2));
        u(i*2+2,1) = 0.5*(-v(i+1)+v(i+2));
    end
elseif p==2,
    tmpE = 1/sqrt(6);
    for i=0:N-1,
        ptmp = 3*i;
        u(ptmp+1,1) = 0.5*(v(i+1)+v(i+2))-tmpE*b(i+1);
        u(ptmp+2,1) = 0.5*(-v(i+1)+v(i+2));
        u(ptmp+3,1) = b(i+1)*tmpE;
    end
elseif p==3,
        tmpD = [1/sqrt(6), 1/sqrt(10)];
    for i=0:N-1,
        ptmp1 = 4*i;
        ptmp2 = 2*i;
        u(ptmp1+1) = 0.5*(v(i+1)+v(i+2))-b(ptmp2+1)*tmpD(1);
        u(ptmp1+2) = 0.5*(-v(i+1)+v(i+2))-b(ptmp2+2)*tmpD(2);
        u(ptmp1+3) = b(ptmp2+1)*tmpD(1);
        u(ptmp1+4) = b(ptmp2+2)*tmpD(2);
    end
else
    tmpA = [1/sqrt(6), 1/sqrt(10), 1/sqrt(2*(2*p-3)), 1/sqrt(2*(2*p-1))];
    tmpB = zeros(1,p-3);
    tmpC = zeros(1,p-3);
    for j=1:p-3,
            tmpB(j) = 1/sqrt(2*(2*j+1));
            tmpC(j) = -1/sqrt(2*(2*j+5));
    end
    for i=0:N-1,
        ptmp1 = i*(p+1);
        ptmp2 = i*(p-1);
        u(ptmp1+1) = 0.5*(v(i+1)+v(i+2))-b(ptmp2+1)*tmpA(1);
        u(ptmp1+2) = 0.5*(-v(i+1)+v(i+2))-b(ptmp2+2)*tmpA(2);
        for j=1:p-3,
            u(ptmp1+j+2) = b(ptmp2+j)*tmpB(j)+b(ptmp2+j+2)*tmpC(j);
        end
        u(ptmp1+p)    = b(ptmp2+p-2)*tmpA(3);
        u(ptmp1+p+1)  = b(ptmp2+p-1)*tmpA(4);
    end
end
