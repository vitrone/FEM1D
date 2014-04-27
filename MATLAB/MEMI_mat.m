clc
pp = 10;
syms i;
% MEMI = sym('MEMI',[pp+1 pp+1]);
MEMI = sym(zeros(pp+1,pp+1));
M1 = sym([         2/3             1/3   -1/sqrt(6)   1/sqrt(10)/3;
              1/3             2/3   -1/sqrt(6)  -1/sqrt(10)/3;
       -1/sqrt(6)      -1/sqrt(6)          2/5              0;
     1/sqrt(10)/3   -1/sqrt(10)/3            0           2/21;
     ]);

MEMI(1:4,1:4) = M1;
i = 2;
MEMI(i+1,i+1)   = 2/((2*i+1)*(2*i-3));
MEMI(i+1,i+3)   = -1/((2*i+1)*sqrt((2*i-1)*(2*i+3)));
i = 3;
MEMI(i+1,i+1)   = 2/((2*i+1)*(2*i-3));
MEMI(i+1,i+3)   = -1/((2*i+1)*sqrt((2*i-1)*(2*i+3)));
for i=4:pp-2,
    MEMI(i+1,i-1)   = -1/((2*i-3)*sqrt((2*i-1)*(2*i-5))); %#ok<*SPRIX>
    MEMI(i+1,i+1)   = 2/((2*i+1)*(2*i-3));
    MEMI(i+1,i+3)   = -1/((2*i+1)*sqrt((2*i-1)*(2*i+3)));
end
for i=pp-1:pp,
    MEMI(i+1,i-1)   = -1/((2*i-3)*sqrt((2*i-1)*(2*i-5)));
    MEMI(i+1,i+1)   = 2/((2*i+1)*(2*i-3));
end

%%
dlmwrite('MEMI.dat', double(MEMI), 'delimiter', ',', 'precision', '% 0.20f');
%ccode(MEMI, 'file', 'MEMI.dat')
%spy(double(MEMI));
