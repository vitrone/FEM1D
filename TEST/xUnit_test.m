clc
clear
set_pathvariables;
tc = xUnit_fem1d

run(tc, 'test_interpolation')
%run(tc, 'test_GP_ABC1a_CQ')
