function jac = jacobian_pattern()

%This function returns the Jacobian pattern for the clinical ODE, for use
%in simulating the ODE with ode15s().

%%% order of model states
% 1. pRAS
% 2. pALK
% 3. pERK
% 4. pAKT
% 5. pS6
% 6. Nprolif
% 7. Nnecrotic
% 8. Nnec1
% 9. Nnec2
% 10.Nnec3
% 11.Nnec4
% 12.Nkilled

%jac = zeros(15,15);
jac = zeros(12,12);

jac(1,1) = 1;
jac(1,2) = 1;
jac(2,2) = 1;
jac(3,1) = 1;
jac(3,3) = 1;
jac(4,1) = 1;
jac(4,4) = 1;
jac(5,3) = 1;
jac(5,4) = 1;
jac(5,5) = 1;
jac(6,4) = 1;
jac(6,5) = 1;
jac(6,6) = 1;
jac(6,7) = 1;
jac(6,8) = 1;
jac(6,9) = 1;
jac(6,10) = 1;
jac(6,11) = 1;
jac(7,4) = 1;
jac(7,5) = 1;
jac(7,6) = 1;
jac(7,7) = 1;
jac(7,8) = 1;
jac(7,9) = 1;
jac(7,10) = 1;
jac(7,11) = 1;
jac(8,7) = 1;
jac(8,8) = 1;
jac(9,8) = 1;
jac(9,9) = 1;
jac(10,9) = 1;
jac(10,10) = 1;
jac(11,10) = 1;
jac(11,11) = 1;
jac(12,4) = 1;
jac(12,5) = 1;
jac(12,6) = 1;
jac(12,12) = 1;

end