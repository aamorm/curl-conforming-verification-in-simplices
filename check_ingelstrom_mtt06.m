% Only the rotational part is checked.
% See notes in 17.03.2019
syms L1 L2 L3 L4
assume(L1, 'real');
assumeAlso(L1 >= 0 & L1 <= 1);
assume(L2, 'real');
assumeAlso(L2 >= 0 & L2 <= 1);
assume(L3, 'real');
assumeAlso(L3 >= 0 & L3 <= 1);
assume(L4, 'real');
assumeAlso(L4 >= 0 & L4 <= 1);

affine_coordinates = [L1, L2, L3, L4];
gL1 = gradient(L1, affine_coordinates);
gL2 = gradient(L2, affine_coordinates);
gL3 = gradient(L3, affine_coordinates);
gL4 = gradient(L4, affine_coordinates);
tolerance = 1e-10;
detailed_output = false;

% First order functions
disp('Order 1')
disp('Rotational edge functions')
homogeneous_order = 1;
A1_function = -L2*gL1+L1*gL2;
good_function = check_space(A1_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("A1 function is in Nedelec space")
else
    disp("A1 function is not in Nedelec space")
end

% Second order functions
disp('Order 2')
homogeneous_order = 2;

disp('Rotational face functions')
A2_function_1 = 3*L2*L3*gL1 - gradient(L1*L2*L3,affine_coordinates);
good_function = check_space(A2_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First A2 function is in Nedelec space")
else
    disp("First A2 function is not in Nedelec space")
end
A2_function_2 = 3*L3*L1*gL2 - gradient(L1*L2*L3,affine_coordinates);
good_function = check_space(A2_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second A2 function is in Nedelec space")
else
    disp("Second A2 function is not in Nedelec space")
end

% Third order functions
disp('Order 3')
homogeneous_order = 3;

disp('Rotational face functions')
A3_function_1 = 4*L2*L3*(L2-L3)*gL1 - gradient(L1*L2*L3*(L2-L3),affine_coordinates);
good_function = check_space(A3_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First A3 function is in Nedelec space")
else
    disp("First A3 function is not in Nedelec space")
end
A3_function_2 = 4*L3*L1*(L3-L1)*gL2 - gradient(L2*L3*L1*(L3-L1),affine_coordinates);
good_function = check_space(A3_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second A3 function is in Nedelec space")
else
    disp("Second A3 function is not in Nedelec space")
end
A3_function_3 = 4*L1*L2*(L1-L2)*gL3 - gradient(L3*L1*L2*(L1-L2),affine_coordinates);
good_function = check_space(A3_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third A3 function is in Nedelec space")
else
    disp("Third A3 function is not in Nedelec space")
end

disp('Rotational volume functions')
A3_function_1 = 4*L2*L3*L4*gL1 - gradient(L1*L2*L3*L4,affine_coordinates);
good_function = check_space(A3_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First A3 function is in Nedelec space")
else
    disp("First A3 function is not in Nedelec space")
end
A3_function_2 = 4*L3*L4*L1*gL2 - gradient(L1*L2*L3*L4,affine_coordinates);
good_function = check_space(A3_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second A3 function is in Nedelec space")
else
    disp("Second A3 function is not in Nedelec space")
end
A3_function_3 = 4*L4*L1*L2*gL3 - gradient(L1*L2*L3*L4,affine_coordinates);
good_function = check_space(A3_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third A3 function is in Nedelec space")
else
    disp("Third A3 function is not in Nedelec space")
end

% Fourth order functions
disp('Order 4')
homogeneous_order = 4;

disp('Rotational face functions')
A4_function_1 = 5*L2*L3*(L2^2-3*L2*L3+L3^2)*gL1 - gradient(L1*L2*L3*(L2^2-3*L2*L3+L3^2),affine_coordinates);
good_function = check_space(A4_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First A4 function is in Nedelec space")
else
    disp("First A4 function is not in Nedelec space")
end
A4_function_2 = 5*L3*L1*(L3^2-3*L3*L1+L1^2)*gL2 - gradient(L1*L2*L3*(L3^2-3*L3*L1+L1^2),affine_coordinates);
good_function = check_space(A4_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second A4 function is in Nedelec space")
else
    disp("Second A4 function is not in Nedelec space")
end
A4_function_3 = 5*L1*L2*(L1^2-3*L1*L2+L2^2)*gL3 - gradient(L1*L2*L3*(L1^2-3*L1*L2+L2^2),affine_coordinates);
good_function = check_space(A4_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third A4 function is in Nedelec space")
else
    disp("Third A4 function is not in Nedelec space")
end
A4_function_4 = (6*L1-L2-L3)*(L2-L3)*L2*L3*gL1 + (6*L2-L3-L1)*(L3-L1)*L3*L1*gL2 + (6*L3-L1-L2)*(L1-L2)*L1*L2*gL3;
good_function = check_space(A4_function_4, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Fourth A4 function is in Nedelec space")
else
    disp("Fourth A4 function is not in Nedelec space")
end

disp('Rotational volume functions')
A4_function_1 = 5*L2*L3*L4*(L2-L3)*gL1 - gradient(L1*L2*L3*L4*(L2-L3),affine_coordinates);
good_function = check_space(A4_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First A4 function is in Nedelec space")
else
    disp("First A4 function is not in Nedelec space")
end
A4_function_2 = 5*L2*L3*L4*(L3-L4)*gL1 - gradient(L1*L2*L3*L4*(L3-L4),affine_coordinates);
good_function = check_space(A4_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second A4 function is in Nedelec space")
else
    disp("Second A4 function is not in Nedelec space")
end
A4_function_3 = 5*L3*L4*L1*(L3-L4)*gL2 - gradient(L1*L2*L3*L4*(L3-L4),affine_coordinates);
good_function = check_space(A4_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third A4 function is in Nedelec space")
else
    disp("Third A4 function is not in Nedelec space")
end
A4_function_4 = 5*L3*L4*L1*(L4-L1)*gL2 - gradient(L1*L2*L3*L4*(L4-L1),affine_coordinates);
good_function = check_space(A4_function_4, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Fourth A4 function is in Nedelec space")
else
    disp("Fourth A4 function is not in Nedelec space")
end
A4_function_5 = 5*L4*L1*L2*(L4-L1)*gL3 - gradient(L1*L2*L3*L4*(L4-L1),affine_coordinates);
good_function = check_space(A4_function_5, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Fifth A4 function is in Nedelec space")
else
    disp("Fifth A4 function is not in Nedelec space")
end
A4_function_6 = 5*L4*L1*L2*(L1-L2)*gL3 - gradient(L1*L2*L3*L4*(L1-L2),affine_coordinates);
good_function = check_space(A4_function_6, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Sixth A4 function is in Nedelec space")
else
    disp("Sixth A4 function is not in Nedelec space")
end
A4_function_7 = 5*L1*L2*L3*(L1-L2)*gL4 - gradient(L1*L2*L3*L4*(L1-L2),affine_coordinates);
good_function = check_space(A4_function_7, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Seventh A4 function is in Nedelec space")
else
    disp("Seventh A4 function is not in Nedelec space")
end
A4_function_8 = 5*L1*L2*L3*(L2-L3)*gL4 - gradient(L1*L2*L3*L4*(L2-L3),affine_coordinates);
good_function = check_space(A4_function_8, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Eigth A4 function is in Nedelec space")
else
    disp("Eigth A4 function is not in Nedelec space")
end
