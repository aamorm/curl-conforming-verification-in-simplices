% The gradient functions are included in the order k+1, e.g.,
% G1 belong to the second-order functions.
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
detailed_output = true;

% First order functions
disp('Order 1')
disp('Edge functions')
homogeneous_order = 1;
W1_function = -L2*gL1+L1*gL2;
good_function = check_space(W1_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("W1 function is in Nedelec space")
else
    disp("W1 function is not in Nedelec space")
end

% Second order functions
disp('Order 2')
homogeneous_order = 2;

disp('Edge functions')
G1_function = L2*gL1+L1*gL2;
good_function = check_space(G1_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G1 function is in Nedelec space")
else
    disp("G1 function is not in Nedelec space")
end

disp('Face functions')
R20_function_1 = L2*L3*gL1 + L1*L3*gL2 - 2*L1*L2*gL3;
good_function = check_space(R20_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First R20 function is in Nedelec space")
else
    disp("First R20 function is not in Nedelec space")
end
R20_function_2 = L1*L3*gL2 + L1*L2*gL3 - 2*L2*L3*gL1;
good_function = check_space(R20_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second R20 function is in Nedelec space")
else
    disp("Second R20 function is not in Nedelec space")
end

detailed_output = true;
% Third order functions
disp('Order 3')
homogeneous_order = 3;

disp('Edge functions')
G2_function = (2*L1-L2)*L2*gL1 - (2*L2-L1)*L1*gL2;
good_function = check_space(G2_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G2 function is in Nedelec space")
else
    disp("G2 function is not in Nedelec space")
end

disp('Face functions')
G20_function = L2*L3*gL1 + L1*L3*gL2 + L1*L2*gL3;
good_function = check_space(G20_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G20 function is in Nedelec space")
else
    disp("G20 function is not in Nedelec space")
end
R30_function_1 = L1*L2*(L1-L2)*gL3;
good_function = check_space(R30_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First non-orthogonalized R30 function is in Nedelec space")
else
    disp("First non-orthogonalized R30 function is not in Nedelec space")
end
R30_function_2 = L2*L3*(L2-L3)*gL1;
good_function = check_space(R30_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second non-orthogonalized R30 function is in Nedelec space")
else
    disp("Second non-orthogonalized R30 function is not in Nedelec space")
end
R30_function_3 = L3*L1*(L3-L1)*gL2;
good_function = check_space(R30_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third non-orthogonalized R30 function is in Nedelec space")
else
    disp("Third non-orthogonalized R30 function is not in Nedelec space")
end

% Orthogonalized version. You can generate the other three commuting the indices.
R30_function_orth_1 = -L2*L3*(1398*L1-699*L2+86)*gL1 + L1*L3*(1398*L2-699*L1+86)*gL2 + 2020*L1*L2*(L1-L2)*gL3;
good_function = check_space(R30_function_orth_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Orthogonalized R30 function is in Nedelec space")
else
    disp("Orthogonalized R30 function is not in Nedelec space")
end

disp('Volume functions')
R300_function = L1*L2*L3*gL1+L1*L2*L3*gL2+L1*L2*L3*gL3;
good_function = check_space(R300_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Non-orthogonalized R300 function is in Nedelec space")
else
    disp("Non-orthogonalized R300 function is not in Nedelec space")
end
R301_function = L2*L3*L4*gL1;
good_function = check_space(R301_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Non-orthogonalized R301 function is in Nedelec space")
else
    disp("Non-orthogonalized R301 function is not in Nedelec space")
end
R302_function = L3*L1*L4*gL2;
good_function = check_space(R302_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Non-orthogonalized R302 function is in Nedelec space")
else
    disp("Non-orthogonalized R302 function is not in Nedelec space")
end

R300_function = L2*L3*(L4+3*L1)*gL1+L1*L3*(L4+3*L2)*gL2+L1*L2*(L4+3*L3)*gL3;
good_function = check_space(R300_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Orthogonalized R300 function is in Nedelec space")
else
    disp("Orthogonalized R300 function is not in Nedelec space")
end
R301_function = 2*L2*L3*L4*gL1-L3*L1*L4*gL2-L1*L2*L4*gL3;
good_function = check_space(R301_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Orthogonalized R301 function is in Nedelec space")
else
    disp("Orthogonalized R301 function is not in Nedelec space")
end
R302_function = L3*L1*L4*gL2-L1*L2*L4*gL3;
good_function = check_space(R302_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Orthogonalized R302 function is in Nedelec space")
else
    disp("Orthogonalized R302 function is not in Nedelec space")
end