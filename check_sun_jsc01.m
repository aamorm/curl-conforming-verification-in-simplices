% These functions are located in tables 1,2,3.
% The gradient functions are included in the order k+1, e.g.,
% G1 belong to the second-order functions.
syms L0 L1 L2 L3
assume(L0, 'real');
assumeAlso(L0 >= 0 & L0 <= 1);
assume(L1, 'real');
assumeAlso(L1 >= 0 & L1 <= 1);
assume(L2, 'real');
assumeAlso(L2 >= 0 & L2 <= 1);
assume(L3, 'real');
assumeAlso(L3 >= 0 & L3 <= 1);

affine_coordinates = [L0, L1, L2, L3];
gL0 = gradient(L0, affine_coordinates);
gL1 = gradient(L1, affine_coordinates);
gL2 = gradient(L2, affine_coordinates);
gL3 = gradient(L3, affine_coordinates);
term_edge_L0 = L1*gL0;
term_edge_L1 = L0*gL1;
term_face_L0 = L1*L2*gL0;
term_face_L1 = L0*L2*gL1;
term_face_L2 = L0*L1*gL2;
tolerance = 1e-10;
detailed_output = true;

% First order functions
disp('Order 1')
disp('Edge functions')
homogeneous_order = 1;
R1_function = term_edge_L0-term_edge_L1;
good_function = check_space(R1_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("R1 function is in Nedelec space")
else
    disp("R1 function is not in Nedelec space")
end

% Second order functions
disp('Order 2')
homogeneous_order = 2;

disp('Edge functions')
G1_function = term_edge_L0+term_edge_L1;
good_function = check_space(G1_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G1 function is in Nedelec space")
else
    disp("G1 function is not in Nedelec space")
end

disp('Face functions')
R1_function_1 = term_face_L0+term_face_L1-2*term_face_L2;
good_function = check_space(R1_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First R1 function is in Nedelec space")
else
    disp("First R1 function is not in Nedelec space")
end
R1_function_2 = term_face_L0-term_face_L1;
good_function = check_space(R1_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second R1 function is in Nedelec space")
else
    disp("Second R1 function is not in Nedelec space")
end

% Third order functions
disp('Order 3')
homogeneous_order = 3;

disp('Edge functions')
G2_function = (2*L0-L1)*term_edge_L0-(2*L1-L0)*term_edge_L1;
good_function = check_space(G2_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G2 function is in Nedelec space")
else
    disp("G2 function is not in Nedelec space")
end

disp('Face functions')
G2_function = term_face_L0 + term_face_L1 + term_face_L2;
good_function = check_space(G2_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("G2 function is in Nedelec space")
else
    disp("G2 function is not in Nedelec space")
end
R3_function_1 = (L1-L2)*term_face_L0 - (L0-L2)*term_face_L1 + (L0-L1)*term_face_L2;
good_function = check_space(R3_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First R3 function is in Nedelec space")
else
    disp("First R3 function is not in Nedelec space")
end
R3_function_2 = (393*L0+80*L1-212*L2)*term_face_L0 - (393*L1+80*L0-212*L2)*term_face_L1 + (-292*L0+292*L1)*term_face_L2;
good_function = check_space(R3_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second R3 function is in Nedelec space")
else
    disp("Second R3 function is not in Nedelec space")
end
R3_function_3 = (-131*L0+168*L1-124*L2)*term_face_L0 + (-131*L1+168*L0-124*L2)*term_face_L1 + (-44*L0-44*L1+262*L2)*term_face_L2;
good_function = check_space(R3_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third R3 function is in Nedelec space")
else
    disp("Third R3 function is not in Nedelec space")
end

disp('Volume functions')
term_volume_L0 = L1*L2*L3*gL0;
term_volume_L1 = L0*L2*L3*gL1;
term_volume_L2 = L0*L1*L3*gL2;
term_volume_L3 = L0*L1*L2*gL3;
R3_function_1 = term_volume_L0 + term_volume_L1 - term_volume_L2 - term_volume_L3;
good_function = check_space(R3_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("First R3 function is in Nedelec space")
else
    disp("First R3 function is not in Nedelec space")
end
R3_function_2 = term_volume_L2 - term_volume_L3;
good_function = check_space(R3_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Second R3 function is in Nedelec space")
else
    disp("Second R3 function is not in Nedelec space")
end
R3_function_3 = term_volume_L0 - term_volume_L1;
good_function = check_space(R3_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("Third R3 function is in Nedelec space")
else
    disp("Third R3 function is not in Nedelec space")
end