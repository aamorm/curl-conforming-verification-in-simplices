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

% From table 2. i,j,k,l are set to 1, 2, 3, and 4 (for different edges and faces
% the indices need to be changed.)

disp('Order 1')
homogeneous_order = 1;
edge_function = -L2*gL1+L1*gL2;
good_function = check_space(edge_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", edge function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", edge function is not in Nedelec space")
end

disp('Order 2')
homogeneous_order = 2;
edge_function = 1/14*((2*L2+9*L1*L2-9*L2^2)*gL1+(2*L1+9*L1*L2-9*L1^2)*gL2);
good_function = check_space(edge_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", edge function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", edge function is not in Nedelec space")
end
face_function_1 = L1*L3*gL2-L1*L2*gL3;
good_function = check_space(face_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", first face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", first face function is not in Nedelec space")
end
face_function_2 = -L2*L3*gL1+1/4*3*L1*L3*gL2+1/4*L1*L2*gL3;
good_function = check_space(face_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", second face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", second face function is not in Nedelec space")
end

disp('Order 3')
homogeneous_order = 3;
edge_function = 1/276*((-9*L2+3*L1*L2-78*L1^2*L2+41*L2^2+196*L1*L2^2-78*L2^3)*gL1+(9*L1-41*L1^2+78*L1^3-3*L1*L2-196*L1^2*L2+78*L1*L2^2)*gL2);
good_function = check_space(edge_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", edge function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", edge function is not in Nedelec space")
end
face_function_1 = 1/91*((-6*L2*L3+4*L1*L2*L3+3*L2*L3^2)*gL1+(37*L1*L3-82*L1^2*L3-L1*L3^2)*gL2+(-33*L1*L2+78*L1^2*L2-2*L1*L2*L3)*gL3);
good_function = check_space(face_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", first face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", first face function is not in Nedelec space")
end
face_function_2 = 1/580*(498*L2*L3-719*L1*L2*L3-539*L2*L3^2)*gL1+1/1259*(-266*L1*L3+1135*L1^2*L3+390*L1*L3^2)*gL2+1/870*(-24*L1*L2+294*L1^2*L2+539*L1*L2*L3)*gL3;
good_function = check_space(face_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", second face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", second face function is not in Nedelec space")
end
face_function_3 = 1/448*(-89*L2*L3-163*L1*L2*L3+213*L2*L3^2)*gL1+1/384*(-22*L1*L3+80*L1^2*L3-61*L1*L3^2)*gL2+1/192*(-12*L1*L2+30*L1^2*L2-61*L1*L2*L3)*gL3;
good_function = check_space(face_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", third face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", third face function is not in Nedelec space")
end
face_function_4 = 1/314*(207*L2*L3-1176*L1*L2*L3-441*L2*L3^2)*gL1+1/314*(-465*L1*L3+672*L1^2*L3+147*L1*L3^2)*gL2+1/157*(-195*L1*L2+252*L1^2*L2+147*L1*L2*L3)*gL3;
good_function = check_space(face_function_4, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", fourth face function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", fourth face function is not in Nedelec space")
end
vol_function_1 = ((3*L2*L3-2*L1*L2*L3-3*L2^2*L3-3*L2*L3^2)*gL1+(-1*L1*L3+1*L1^2*L3+2*L1*L2*L3+1*L1*L3^2)*gL2+(-1*L1*L2+1*L1^2*L2+1*L1*L2^2+2*L1*L2*L3)*gL3);
good_function = check_space(vol_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", first vol function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", first vol function is not in Nedelec space")
end
vol_function_2 = 1/15*((24*L2*L3-56*L1*L2*L3-24*L2^2*L3-24*L2*L3^2)*gL1+(-48*L1*L3+48*L1^2*L3+16*L1*L2*L3+48*L1*L3^2)*gL2+(-8*L1*L2+8*L1^2*L2+8*L1*L2^2-24*L1*L2*L3)*gL3);
good_function = check_space(vol_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", second vol function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", second vol function is not in Nedelec space")
end
vol_function_3 = 1/33*((-12*L2*L3+28*L1*L2*L3+12*L2^2*L3+12*L2*L3^2)*gL1+(-64*L1*L3+64*L1^2*L3+80*L1*L2*L3+64*L1*L3^2)*gL2+(92*L1*L2-92*L1^2*L2-92*L1*L2^2-76*L1*L2*L3)*gL3);
good_function = check_space(vol_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("For order "+num2str(homogeneous_order)+", third vol function is in Nedelec space")
else
    disp("For order "+num2str(homogeneous_order)+", third vol function is not in Nedelec space")
end