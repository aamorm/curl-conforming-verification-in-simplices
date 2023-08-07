% See notes in 07.04.2019
% The basis functions are full order. For simplices, the only
% different is in the gradient (this means that the J operator
% only affects to the rotational part and it can be used
% for mixed-order and full-order functions)
% Here, the gradient of V_{p+1} is discarded.
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
whitney12 = L1*gL2-L2*gL1;
tolerance = 1e-10;
detailed_output = true;

maximum_order_checked = 2;

% Edge functions are discarded since it is composed of gradient of p+1 order.
% The two families to be considered are phi_{1j}^{F,1} and phi_{ij}^{F,1}

disp('Face functions 1j')
for homogeneous_order = 2:maximum_order_checked
    i = 1;
    for j = 1:homogeneous_order-1
        if (detailed_output)
            disp("Checking i="+num2str(i)+", j="+num2str(j))
        end
        v_face = get_v_ij_face(i, j, L1, L2, L3);
        face_function = v_face*whitney12;
        good_function = check_space(face_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", face 1j function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", face 1j function is not in Nedelec space")
        end
    end
end

disp('Face functions ij')
for homogeneous_order = 2:maximum_order_checked
    for i = 2:homogeneous_order
        for j = 1:homogeneous_order
            if ((i+j)<=(homogeneous_order+1))
                if (detailed_output)
                    disp("Checking i="+num2str(i)+", j="+num2str(j))
                end
                u_i = get_u_i (i, L1, L2);
                v_ij = get_v_ij_face (i, j, L1, L2, L3);
                face_function = gradient(u_i, affine_coordinates)*v_ij - gradient(v_ij, affine_coordinates)*u_i;
                good_function = check_space(face_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", face ij function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", face ij function is not in Nedelec space")
                end
            end
        end
    end
end

% Non curl-free volume function (rotational part)
disp('Volume function 1jk (a)')
for homogeneous_order = 3:maximum_order_checked
    i = 1;
    for j = 1:homogeneous_order
        for k = 1:homogeneous_order
            if ((j+k)<=(homogeneous_order-1))
                if (detailed_output)
                    disp("Checking i="+num2str(i)+", j="+num2str(j)+", k="+num2str(k))
                end
                v_1j = get_v_ij_vol (i, j, affine_coordinates);
                w_1jk = get_w_ijk_vol (i, j, k, affine_coordinates);
                vol_function = whitney12*v_1j*w_1jk;
                good_function = check_space(vol_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", volume 1jk(a) function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", volume 1jk(a) function is not in Nedelec space")
                end
            end
        end
    end
end
disp('Volume function ijk (b)')
for homogeneous_order = 3:maximum_order_checked
    for i = 2:homogeneous_order
        for j = 1:homogeneous_order
            for k = 1:homogeneous_order
                if ((i+j+k)<=(homogeneous_order+1))
                    if (detailed_output)
                        disp("Checking i="+num2str(i)+", j="+num2str(j)+", k="+num2str(k))
                    end
                    u_i = get_u_i (i, L1, L2);
                    v_ij = get_v_ij_vol (i, j, affine_coordinates);
                    w_ijk = get_w_ijk_vol (i, j, k, affine_coordinates);
                    vol_function = v_ij*w_ijk*gradient(u_i,affine_coordinates);
                    good_function = check_space(vol_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", volume ijk(b) function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", volume ijk(b) function is not in Nedelec space")
                    end
                end
            end
        end
    end
end
disp('Volume function ijk (c)')
for homogeneous_order = 3:maximum_order_checked
    for i = 2:homogeneous_order
        for j = 1:homogeneous_order
            for k = 1:homogeneous_order
                if ((i+j+k)<=(homogeneous_order+1))
                    if (detailed_output)
                        disp("Checking i="+num2str(i)+", j="+num2str(j)+", k="+num2str(k))
                    end
                    u_i = get_u_i (i, L1, L2);
                    v_ij = get_v_ij_vol (i, j, affine_coordinates);
                    w_ijk = get_w_ijk_vol (i, j, k, affine_coordinates);
                    vol_function = u_i*v_ij*gradient(w_ijk,affine_coordinates);
                    good_function = check_space(vol_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", volume ijk(c) function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", volume ijk(c) function is not in Nedelec space")
                    end
                end
            end
        end
    end
end


function w_vol = get_w_ijk_vol (i, j, k, phi)

    w_vol = get_int_jacobi(2*phi(4)-1,k,2*i+2*j-2);

end

function v_vol = get_v_ij_vol (i, j, phi)

    v_vol = get_int_jacobi((2*phi(3)-(1-phi(4)))/(1-phi(4)),j,2*i-1)*(1-phi(4))^j;

end

function u = get_u_i (i, x1, x2)

    u = get_int_jacobi((x2-x1)/(x2+x1),i,0)*(x2+x1)^i;

end

function v_face = get_v_ij_face (i, j, x1, x2, x3)

    v_face = get_int_jacobi((x3-x2-x1),j,2*i-1);

end

function int_jacobi = get_int_jacobi (x, n, alpha)

    syms y
    assume(y,'real');

    sym_jacobi = get_sym_jacobi_beta_0 (y, n-1, alpha);
    int_jacobi = int(sym_jacobi, y, -1, x);

end

% x is symbolic variable
function sym_jacobi = get_sym_jacobi_beta_0 (x, n, alpha)

    constant = (-1)^n/(2^n*factorial(n)*(1-x)^alpha);

    diff_part = diff((1-x)^(n+alpha)*(1+x)^n,n);

    sym_jacobi = constant*diff_part;

end