% See notes in 21.03.2019
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
detailed_output = false;

maximum_order_checked = 4;

disp('Edge functions')
for homogeneous_order = 1:maximum_order_checked
    edge_function = get_hom_int_legendre_poly(homogeneous_order-1, L1, L2)*whitney12;
    good_function = check_space(edge_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
    if (good_function)
        disp("For order "+num2str(homogeneous_order)+", edge function is in Nedelec space")
    else
        disp("For order "+num2str(homogeneous_order)+", edge function is not in Nedelec space")
    end
end

disp('Face functions')
% Only one family (the other one is a permutation)
for homogeneous_order = 2:maximum_order_checked
    for i = 0:homogeneous_order
        for j = 1:homogeneous_order
            if ((i+j) <= (homogeneous_order-1))
                if (detailed_output)
                    disp("Checking i="+num2str(i)+", j="+num2str(j))
                end
                E_i = get_hom_int_legendre_poly(i, L1, L2)*whitney12;
                face_function = get_hom_int_jacobi_poly(j, 2*i+1, L1+L2, L3)*E_i;
                good_function = check_space(face_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", face function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", face function is not in Nedelec space")
                end
            end
        end
    end
end

disp('Volume functions')
% Only one family (the other two are permutations)
for homogeneous_order = 3:maximum_order_checked
    for i = 0:homogeneous_order
        for j = 1:homogeneous_order
            for k = 1:homogeneous_order
                if ((i+j+k) <= (homogeneous_order-1))
                    if (detailed_output)
                        disp("Checking i="+num2str(i)+", j="+num2str(j)+", k="+num2str(k))
                    end
                    E_i = get_hom_int_legendre_poly(i, L1, L2)*whitney12;
                    face_function = get_hom_int_jacobi_poly(j, 2*i+1, L1+L2, L3)*E_i;
                    vol_function = get_hom_int_jacobi_poly(k, 2*(i+j), 1-L4, L4)*face_function;
                    good_function = check_space(vol_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", vol function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", vol function is not in Nedelec space")
                    end
                end
            end
        end
    end
end

function [hom_int_legendre_poly] = get_hom_int_legendre_poly (i, s0, s1)
    hom_int_legendre_poly = get_int_legendre_poly (i, s1, s0+s1);
end

function int_legendre_poly = get_int_legendre_poly(i, x, t)
    if (i == 0)
        int_legendre_poly = 1;
    elseif (i == 1)
        int_legendre_poly = 2*x-t;
    elseif (i>=2)
        int_legendre_poly = (1/i)*((2*i-1)*(2*x-t)*get_int_legendre_poly(i-1,x,t) - (i-1)*t^2*get_int_legendre_poly(i-2,x,t));
    end
end

function [hom_int_jacob_poly] = get_hom_int_jacobi_poly (j, alpha, t0, t1)
    hom_int_jacob_poly = get_int_jacobi_poly (j, alpha, t1, t0+t1);
end

function [int_jacob_poly] = get_int_jacobi_poly (j, alpha, x, t)

    [a_j, b_j, c_j] = get_constants_for_int_jacobi (j, alpha);

    if (j == 1)
        int_jacob_poly = x;
    elseif (j>=2)
        int_jacob_poly = a_j*get_jacobi_poly (j, alpha, x, t) + ...
                         b_j*t*get_jacobi_poly (j-1, alpha, x, t) - ...
                         c_j*t^2*get_jacobi_poly (j-2, alpha, x, t);
    end

end

function [jacob_poly] = get_jacobi_poly (j, alpha, x, t)

    [a_j, b_j, c_j, d_j] = get_constants_for_jacobi (j, alpha);

    if (j == 0)
        jacob_poly = 1;
    elseif (j==1)
        jacob_poly = 2*x-t+alpha*x;
    elseif (j>=2)
        jacob_poly = (1/a_j)*b_j*(c_j*(2*x-t)+alpha^2*t)*get_jacobi_poly(j-1, alpha, x, t) - ...
                        d_j*t^2*get_jacobi_poly(j-2, alpha, x, t);
    end

end

function [a_j, b_j, c_j] = get_constants_for_int_jacobi (j, alpha)

    a_j = (j+alpha)/((2*j+alpha-1)*(2*j+alpha));
    b_j = alpha/((2*j+alpha-2)*(2*j+alpha));
    c_j = (j-1)/((2*j+alpha-2)*(2*j+alpha-1));

end

function [a_j, b_j, c_j, d_j] = get_constants_for_jacobi (j, alpha)

    a_j = 2*j*(j+alpha)*(2*j+alpha-2);
    b_j = 2*j+alpha-1;
    c_j = (2*j+alpha)*(2*j+alpha-2);
    d_j = 2*(j+alpha-1)*(j-1)*(2*j+alpha);

end