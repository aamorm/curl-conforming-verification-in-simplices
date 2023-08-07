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

maximum_order_checked = 3;
% Higher orders are rational functions that the symbolic implementation of MATLAb does not
% cancel out with (1-L_{1,2,3,4})^i (although they should)

% The edge functions are Whitney functions (the test holds) and curl-free
% functions (they belong to the irrotational space).

% In face 123, so j1=4, j2=1, j3=2, j4=3.
disp('Edge-based face functions')
for homogeneous_order = 2:maximum_order_checked
    for i = 0:homogeneous_order-2
        c_i = (i+3)*sqrt((2*i+4)*(2*i+5)*(2*i+7)/(i+1));

        % Edge 12. k1 = 1, k2 = 2, k3 = 3.
        jacob_poly = get_jacobi_poly(i,1,2,(2*L2/(1-L1)-1));
        basis_function_edge12 = simplify(c_i*L1*L2*(1-L1)^i*jacob_poly*gL3);

        good_function = check_space(basis_function_edge12, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", first edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", first edge-based face function is not in Nedelec space")
        end

        % Edge 13. k1 = 1, k2 = 3, k3 = 2.
        jacob_poly = get_jacobi_poly(i,1,2,(2*L3/(1-L1)-1));
        basis_function_edge13 = simplify(c_i*L1*L3*(1-L1)^i*jacob_poly*gL2);
        good_function = check_space(basis_function_edge13, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", second edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", second edge-based face function is not in Nedelec space")
        end

        % Edge 23. k1 = 2, k2 = 3, k3 = 1.
        jacob_poly = get_jacobi_poly(i,1,2,(2*L3/(1-L2)-1));
        basis_function_edge23 = simplify(c_i*L2*L3*(1-L2)^i*jacob_poly*gL1);
        good_function = check_space(basis_function_edge23, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", third edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", third edge-based face function is not in Nedelec space")
        end
    end
end

disp('Face bubble functions')
% In face 123, so j1=4, j2=1, j3=2, j4=3.
for homogeneous_order = 3:maximum_order_checked
    for m = 0:homogeneous_order-3
        for n = 0:homogeneous_order-3
            if ((m+n) <= (homogeneous_order-3))
                if (detailed_output)
                    disp("Checking m="+num2str(m)+", n="+num2str(n))
                end
                c_mn1 = sqrt((2*n+3)*(m+n+3)*(m+2*n+4)*(m+2*n+5));
                c_mn2 = sqrt((2*m+2*n+7)*(2*m+2*n+8)*(2*m+2*n+9)/((m+1)*(m+2)));
                upsilon = c_mn1*c_mn2*L1*L2*L3;

                affine_factors = (1-L1)^m*(1-L1-L2)^n;

                jacobi_polynomial_m = get_jacobi_poly(m,2*n+3,2,(2*L2/(1-L1)-1));
                jacobi_polynomial_n = get_jacobi_poly(n,0,2,(2*L3/(1-L1-L2)-1));

                tau_12 = gL2-gL1;
                tau_13 = gL3-gL1;
                face_function_42 = upsilon*affine_factors*jacobi_polynomial_m*jacobi_polynomial_n*tau_12;
                good_function = check_space(face_function_42, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", first face function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", first face function is not in Nedelec space")
                end
                % I think that this should be psi_m_12 and psi_l_13, but...
                face_function_43 = upsilon*affine_factors*jacobi_polynomial_m*jacobi_polynomial_n*tau_13;
                good_function = check_space(face_function_43, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", second face function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", second face function is not in Nedelec space")
                end
            end
        end
    end
end

disp('Face-based interior functions')
% In face 123, so j1=4, j2=1, j3=2, j4=3.
for homogeneous_order = 3:maximum_order_checked
    for l = 0:homogeneous_order-3
        for m = 0:homogeneous_order-3
            if ((l+m) <= (homogeneous_order-3))
                if (detailed_output)
                    disp("Checking l="+num2str(l)+", m="+num2str(m))
                end
                upsilon = c_mn1*c_mn2*L1*L2*L3;
                affine_factors = (1-L1)^m*(1-L1-L2)^n;
                jacobi_polynomial_m = get_jacobi_poly(m,2*n+3,2,(2*L2/(1-L1)-1));
                jacobi_polynomial_n = get_jacobi_poly(n,0,2,(2*L3/(1-L1-L2)-1));

                vol_function = upsilon*affine_factors*jacobi_polynomial_m*jacobi_polynomial_n*gL4;
                good_function = check_space(vol_function, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", face-based interior function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", face-based interior function is not in Nedelec space")
                end
            end
        end
    end
end

disp('Interior bubble function')
for homogeneous_order = 4:maximum_order_checked
    for l = 0:homogeneous_order-4
        for m = 0:homogeneous_order-4
            for n = 0:homogeneous_order-4
                if ((l+m+n) <= (homogeneous_order-4))
                    if (detailed_output)
                        disp("Checking l="+num2str(l)+", m="+num2str(m)+", n="+num2str(n))
                    end
                    clmn1 = sqrt((l+2*m+2*n+9)*(l+2*m+2*n+10)*(2*l+2*m+2*n+11)*(m+2*n+6)/((l+1)*(m+1)*(n+1)));
                    clmn2 = sqrt((m+2*n+7)*(2*m+2*n+8)*(n+3)*(n+4)*(2*n+5)/((l+2)*(m+2)*(n+2)));
                    clmn = clmn1*clmn2;
                    % L0 here is L4.
                    psi_factor = clmn*L1*L2*L3*L4*(1-L1)^m*(1-L1-L2)^n;

                    jacobi_polynomial_l = get_jacobi_poly(l,2*m+2*n+8,2,(2*L1-1));
                    jacobi_polynomial_m = get_jacobi_poly(m,2*n+5,2,(2*L2/(1-L1)-1));
                    jacobi_polynomial_n = get_jacobi_poly(n,2,2,(2*L3/(1-L1-L2)-1));

                    vol_function_1 = psi_factor*jacobi_polynomial_l*jacobi_polynomial_m*jacobi_polynomial_n*gL1;
                    good_function = check_space(vol_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", first bubble interior function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", first bubble interior function is not in Nedelec space")
                    end
                    vol_function_2 = psi_factor*jacobi_polynomial_l*jacobi_polynomial_m*jacobi_polynomial_n*gL2;
                    good_function = check_space(vol_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", second bubble interior function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", second bubble interior function is not in Nedelec space")
                    end
                    vol_function_3 = psi_factor*jacobi_polynomial_l*jacobi_polynomial_m*jacobi_polynomial_n*gL3;
                    good_function = check_space(vol_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", third bubble interior function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", third bubble interior function is not in Nedelec space")
                    end
                end
            end
        end
    end
end

function [jacob_poly] = get_jacobi_poly (n, alpha, beta, x)

    if (n == 0)
        jacob_poly = 1;
    elseif (n==1)
        jacob_poly = (alpha+1) + (alpha+beta+2)*(x-1)/2;
    elseif (n>=2)
        jacob_poly = (2*n+alpha+beta-1)*((2*n+alpha+beta)*(2*n+alpha+beta-2)*x+alpha^2-beta^2)*get_jacobi_poly(n-1, alpha, beta, x) - ...
                      2*(n+alpha-1)*(n+beta-1)*(2*n+alpha+beta)*get_jacobi_poly(n-2, alpha, beta, x);
        jacob_poly = jacob_poly/(2*n*(n+alpha+beta)*(2*n+alpha+beta-2));
    end

end