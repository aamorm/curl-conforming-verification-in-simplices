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

% In face 123
disp('Edge-based face functions')
for homogeneous_order = 2:maximum_order_checked
    for l = 0:homogeneous_order-2
        edge_function = whitney12;
        % Edge 12
        beta_e12 = L1*L2;
        xi_12 = L2-L1;
        psi_l_12 = get_legendre_poly(l, xi_12);
        basis_function_edge12 = beta_e12*psi_l_12*gL3;
        good_function = check_space(basis_function_edge12, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", first edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", first edge-based face function is not in Nedelec space")
        end

        % Edge 13
        beta_e13 = L1*L3;
        xi_13 = L3-L1;
        psi_l_13 = get_legendre_poly(l, xi_13);
        basis_function_edge13 = beta_e13*psi_l_13*gL2;
        good_function = check_space(basis_function_edge13, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", second edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", second edge-based face function is not in Nedelec space")
        end

        % Edge 23
        beta_e23 = L2*L3;
        xi_23 = L3-L2;
        psi_l_23 = get_legendre_poly(l, xi_23);
        basis_function_edge23 = beta_e23*psi_l_23*gL1;
        good_function = check_space(basis_function_edge23, homogeneous_order, affine_coordinates, tolerance, detailed_output);
        if (good_function)
            disp("For order "+num2str(homogeneous_order)+", third edge-based face function is in Nedelec space")
        else
            disp("For order "+num2str(homogeneous_order)+", third edge-based face function is not in Nedelec space")
        end
    end
end

disp('Face bubble functions')
% Only face 123
for homogeneous_order = 3:maximum_order_checked
    for l = 0:homogeneous_order-3
        for m = 0:homogeneous_order-3
            if ((l+m) <= (homogeneous_order-3))
                if (detailed_output)
                    disp("Checking l="+num2str(l)+", m="+num2str(m))
                end
                beta_face = L1*L2*L3;
                xi_12 = L2-L1;
                xi_13 = L3-L1;
                psi_l_12 = get_legendre_poly(l, xi_12);
                psi_m_13 = get_legendre_poly(m, xi_13);
                tau_12 = gL2-gL1;
                tau_13 = gL3-gL1;
                face_function_1 = beta_face*psi_l_12*psi_m_13*tau_12;
                good_function = check_space(face_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                if (good_function)
                    disp("For order "+num2str(homogeneous_order)+", first face function is in Nedelec space")
                else
                    disp("For order "+num2str(homogeneous_order)+", first face function is not in Nedelec space")
                end
                % I think that this should be psi_m_12 and psi_l_13, but...
                face_function_3 = beta_face*psi_l_12*psi_m_13*tau_13;
                good_function = check_space(face_function_3, homogeneous_order, affine_coordinates, tolerance, detailed_output);
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
% Only face 123
for homogeneous_order = 3:maximum_order_checked
    for l = 0:homogeneous_order-3
        for m = 0:homogeneous_order-3
            if ((l+m) <= (homogeneous_order-3))
                if (detailed_output)
                    disp("Checking l="+num2str(l)+", m="+num2str(m))
                end
                beta_face = L1*L2*L3;
                xi_12 = L2-L1;
                xi_13 = L3-L1;
                psi_l_12 = get_legendre_poly(l, xi_12);
                psi_m_13 = get_legendre_poly(m, xi_13);
                vol_function = beta_face*psi_l_12*psi_m_13*gL4;
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
                    beta_vol = L1*L2*L3*L4;
                    xi_12 = L2-L1;
                    xi_13 = L3-L1;
                    xi_14 = L4-L1;
                    psi_l_12 = get_legendre_poly(l, xi_12);
                    psi_m_13 = get_legendre_poly(m, xi_13);
                    psi_n_14 = get_legendre_poly(n, xi_14);
                    vol_function_1 = beta_face*psi_l_12*psi_m_13*psi_n_14*gL1;
                    good_function = check_space(vol_function_1, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", first bubble interior function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", first bubble interior function is not in Nedelec space")
                    end
                    vol_function_2 = beta_face*psi_l_12*psi_m_13*psi_n_14*gL2;
                    good_function = check_space(vol_function_2, homogeneous_order, affine_coordinates, tolerance, detailed_output);
                    if (good_function)
                        disp("For order "+num2str(homogeneous_order)+", second bubble interior function is in Nedelec space")
                    else
                        disp("For order "+num2str(homogeneous_order)+", second bubble interior function is not in Nedelec space")
                    end
                    vol_function_3 = beta_face*psi_l_12*psi_m_13*psi_n_14*gL3;
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

function legendre_poly = get_legendre_poly(i, x)
    if (i == 0)
        legendre_poly = 1;
    elseif (i == 1)
        legendre_poly = x;
    elseif (i>=2)
        legendre_poly = (1/i)*((2*i-1)*x*get_legendre_poly(i-1,x) - (i-1)*get_legendre_poly(i-2,x));
    end
end

% function gegenbauer_poly = get_gegenbauer_poly(i, alpha, x)
%     if (i == 0)
%         gegenbauer_poly = 1;
%     elseif (i == 1)
%         gegenbauer_poly = 2*alpha*x;
%     elseif (i>=2)
%         gegenbauer_poly = (1/i)*(2*(alpha+i-1)*x*get_gegenbauer_poly(i-1,alpha,x) - (2*alpha+i-2)*get_gegenbauer_poly(i-2,alpha,x));
%     end
% end
