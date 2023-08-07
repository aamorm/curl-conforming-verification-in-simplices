function good_function = check_space (basis_function, homogeneous_order, sym_coordinates, tolerance, detailed_output)

    if (nargin == 3)
        tolerance = 1e-10;
        detailed_output = false;
    elseif (nargin == 4)
        detailed_output = false;
    end

    cartesian_dimension = length(basis_function);
    % Generate possible multi_indices
    possible_multi_indices = generate_all_possible_multi_indices (homogeneous_order, cartesian_dimension);

    coefficients_multi_indices = zeros(length(possible_multi_indices),1);
    for index_component = 1:cartesian_dimension
        [coefficients, exponents] = coeffs(basis_function(index_component), sym_coordinates);
        if (detailed_output)
            if (~isempty(coefficients))
                disp("For component #"+num2str(index_component));
            end
        end
        for index_monomial = 1:length(coefficients)
            % Working from R2018a
            polynomial_degrees = zeros(1,4);
            for index_coord = 1:length(polynomial_degrees)
                polynomial_degrees(index_coord) = polynomialDegree(exponents(index_monomial),sym_coordinates(index_coord));
            end
            multi_index_in_function = polynomial_degrees;
            multi_index_in_function(index_component) = multi_index_in_function(index_component) + 1;
            [is_within_possible, location_multi_index] = ismember(multi_index_in_function, possible_multi_indices, 'rows');
            if (~is_within_possible)
                if (detailed_output)
                    disp("Monomial "+char(exponents(index_monomial))+", with multi-index ["+num2str(multi_index_in_function)+"] not found")
                end
                if (sum(multi_index_in_function) < (homogeneous_order+1))
                    if (detailed_output)
                        disp("This monomial is in a lower-order space")
                    end
                else
                    if (detailed_output)
                        disp("This monomial is in a higher-order space")
                    end
                    good_function = false;
                    return
                end
            else
                coefficients_multi_indices(location_multi_index) = coefficients_multi_indices(location_multi_index) + double(coefficients(index_monomial));
                if (detailed_output)
                    disp("For monomial "+char(exponents(index_monomial))+"...")
                    disp("the multi-index is ["+num2str(multi_index_in_function)+"],")
                    disp("and its coefficient is "+num2str(double(coefficients(index_monomial)))+".")
                end
            end
        end
    end

    good_function = true;
    for index_multi_indices = 1:length(possible_multi_indices)
        if (coefficients_multi_indices(index_multi_indices) ~= 0)
            if (detailed_output)
                disp("The multi-index ["+num2str(possible_multi_indices(index_multi_indices,:))+"] is not zero, "+num2str(coefficients_multi_indices(index_multi_indices))+".");
            end
            if (abs(coefficients_multi_indices(index_multi_indices)) > tolerance)
                good_function = false;
            end
        end
    end
end

function multi_indices = generate_all_possible_multi_indices (homogeneous_order, cartesian_dimension)

    % Some of them are repeated.
    num_possible_multi_indices = nchoosek(cartesian_dimension+homogeneous_order, cartesian_dimension-1);
    multi_indices = zeros(num_possible_multi_indices, cartesian_dimension);
    counter_multi_index = 0;
    % Since the maximum value in the multi-index is the homogeneous_order,
    % we can convert
    num_max_possibilities = (homogeneous_order+1)^cartesian_dimension;
    for possible_counter = 1:num_max_possibilities
        possible_multi_index = dec2base(possible_counter,homogeneous_order+1,cartesian_dimension) - '0';
        if (sum(possible_multi_index) == (homogeneous_order+1))
            multi_indices(counter_multi_index+1,:) = possible_multi_index;
            counter_multi_index = counter_multi_index + 1;
        end
    end
    multi_indices = multi_indices(1:counter_multi_index,:);

end