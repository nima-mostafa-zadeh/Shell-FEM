function main(tt::Float64)

    cylinder_radius::Float64 = 1.0
    cylinder_height::Float64 = 1.0
    vesicle_elements_num::Int64 = 10
    vesicle_nodes_num::Int64 = vesicle_elements_num + 1
    vesicle_nodes::Matrix{Float64} = begin
        Y = range(cylinder_height, 0.0, length = vesicle_nodes_num)
        [cylinder_radius * ones(Float64, vesicle_nodes_num, 1) Y]
    end
    vesicle_elements::Matrix{Int} = begin
        [collect(1:(vesicle_nodes_num - 1)) collect(2:vesicle_nodes_num)]
    end

    material_density::Float64 = 10.0
    material_damping::Float64 = 100.0
    material_rotational_inertia_scalar::Float64 = 500.0

    material_shear_modulus::Float64 = 500.0
    material_bulk_modulus::Float64 = 10000.0
    material_shear_modulus_transverse::Float64 = 1.0

    material_poisson_ratio::Float64 = (3 * material_bulk_modulus - 2 * material_shear_modulus) /
            (6 * material_bulk_modulus + 2 * material_shear_modulus)
    material_youngs_modulus::Float64 = 2 * material_shear_modulus * (1 + material_poisson_ratio)
    material_factor = material_youngs_modulus / (1 - material_poisson_ratio^2)

    job_name::String = "Job_DOF_C4"
    
    section_thickness::Float64 = 0.01
    section_integration_points_num::Int64 = 4
    t::Vector{Float64} = [-1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), 0.0, 0.0]
    s::Vector{Float64} = [-1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3)]
    w::Vector{Float64} = [1.0, 1.0, 1.0, 1.0, 2.0, 2.0]

    step_dt::Float64 = 5.0e-7
    step_total_time::Float64 = tt
    step_steps_num::Int64 = floor(step_total_time / step_dt)
    step_print::Int64 = 1000
    step_time::Float64 = 0.0

    force_external_magnitude::Float64 = 0.01
    step_load = (0.0, 1.0, 10.0)

    coordinates_initial::Vector{Float64} = reshape(vesicle_nodes, 2 * vesicle_nodes_num)
    velocity_initial::Vector{Float64} = zeros(2 * vesicle_nodes_num)
    acceleration_initial::Vector{Float64} = zeros(2 * vesicle_nodes_num)
    directors_initial::Vector{Float64} = begin
        [ones(Float64, vesicle_nodes_num); zeros(Float64, vesicle_nodes_num)]
    end
    angular_velocity_initial::Vector{Float64} = zeros(vesicle_nodes_num)
    angular_acceleration_initial::Vector{Float64} = zeros(vesicle_nodes_num)

    coordinates_current::Vector{Float64} = reshape(vesicle_nodes, 2 * vesicle_nodes_num)
    velocity_current::Vector{Float64} = zeros(2 * vesicle_nodes_num)
    acceleration_current::Vector{Float64} = zeros(2 * vesicle_nodes_num)
    directors_current::Vector{Float64} = begin
        [ones(Float64, vesicle_nodes_num); zeros(Float64, vesicle_nodes_num)]
    end
    angular_velocity_current::Vector{Float64} = zeros(vesicle_nodes_num)
    angular_acceleration_current::Vector{Float64} = zeros(vesicle_nodes_num)

    boundary_condition = [
        (type = "TRANSLATION FIX",
        index = 2 * vesicle_nodes_num,
        value = coordinates_initial[2 * vesicle_nodes_num]),
        (type = "ROTATION FIX",
        index = vesicle_nodes_num,
        value = directors_initial[vesicle_nodes_num]),
        (type = "ROTATION FIX",
        index = 2 * vesicle_nodes_num,
        value = directors_initial[2 * vesicle_nodes_num]),
    ]

    mass_matrix_global::Vector{Float64} = zeros(2 * vesicle_nodes_num)
    rotational_inertia_global::Vector{Float64} = zeros(vesicle_nodes_num)

    mass_matrix_local::Matrix{Float64} = zeros(vesicle_elements_num, 4)
    rotational_inertia_local::Matrix{Float64} = zeros(vesicle_elements_num, 4)
    rotational_inertia_local_compacted::Matrix{Float64} = zeros(vesicle_elements_num, 2)

    for i = 1:vesicle_elements_num
        x0_1 = coordinates_initial[i]
        x0_2 = coordinates_initial[i + 1]
        y0_1 = coordinates_initial[vesicle_nodes_num + i]
        y0_2 = coordinates_initial[vesicle_nodes_num + i + 1]

        px0_1 = directors_initial[i]
        px0_2 = directors_initial[i + 1]
        py0_1 = directors_initial[vesicle_nodes_num + i]
        py0_2 = directors_initial[vesicle_nodes_num + i + 1]

        for ip = 1:section_integration_points_num
            n1 = 0.5 * (1.0 - t[ip])
            n2 = 0.5 * (1.0 + t[ip])
            n3 = 0.5 * (1.0 - t[ip]) * section_thickness * s[ip] / 2.0
            n4 = 0.5 * (1.0 + t[ip]) * section_thickness * s[ip] / 2.0

            dn1_dt = -0.5
            dn2_dt = 0.5
            dn3_dt = -0.5 * section_thickness * s[ip] / 2.0
            dn4_dt = 0.5 * section_thickness * s[ip] / 2.0

            dn1_ds = 0.0
            dn2_ds = 0.0
            dn3_ds = 0.5 * (1.0 - t[ip]) * section_thickness / 2.0
            dn4_ds = 0.5 * (1.0 + t[ip]) * section_thickness / 2.0

            x_ip = n1 * x0_1 + 
                   n2 * x0_2 + 
                   n3 * px0_1 + 
                   n4 * px0_2

            dx_dt = dn1_dt * x0_1 + 
                    dn2_dt * x0_2 + 
                    dn3_dt * px0_1 + 
                    dn4_dt * px0_2
            
            dy_dt = dn1_dt * y0_1 + 
                    dn2_dt * y0_2 + 
                    dn3_dt * py0_1 + 
                    dn4_dt * py0_2

            dx_ds = dn1_ds * x0_1 + 
                    dn2_ds * x0_2 + 
                    dn3_ds * px0_1 + 
                    dn4_ds * px0_2

            dy_ds = dn1_ds * y0_1 + 
                    dn2_ds * y0_2 + 
                    dn3_ds * py0_1 + 
                    dn4_ds * py0_2

            det_j = dx_dt * dy_ds - dx_ds * dy_dt

            mass_factor = (2 * pi * material_density * w[ip] * det_j * x_ip)
            mass_matrix_local[i, 1] += mass_factor * n1
            mass_matrix_local[i, 2] += mass_factor * n2
            mass_matrix_local[i, 3] += mass_factor * n1
            mass_matrix_local[i, 4] += mass_factor * n2

            rotational_inertia_local[i, 1] += mass_factor * n3
            rotational_inertia_local[i, 2] += mass_factor * n4
            rotational_inertia_local[i, 3] += mass_factor * n3
            rotational_inertia_local[i, 4] += mass_factor * n4
        end

        rotational_inertia_local_compacted[i, 1] = 
            rotational_inertia_local[i, 1] * py0_1^2 + 
            rotational_inertia_local[i, 3] * px0_1^2
        rotational_inertia_local_compacted[i, 2] = 
            rotational_inertia_local[i, 2] * py0_2^2 + 
            rotational_inertia_local[i, 4] * px0_2^2        
    end
    
    for i = 1:vesicle_elements_num
        mass_matrix_global[i] += mass_matrix_local[i, 1]
        mass_matrix_global[i + 1] += mass_matrix_local[i, 2]
        mass_matrix_global[i + vesicle_nodes_num] += mass_matrix_local[i, 3]
        mass_matrix_global[i + vesicle_nodes_num + 1] += mass_matrix_local[i, 4]
        rotational_inertia_global[i] += rotational_inertia_local_compacted[i, 1]
        rotational_inertia_global[i + 1] += rotational_inertia_local_compacted[i, 2]
    end
    
    rotational_inertia_global *= material_rotational_inertia_scalar

    force_internal_elements = Matrix{Float64}(undef, vesicle_elements_num, 4)
    rotational_force_internal_elements = Matrix{Float64}(undef, vesicle_elements_num, 4)
    torque_internal_elements = Matrix{Float64}(undef, vesicle_elements_num, 2)

    force_internal_global = Vector{Float64}(undef, 2 * vesicle_nodes_num)
    force_external_global = Vector{Float64}(undef, 2 * vesicle_nodes_num)
    torque_global = Vector{Float64}(undef, vesicle_nodes_num)

    b_matrix = Matrix{Float64}(undef, 5, 8)

    open(job_name * ".csv", "w") do io
        println(io, "t,x,y")  # Write header
    end
    
    
    for step = 1:step_steps_num
        fill!(force_internal_elements, 0.0)
        fill!(rotational_force_internal_elements, 0.0)
        fill!(torque_internal_elements, 0.0)
        fill!(force_internal_global, 0.0)
        fill!(force_external_global, 0.0)
        fill!(torque_global, 0.0)

        for i = 1:vesicle_elements_num

            x_1 = coordinates_current[i]
            x_2 = coordinates_current[i + 1]
            y_1 = coordinates_current[vesicle_nodes_num + i]
            y_2 = coordinates_current[vesicle_nodes_num + i + 1]

            px_1 = directors_current[i]
            px_2 = directors_current[i + 1]
            py_1 = directors_current[vesicle_nodes_num + i]
            py_2 = directors_current[vesicle_nodes_num + i + 1]

            x0_1 = coordinates_initial[i]
            x0_2 = coordinates_initial[i + 1]
            y0_1 = coordinates_initial[vesicle_nodes_num + i]
            y0_2 = coordinates_initial[vesicle_nodes_num + i + 1]

            px0_1 = directors_initial[i]
            px0_2 = directors_initial[i + 1]
            py0_1 = directors_initial[vesicle_nodes_num + i]
            py0_2 = directors_initial[vesicle_nodes_num + i + 1]

            for ip = 1:section_integration_points_num
                n1 = 0.5 * (1.0 - t[ip])
                n2 = 0.5 * (1.0 + t[ip])
                n3 = 0.5 * (1.0 - t[ip]) * section_thickness * s[ip] / 2.0
                n4 = 0.5 * (1.0 + t[ip]) * section_thickness * s[ip] / 2.0

                dn1_dt = -0.5
                dn2_dt = 0.5
                dn3_dt = -0.5 * section_thickness * s[ip] / 2.0
                dn4_dt = 0.5 * section_thickness * s[ip] / 2.0

                dn1_ds = 0.0
                dn2_ds = 0.0
                dn3_ds = 0.5 * (1.0 - t[ip]) * section_thickness / 2.0
                dn4_ds = 0.5 * (1.0 + t[ip]) * section_thickness / 2.0

                x_ip = n1 * x0_1 + 
                       n2 * x0_2 + 
                       n3 * px0_1 + 
                       n4 * px0_2

                dx_dt = dn1_dt * x0_1 + 
                        dn2_dt * x0_2 + 
                        dn3_dt * px0_1 + 
                        dn4_dt * px0_2
                
                dy_dt = dn1_dt * y0_1 + 
                        dn2_dt * y0_2 + 
                        dn3_dt * py0_1 + 
                        dn4_dt * py0_2

                #println("$(x_ip), $(dx_dt), $(dy_dt)")

                q11 = dx_dt / sqrt(dx_dt^2 + dy_dt^2)
                q12 = dy_dt / sqrt(dx_dt^2 + dy_dt^2)
                q21 = -q12
                q22 = q11

                x_1_rot = q11 * x_1 + q12 * y_1
                y_1_rot = q21 * x_1 + q22 * y_1
                
                x_2_rot = q11 * x_2 + q12 * y_2
                y_2_rot = q21 * x_2 + q22 * y_2

                x0_1_rot = q11 * x0_1 + q12 * y0_1
                y0_1_rot = q21 * x0_1 + q22 * y0_1
                
                x0_2_rot = q11 * x0_2 + q12 * y0_2
                y0_2_rot = q21 * x0_2 + q22 * y0_2

                px_1_rot = q11 * px_1 + q12 * py_1
                py_1_rot = q21 * px_1 + q22 * py_1
                
                px_2_rot = q11 * px_2 + q12 * py_2
                py_2_rot = q21 * px_2 + q22 * py_2

                px0_1_rot = q11 * px0_1 + q12 * py0_1
                py0_1_rot = q21 * px0_1 + q22 * py0_1
                
                px0_2_rot = q11 * px0_2 + q12 * py0_2
                py0_2_rot = q21 * px0_2 + q22 * py0_2

                dx_dt = dn1_dt * x0_1_rot + 
                        dn2_dt * x0_2_rot + 
                        dn3_dt * px0_1_rot + 
                        dn4_dt * px0_2_rot
                
                dy_dt = dn1_dt * y0_1_rot + 
                        dn2_dt * y0_2_rot + 
                        dn3_dt * py0_1_rot + 
                        dn4_dt * py0_2_rot

                dx_ds = dn1_ds * x0_1_rot + 
                        dn2_ds * x0_2_rot + 
                        dn3_ds * px0_1_rot + 
                        dn4_ds * px0_2_rot

                dy_ds = dn1_ds * y0_1_rot + 
                        dn2_ds * y0_2_rot + 
                        dn3_ds * py0_1_rot + 
                        dn4_ds * py0_2_rot

                #println("$(dx_dt), $(dy_dt), $(dx_ds), $(dy_ds)")

                det_j = dx_dt * dy_ds - dx_ds * dy_dt

                dt_dx = dy_ds / det_j
                dt_dy = -dx_ds / det_j
                ds_dx = -dy_dt / det_j
                ds_dy = dx_dt / det_j

                dn1_dx = dt_dx * dn1_dt + ds_dx * dn1_ds
                dn2_dx = dt_dx * dn2_dt + ds_dx * dn2_ds
                dn3_dx = dt_dx * dn3_dt + ds_dx * dn3_ds
                dn4_dx = dt_dx * dn4_dt + ds_dx * dn4_ds

                dn1_dy = dt_dy * dn1_dt + ds_dy * dn1_ds
                dn2_dy = dt_dy * dn2_dt + ds_dy * dn2_ds
                dn3_dy = dt_dy * dn3_dt + ds_dy * dn3_ds
                dn4_dy = dt_dy * dn4_dt + ds_dy * dn4_ds

                b11 = dn1_dx; b12 = dn2_dx; b15 = dn3_dx; b16 = dn4_dx;
                b21 = dn1_dy; b22 = dn2_dy; b25 = dn3_dy; b26 = dn4_dy;
                b33 = dn1_dx; b34 = dn2_dx; b37 = dn3_dx; b38 = dn4_dx;
                b43 = dn1_dy; b44 = dn2_dy; b47 = dn3_dy; b48 = dn4_dy;

                b51 = q11 * n1 / x_ip
                b52 = q11 * n2 / x_ip
                b53 = q21 * n1 / x_ip
                b54 = q21 * n2 / x_ip
                b55 = q11 * n3 / x_ip
                b56 = q11 * n4 / x_ip
                b57 = q21 * n3 / x_ip
                b58 = q21 * n4 / x_ip

                f_inv_11 = b11 * x0_1_rot + b12 * x0_2_rot + b15 * px0_1_rot + b16 * px0_2_rot
                f_inv_12 = b21 * x0_1_rot + b22 * x0_2_rot + b25 * px0_1_rot + b26 * px0_2_rot
                f_inv_21 = b33 * y0_1_rot + b34 * y0_2_rot + b37 * py0_1_rot + b38 * py0_2_rot
                f_inv_22 = b43 * y0_1_rot + b44 * y0_2_rot + b47 * py0_1_rot + b48 * py0_2_rot 
                f_inv_33 = b51 * x0_1_rot +
                           b52 * x0_2_rot +
                           b53 * y0_1_rot +
                           b54 * y0_2_rot + 
                           b55 * px0_1_rot + 
                           b56 * px0_2_rot +
                           b57 * py0_1_rot + 
                           b58 * py0_2_rot

                det_f_inv = f_inv_11 * f_inv_22 * f_inv_33 - 
                            f_inv_12 * f_inv_21 * f_inv_33
                
                f_grad_11 = f_inv_22 * f_inv_33 / det_f_inv;
                f_grad_12 = -f_inv_12 * f_inv_33 / det_f_inv;
                f_grad_21 = -f_inv_21 * f_inv_33 / det_f_inv;
                f_grad_22 = f_inv_11 * f_inv_33 / det_f_inv;
                f_grad_33 = 1 / f_inv_33;

                #println("$(f_grad_11), $(f_grad_12), $(f_grad_21), $(f_grad_22), $(f_grad_33)")

                epsilon_11 = f_grad_11 - 1.0
                epsilon_33 = f_grad_33 - 1.0

                if ip < 5
                    sigma_11 = material_factor * (epsilon_11 + material_poisson_ratio * epsilon_33)
                    sigma_12 = 0;
                    sigma_21 = 0;
                    sigma_22 = 0;
                    sigma_33 = material_factor * (epsilon_33 + material_poisson_ratio * epsilon_11)
                end
                
                force_factor = det_j * x_ip * 2 * pi * w[ip]
                force_local_rot_1 = (sigma_11 * b11 + sigma_12 * b21 + sigma_33 * b51) * force_factor
                force_local_rot_2 = (sigma_11 * b12 + sigma_12 * b22 + sigma_33 * b52) * force_factor
                force_local_rot_3 = (sigma_21 * b33 + sigma_22 * b43 + sigma_33 * b53) * force_factor
                force_local_rot_4 = (sigma_21 * b34 + sigma_22 * b44 + sigma_33 * b54) * force_factor
                force_local_rot_5 = (sigma_11 * b15 + sigma_12 * b25 + sigma_33 * b55) * force_factor
                force_local_rot_6 = (sigma_11 * b16 + sigma_12 * b26 + sigma_33 * b56) * force_factor
                force_local_rot_7 = (sigma_21 * b37 + sigma_22 * b47 + sigma_33 * b57) * force_factor
                force_local_rot_8 = (sigma_21 * b38 + sigma_22 * b48 + sigma_33 * b58) * force_factor
                
                force_local_1 = q11 * force_local_rot_1 + q21 * force_local_rot_3
                force_local_3 = q12 * force_local_rot_1 + q22 * force_local_rot_3
                force_local_2 = q11 * force_local_rot_2 + q21 * force_local_rot_4
                force_local_4 = q12 * force_local_rot_2 + q22 * force_local_rot_4
                force_local_5 = q11 * force_local_rot_5 + q21 * force_local_rot_7
                force_local_7 = q12 * force_local_rot_5 + q22 * force_local_rot_7
                force_local_6 = q11 * force_local_rot_6 + q21 * force_local_rot_8
                force_local_8 = q12 * force_local_rot_6 + q22 * force_local_rot_8

                
                force_internal_elements[i, 1] += force_local_1
                force_internal_elements[i, 2] += force_local_2
                force_internal_elements[i, 3] += force_local_3
                force_internal_elements[i, 4] += force_local_4

                rotational_force_internal_elements[i, 1] += force_local_5
                rotational_force_internal_elements[i, 2] += force_local_6
                rotational_force_internal_elements[i, 3] += force_local_7
                rotational_force_internal_elements[i, 4] += force_local_8
                
            end

            torque_internal_elements[i, 1] = -rotational_force_internal_elements[i, 1] * py_1 + 
                            rotational_force_internal_elements[i, 3] * px_1;
            torque_internal_elements[i, 2] = -rotational_force_internal_elements[i, 2] * py_2 +
                            rotational_force_internal_elements[i, 4] * px_2;
        end
        for i = 1:vesicle_elements_num
            force_internal_global[i] += force_internal_elements[i, 1]
            force_internal_global[i + 1] += force_internal_elements[i, 2]
            force_internal_global[i + vesicle_nodes_num] += force_internal_elements[i, 3]
            force_internal_global[i + vesicle_nodes_num + 1] += force_internal_elements[i, 4]
            torque_global[i] += torque_internal_elements[i, 1]
            torque_global[i + 1] += torque_internal_elements[i, 2]
        end
        if step_time >= step_load[1] && step_time < step_load[2]
            force_external_global[vesicle_nodes_num + 1] = force_external_magnitude * step_time / step_load[2]
        elseif step_time >= step_load[2] && step_time < step_load[3]
            force_external_global[vesicle_nodes_num + 1] = force_external_magnitude
        else
            force_external_global[vesicle_nodes_num + 1] = 0.0
        end

        for i = 1:(2 * vesicle_nodes_num)
            force_internal_global[i] =- material_damping * mass_matrix_global[i] * velocity_current[i]
            acceleration_previous = acceleration_current[i]
            acceleration_current[i] = (force_external_global[i] - force_internal_global[i]) / mass_matrix_global[i]
            velocity_current[i] += 0.5 * step_dt * (acceleration_current[i] + acceleration_previous);
            coordinates_current[i] += velocity_current[i] * step_dt + 0.5 * acceleration_previous * step_dt^2
        end
        
        for i = 1:vesicle_nodes_num
            torque_global[i] -= material_damping * rotational_inertia_global[i] * angular_velocity_current[i]
            angular_acceleration_previous = angular_acceleration_current[i]
            angular_acceleration_current[i] = torque_global[i] / rotational_inertia_global[i]
            angular_velocity_current[i] += 0.5 * step_dt * (angular_acceleration_current[i] + angular_acceleration_previous)
            rotation_angle = angular_velocity_current[i] * step_dt + 0.5 * angular_acceleration_previous * step_dt^2
            
            px = directors_current[i];
            py = directors_current[i + vesicle_nodes_num]
            directors_current[i] = cos(rotation_angle) * px - sin(rotation_angle) * py
            directors_current[i + vesicle_nodes_num] = sin(rotation_angle) * px + cos(rotation_angle) * py
            director_normalization_factor = sqrt(directors_current[i]^2 + directors_current[i + vesicle_nodes_num]^2)
            directors_current[i] /= director_normalization_factor
            directors_current[i + vesicle_nodes_num] /= director_normalization_factor
        end
    
        for i = 1:length(boundary_condition)
            if boundary_condition[i].type == "TRANSLATION FIX"
                coordinates_current[boundary_condition[i].index] = boundary_condition[i].value;
            elseif boundary_condition[i].type == "ROTATION FIX"
                directors_current[boundary_condition[i].index] = boundary_condition[i].value;
            end
        end
        if mod(step, step_print) == 0
            open(job_name * ".csv", "a") do io
            println(io, step_time, ",", 
                        coordinates_current[1], ",", 
                        coordinates_current[vesicle_nodes_num + 1])
            end
            println("$(step_time),\t$(coordinates_current[1]),\t$(coordinates_current[vesicle_nodes_num + 1])")
        end

        step_time += step_dt
    end

end

main(5.0E-3)