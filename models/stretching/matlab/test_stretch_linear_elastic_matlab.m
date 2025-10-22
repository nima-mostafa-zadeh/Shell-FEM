clear;
clc;


%% PART
CYLINDER_RADIUS = 1;
CYLINDER_HEIGHT = 1;
PART_ELEMENTS_NUM = 10;
PART_NODES_NUM = PART_ELEMENTS_NUM + 1;
PART_NODES = [CYLINDER_RADIUS * ones(PART_NODES_NUM, 1), linspace(CYLINDER_HEIGHT, 0, PART_NODES_NUM)'];
PART_ELEMENTS = [1:(PART_NODES_NUM - 1); 2:PART_NODES_NUM]';

%% MATERIAL
MATERIAL_DENSITY = 10.0;
MATERIAL_DAMPING = 100.0;
MATERIAL_ROTATIONAL_INERTIA_SCALAR = 500;

MATERIAL_SHEAR_MODULUS = 500.0;
MATERIAL_BULK_MODULUS = 10000.0;
MATERIAL_SHEAR_MODULUS_TRANSVERSE = 1.0;

JOB_NAME = "Job_DOF_C4";
header = {'t', 'x', 'y'};
writecell(header, JOB_NAME + ".csv");
%% SECTION
SECTION_THICKNESS = 0.01;
T = [-1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), 0.0, 0.0];
S = [-1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3), -1.0/sqrt(3), 1.0/sqrt(3)];
W = [1.0, 1.0, 1.0, 1.0, 2.0, 2.0];

%% STEPS
STEP_DT = 5E-7;
STEP_TOTAL_TIME = 10.0;
STEP_STEPS_NUM = floor(STEP_TOTAL_TIME / STEP_DT);
STEP_PRINT = 1000;
STEP_TIME = 0;

FORCE_EXTERNAL_VALUE = 0.01;
STEP_LOAD = [0, 10.0];

%% ASSEMBLY
ASSEMBLY_ELEMENTS_NUM = PART_ELEMENTS_NUM;
ASSEMBLY_NODES_NUM = PART_NODES_NUM;
ASSEMBLY_DOF = [PART_ELEMENTS, PART_ELEMENTS + PART_NODES_NUM];
COORDINATES_INITIAL = reshape(PART_NODES, 2*PART_NODES_NUM, 1);
VELOCITY_INITIAL = zeros(2*PART_NODES_NUM, 1);
ACCELERATION_INITIAL = zeros(2*PART_NODES_NUM, 1);
DIRECTORS_INITIAL = [ones(PART_NODES_NUM, 1); zeros(PART_NODES_NUM, 1)];
ANGULAR_VELOCITY_INITIAL = zeros(PART_NODES_NUM, 1);
ANGULAR_ACCELERATION_INITIAL = zeros(PART_NODES_NUM, 1);

COORDINATES_CURRENT = reshape(PART_NODES, 2*PART_NODES_NUM, 1);
VELOCITY_CURRENT = zeros(2*PART_NODES_NUM, 1);
ACCELERATION_CURRENT = zeros(2*PART_NODES_NUM, 1);
DIRECTORS_CURRENT = [ones(PART_NODES_NUM, 1); zeros(PART_NODES_NUM, 1)];
ANGULAR_VELOCITY_CURRENT = zeros(PART_NODES_NUM, 1);
ANGULAR_ACCELERATION_CURRENT = zeros(PART_NODES_NUM, 1);

%% BC
% "TRANSLATION D" == boundary condition over translational degrees of
% freedom
% "ROTATION D" == boundary condition over rotational degrees of
% freedom

BOUNDARY_CONDITION(1).type = "TRANSLATION FIX";
BOUNDARY_CONDITION(1).index = 2 * ASSEMBLY_NODES_NUM;
BOUNDARY_CONDITION(1).value = COORDINATES_INITIAL(2 * ASSEMBLY_NODES_NUM);
%     BOUNDARY_CONDITION(2).type = "ROTATION FIX";
%     BOUNDARY_CONDITION(2).index = ASSEMBLY_NODES_NUM;
%     BOUNDARY_CONDITION(2).value = DIRECTORS_INITIAL(ASSEMBLY_NODES_NUM);
%     BOUNDARY_CONDITION(3).type = "ROTATION FIX";
%     BOUNDARY_CONDITION(3).index = 2 * ASSEMBLY_NODES_NUM;
%     BOUNDARY_CONDITION(3).value = DIRECTORS_INITIAL(2 * ASSEMBLY_NODES_NUM);
%     BOUNDARY_CONDITION(4).type = "TRANSLATION FIX";
%     BOUNDARY_CONDITION(4).index = ASSEMBLY_NODES_NUM;
%     BOUNDARY_CONDITION(4).value = COORDINATES_INITIAL(ASSEMBLY_NODES_NUM);
% *************************************************************************

%% SETUP
% *************************************************************************
%% SHAPE FUNCTIONS
%% MASS MATRIX

MASS_MATRIX_GLOBAL = zeros(2*ASSEMBLY_NODES_NUM, 1);
ROTATIONAL_INERTIA_GLOBAL = zeros(ASSEMBLY_NODES_NUM, 1);

MASS_MATRIX_LOCAL = zeros(ASSEMBLY_NODES_NUM, 4);
ROTATIONAL_INERTIA_LOCAL = zeros(ASSEMBLY_NODES_NUM, 4);
ROTATIONAL_INERTIA_LOCAL_COMPACTED = zeros(ASSEMBLY_NODES_NUM, 2);

for i = 1:ASSEMBLY_ELEMENTS_NUM

    X0_1 = COORDINATES_INITIAL(i);
    X0_2 = COORDINATES_INITIAL(i + 1);
    Y0_1 = COORDINATES_INITIAL(ASSEMBLY_NODES_NUM + i);
    Y0_2 = COORDINATES_INITIAL(ASSEMBLY_NODES_NUM + i + 1);

    PX0_1 = DIRECTORS_INITIAL(i);
    PX0_2 = DIRECTORS_INITIAL(i + 1);
    PY0_1 = DIRECTORS_INITIAL(ASSEMBLY_NODES_NUM + i);
    PY0_2 = DIRECTORS_INITIAL(ASSEMBLY_NODES_NUM + i + 1);

    for ig = 1:6

        N1 = 0.5 * (1 - T(ig));
        N2 = 0.5 * (1 + T(ig));
        N3 = 0.5 * (1 - T(ig)) * SECTION_THICKNESS * S(ig) / 2;
        N4 = 0.5 * (1 + T(ig)) * SECTION_THICKNESS * S(ig) / 2;
    
        DN1_DT = -0.5;
        DN2_DT = 0.5;
        DN3_DT = -0.5 * SECTION_THICKNESS * S(ig) / 2;
        DN4_DT = 0.5 * SECTION_THICKNESS * S(ig) / 2;
    
        DN1_DS = 0;
        DN2_DS = 0;
        DN3_DS = 0.5 * (1 - T(ig)) * SECTION_THICKNESS / 2;
        DN4_DS = 0.5 * (1 + T(ig)) * SECTION_THICKNESS / 2;

        X_IG = N1 * X0_1 + ...
               N2 * X0_2 + ...
               N3 * PX0_1 + ...
               N4 * PX0_2;
    
        DX_DT = DN1_DT * X0_1 + ...
                DN2_DT * X0_2 + ...
                DN3_DT * PX0_1 + ...
                DN4_DT * PX0_2;

        DY_DT = DN1_DT * Y0_1 + ...
                DN2_DT * Y0_2 + ...
                DN3_DT * PY0_1 + ...
                DN4_DT * PY0_2;

        DX_DS = DN1_DS * X0_1 + ...
                DN2_DS * X0_2 + ...
                DN3_DS * PX0_1 + ...
                DN4_DS * PX0_2;

        DY_DS = DN1_DS * Y0_1 + ...
                DN2_DS * Y0_2 + ...
                DN3_DS * PY0_1 + ...
                DN4_DS * PY0_2;
    
        DET_J = DX_DT .* DY_DS - DX_DS .* DY_DT;

        MASS_MATRIX_LOCAL(i, :) = MASS_MATRIX_LOCAL(i, :) + ...
            (2 * pi * MATERIAL_DENSITY * W(ig) * DET_J * X_IG) * [N1, N2, N1, N2];

        ROTATIONAL_INERTIA_LOCAL(i, :) = ROTATIONAL_INERTIA_LOCAL(i, :) + ...
            (2 * pi * MATERIAL_DENSITY * W(ig) * DET_J * X_IG) * [N3, N4, N3, N4];
    end
    ROTATIONAL_INERTIA_LOCAL_COMPACTED(i, 1) = ...
        ROTATIONAL_INERTIA_LOCAL(i, 1) * PY0_1.^2 + ...
        ROTATIONAL_INERTIA_LOCAL(i, 3) * PX0_1.^2;
    ROTATIONAL_INERTIA_LOCAL_COMPACTED(i, 2) = ...
        ROTATIONAL_INERTIA_LOCAL(i, 2) * PY0_2.^2 + ...
        ROTATIONAL_INERTIA_LOCAL(i, 4) * PX0_2.^2;
end
for i = 1:ASSEMBLY_ELEMENTS_NUM
    DOF_INDEX = [i, i + 1, ASSEMBLY_NODES_NUM + i, ASSEMBLY_NODES_NUM + i + 1];
    
    MASS_MATRIX_GLOBAL(DOF_INDEX) = ...
        MASS_MATRIX_GLOBAL(DOF_INDEX) + MASS_MATRIX_LOCAL(i, :)';
    
    ROTATIONAL_INERTIA_GLOBAL(DOF_INDEX(1:2)) = ...
        ROTATIONAL_INERTIA_GLOBAL(DOF_INDEX(1:2)) + ...
        ROTATIONAL_INERTIA_LOCAL_COMPACTED(i, :)';
end

ROTATIONAL_INERTIA_GLOBAL = ROTATIONAL_INERTIA_GLOBAL * MATERIAL_ROTATIONAL_INERTIA_SCALAR;

TEMP_VARS = zeros(100, 2);

for STEP = 1:STEP_STEPS_NUM

    TEST_VAR = 0;

    FORCE_INTERNAL_ELEMENTS = zeros(ASSEMBLY_ELEMENTS_NUM, 4);
    ROTATIONAL_FORCE_INTERNAL_ELEMENTS = zeros(ASSEMBLY_ELEMENTS_NUM, 4);
    TORQUE_INTERNAL_ELEMENTS = zeros(ASSEMBLY_ELEMENTS_NUM, 2);

    FORCE_INTERNAL_GLOBAL = zeros(2 * ASSEMBLY_NODES_NUM, 1);
    FORCE_EXTERNAL_GLOBAL = zeros(2 * ASSEMBLY_NODES_NUM, 1);
    TORQUE_GLOBAL = zeros(ASSEMBLY_NODES_NUM, 1);


    for i = 1:ASSEMBLY_ELEMENTS_NUM

        X_1 = COORDINATES_CURRENT(i);
        X_2 = COORDINATES_CURRENT(i + 1);
        Y_1 = COORDINATES_CURRENT(ASSEMBLY_NODES_NUM + i);
        Y_2 = COORDINATES_CURRENT(ASSEMBLY_NODES_NUM + i + 1);

        PX_1 = DIRECTORS_CURRENT(i);
        PX_2 = DIRECTORS_CURRENT(i + 1);
        PY_1 = DIRECTORS_CURRENT(ASSEMBLY_NODES_NUM + i);
        PY_2 = DIRECTORS_CURRENT(ASSEMBLY_NODES_NUM + i + 1);

        X0_1 = COORDINATES_INITIAL(i);
        X0_2 = COORDINATES_INITIAL(i + 1);
        Y0_1 = COORDINATES_INITIAL(ASSEMBLY_NODES_NUM + i);
        Y0_2 = COORDINATES_INITIAL(ASSEMBLY_NODES_NUM + i + 1);

        PX0_1 = DIRECTORS_INITIAL(i);
        PX0_2 = DIRECTORS_INITIAL(i + 1);
        PY0_1 = DIRECTORS_INITIAL(ASSEMBLY_NODES_NUM + i);
        PY0_2 = DIRECTORS_INITIAL(ASSEMBLY_NODES_NUM + i + 1);

        for ig = 1:4

            N1 = 0.5 * (1 - T(ig));
            N2 = 0.5 * (1 + T(ig));
            N3 = 0.5 * (1 - T(ig)) * SECTION_THICKNESS * S(ig) / 2;
            N4 = 0.5 * (1 + T(ig)) * SECTION_THICKNESS * S(ig) / 2;
        
            DN1_DT = -0.5;
            DN2_DT = 0.5;
            DN3_DT = -0.5 * SECTION_THICKNESS * S(ig) / 2;
            DN4_DT = 0.5 * SECTION_THICKNESS * S(ig) / 2;
        
            DN1_DS = 0;
            DN2_DS = 0;
            DN3_DS = 0.5 * (1 - T(ig)) * SECTION_THICKNESS / 2;
            DN4_DS = 0.5 * (1 + T(ig)) * SECTION_THICKNESS / 2;
    
            X_IG = N1 * X_1 + ...
                   N2 * X_2 + ...
                   N3 * PX_1 + ...
                   N4 * PX_2;
        
            DX_DT = DN1_DT * X_1 + ...
                    DN2_DT * X_2 + ...
                    DN3_DT * PX_1 + ...
                    DN4_DT * PX_2;
    
            DY_DT = DN1_DT * Y_1 + ...
                    DN2_DT * Y_2 + ...
                    DN3_DT * PY_1 + ...
                    DN4_DT * PY_2;

            Q11 = DX_DT / sqrt(DX_DT^2 + DY_DT^2);
            Q12 = DY_DT / sqrt(DX_DT^2 + DY_DT^2);
            Q21 = -Q12;
            Q22 = Q11;

            X_1_ROT = Q11 * X_1 + Q12 * Y_1;
            Y_1_ROT = Q21 * X_1 + Q22 * Y_1;

            X_2_ROT = Q11 * X_2 + Q12 * Y_2;
            Y_2_ROT = Q21 * X_2 + Q22 * Y_2;

            X0_1_ROT = Q11 * X0_1 + Q12 * Y0_1;
            Y0_1_ROT = Q21 * X0_1 + Q22 * Y0_1;

            X0_2_ROT = Q11 * X0_2 + Q12 * Y0_2;
            Y0_2_ROT = Q21 * X0_2 + Q22 * Y0_2;

            PX_1_ROT = Q11 * PX_1 + Q12 * PY_1;
            PY_1_ROT = Q21 * PX_1 + Q22 * PY_1;

            PX_2_ROT = Q11 * PX_2 + Q12 * PY_2;
            PY_2_ROT = Q21 * PX_2 + Q22 * PY_2;

            PX0_1_ROT = Q11 * PX0_1 + Q12 * PY0_1;
            PY0_1_ROT = Q21 * PX0_1 + Q22 * PY0_1;

            PX0_2_ROT = Q11 * PX0_2 + Q12 * PY0_2;
            PY0_2_ROT = Q21 * PX0_2 + Q22 * PY0_2;

            DX_DT = DN1_DT * X_1_ROT + ...
                    DN2_DT * X_2_ROT + ...
                    DN3_DT * PX_1_ROT + ...
                    DN4_DT * PX_2_ROT;
    
            DY_DT = DN1_DT * Y_1_ROT + ...
                    DN2_DT * Y_2_ROT + ...
                    DN3_DT * PY_1_ROT + ...
                    DN4_DT * PY_2_ROT;
    
            DX_DS = DN1_DS * X_1_ROT + ...
                    DN2_DS * X_2_ROT + ...
                    DN3_DS * PX_1_ROT + ...
                    DN4_DS * PX_2_ROT;
    
            DY_DS = DN1_DS * Y_1_ROT + ...
                    DN2_DS * Y_2_ROT + ...
                    DN3_DS * PY_1_ROT + ...
                    DN4_DS * PY_2_ROT;
            
            
            DET_J = DX_DT * DY_DS - DX_DS * DY_DT;

            DT_DX = DY_DS / DET_J;
            DT_DY = -DX_DS / DET_J;
            DS_DX = -DY_DT / DET_J;
            DS_DY = DX_DT / DET_J;
            
            DN1_DX = DT_DX * DN1_DT + DS_DX * DN1_DS;
            DN2_DX = DT_DX * DN2_DT + DS_DX * DN2_DS;
            DN3_DX = DT_DX * DN3_DT + DS_DX * DN3_DS;
            DN4_DX = DT_DX * DN4_DT + DS_DX * DN4_DS;

            DN1_DY = DT_DY * DN1_DT + DS_DY * DN1_DS;
            DN2_DY = DT_DY * DN2_DT + DS_DY * DN2_DS;
            DN3_DY = DT_DY * DN3_DT + DS_DY * DN3_DS;
            DN4_DY = DT_DY * DN4_DT + DS_DY * DN4_DS;
            
            B51 = Q11 * N1 / X_IG;
            B52 = Q11 * N2 / X_IG;
            B53 = Q21 * N1 / X_IG;
            B54 = Q21 * N2 / X_IG;
            B55 = Q11 * N3 / X_IG;
            B56 = Q11 * N4 / X_IG;
            B57 = Q21 * N3 / X_IG;
            B58 = Q21 * N4 / X_IG;

            B = [DN1_DX DN2_DX 0      0      DN3_DX DN4_DX 0      0;
                 DN1_DY DN2_DY 0      0      DN3_DY DN4_DY 0      0;
                 0      0      DN1_DX DN2_DX 0      0      DN3_DX DN4_DX;
                 0      0      DN1_DY DN2_DY 0      0      DN3_DY DN4_DY;
                 B51    B52    B53    B54    B55    B56    B57    B58];
            
            DOF_0 = [X0_1_ROT;
                     X0_2_ROT;
                     Y0_1_ROT;
                     Y0_2_ROT;
                     PX0_1_ROT;
                     PX0_2_ROT;
                     PY0_1_ROT;
                     PY0_2_ROT];

            F_INV_11 = B(1, :) * DOF_0;
            F_INV_12 = B(2, :) * DOF_0;
            F_INV_21 = B(3, :) * DOF_0;
            F_INV_22 = B(4, :) * DOF_0;
            F_INV_33 = B(5, :) * DOF_0;

            DET_F_INV = F_INV_11 * F_INV_22 * F_INV_33 - ...
                        F_INV_12 * F_INV_21 * F_INV_33;

            F_GRAD_11 = F_INV_22 * F_INV_33 / DET_F_INV;
            F_GRAD_12 = -F_INV_12 * F_INV_33 / DET_F_INV;
            F_GRAD_21 = -F_INV_21 * F_INV_33 / DET_F_INV;
            F_GRAD_22 = F_INV_11 * F_INV_33 / DET_F_INV;
            F_GRAD_33 = (-F_INV_12 * F_INV_21 + F_INV_11 * F_INV_22) / DET_F_INV;

%             B_STRAIN_11 = F_GRAD_11^2 + F_GRAD_12^2;
%             B_STRAIN_12 = F_GRAD_11 * F_GRAD_21 + F_GRAD_12 * F_GRAD_22;
%             B_STRAIN_21 = F_GRAD_11 * F_GRAD_21 + F_GRAD_12 * F_GRAD_22;
%             B_STRAIN_22 = F_GRAD_21^2 + F_GRAD_22^2;
%             B_STRAIN_33 = F_GRAD_33^2;
% 
%             VOL_J = sqrt(B_STRAIN_11 * (B_STRAIN_22 * F_INV_33) - B_STRAIN_12 * (B_STRAIN_21 * B_STRAIN_33));
% 
%             B_STAR_11 = VOL_J^(-2.0/3.0)*B_STRAIN_11;
%             B_STAR_12 = VOL_J^(-2.0/3.0)*B_STRAIN_12;
%             B_STAR_21 = VOL_J^(-2.0/3.0)*B_STRAIN_21;
%             B_STAR_22 = VOL_J^(-2.0/3.0)*B_STRAIN_22;
%             B_STAR_33 = VOL_J^(-2.0/3.0)*B_STRAIN_33;
%                 
%                 T1 = (B_STAR_11 + B_STAR_22 + B_STAR_33)/3;
%                 T2 = MATERIAL_SHEAR_MODULUS/VOL_J;
%                 T3 = MATERIAL_BULK_MODULUS*(VOL_J - 1.0);
% 
%                 if ig < 5
%                     SIGMA_11 = T2*(B_STAR_11 - T1) + T3;
%                     SIGMA_12 = 0;
%                     SIGMA_21 = 0;
%                     SIGMA_22 = T2*(B_STAR_22 - T1) + T3;
%                     SIGMA_33 = T2*(B_STAR_33 - T1) + T3;
%                 else
%                     SIGMA_11 = 0;
%                     SIGMA_12 = MATERIAL_SHEAR_MODULUS_TRANSVERSE*B_STAR_12;
%                     SIGMA_21 = MATERIAL_SHEAR_MODULUS_TRANSVERSE*B_STAR_21;
%                     SIGMA_22 = 0;
%                     SIGMA_33 = 0;
%                 end

            EPSILON_11 = F_GRAD_11 - 1.0;  % Axial strain
            EPSILON_33 = F_GRAD_33 - 1.0;  % Hoop strain

            POISSON_RATIO = (3 * MATERIAL_BULK_MODULUS - 2 * MATERIAL_SHEAR_MODULUS) / ...
                (6 * MATERIAL_BULK_MODULUS + 2 * MATERIAL_SHEAR_MODULUS);
            YOUNGS_MODULUS = 2 * MATERIAL_SHEAR_MODULUS*(1 + POISSON_RATIO);
            FACTOR = YOUNGS_MODULUS / (1 - POISSON_RATIO^2);

            if ig < 5
                SIGMA_11 = FACTOR * (EPSILON_11 + POISSON_RATIO * EPSILON_33);  % Axial stress
                SIGMA_12 = 0;
                SIGMA_21 = 0;
                SIGMA_22 = 0;
                SIGMA_33 = FACTOR * (EPSILON_33 + POISSON_RATIO * EPSILON_11);  % Hoop stress
            end


            FORCE_LOCAL_ROT = B' * [SIGMA_11; SIGMA_12; SIGMA_21; SIGMA_22; SIGMA_33] * DET_J * X_IG * 2 * pi * W(ig);

            FORCE_LOCAL = zeros(8, 1);

            FORCE_LOCAL(1) = Q11 * FORCE_LOCAL_ROT(1) + ...
                             Q21 * FORCE_LOCAL_ROT(3);
            FORCE_LOCAL(3) = Q12 * FORCE_LOCAL_ROT(1) + ...
                             Q22 * FORCE_LOCAL_ROT(3);
            FORCE_LOCAL(2) = Q11 * FORCE_LOCAL_ROT(2) + ...
                             Q21 * FORCE_LOCAL_ROT(4);
            FORCE_LOCAL(4) = Q12 * FORCE_LOCAL_ROT(2) + ...
                             Q22 * FORCE_LOCAL_ROT(4);
            FORCE_LOCAL(5) = Q11 * FORCE_LOCAL_ROT(5) + ...
                             Q21 * FORCE_LOCAL_ROT(7);
            FORCE_LOCAL(7) = Q12 * FORCE_LOCAL_ROT(5) + ...
                             Q22 * FORCE_LOCAL_ROT(7);
            FORCE_LOCAL(6) = Q11 * FORCE_LOCAL_ROT(6) + ...
                             Q21 * FORCE_LOCAL_ROT(8);
            FORCE_LOCAL(8) = Q12 * FORCE_LOCAL_ROT(6) + ...
                             Q22 * FORCE_LOCAL_ROT(8);
            
            FORCE_IG = FORCE_LOCAL(1:4);
            TORQUE_IG = FORCE_LOCAL(5:8);

            FORCE_INTERNAL_ELEMENTS(i, :) = ...
                    FORCE_INTERNAL_ELEMENTS(i, :) + FORCE_IG';
            ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, :) = ...
                    ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, :) + TORQUE_IG';
            
            if (mod(STEP, STEP_PRINT) == 0 || STEP == 1) && i == 1 && ig == 1
                fprintf('Step %d: DET_J = %.2e\n', STEP, DET_J);
                fprintf('B-matrix condition number: %.2e\n', cond(B));
                fprintf('Q rotation matrix: [%.3f %.3f; %.3f %.3f]\n', Q11, Q12, Q21, Q22);
                fprintf('F_INV components: [%.3f %.3f %.3f %.3f %.3f]\n', ...
                    F_INV_11, F_INV_12, F_INV_21, F_INV_22, F_INV_33);
            end

        end
        TORQUE_INTERNAL_ELEMENTS(i, 1) = -ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, 1) * PY_1 + ...
                                           ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, 3) * PX_1;
        TORQUE_INTERNAL_ELEMENTS(i, 2) = -ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, 2) * PY_2 + ...
                                           ROTATIONAL_FORCE_INTERNAL_ELEMENTS(i, 4) * PX_2;
    end
    for i = 1:ASSEMBLY_ELEMENTS_NUM
        DOF_INDEX = ASSEMBLY_DOF(i, :);
        FORCE_INTERNAL_GLOBAL(DOF_INDEX) = FORCE_INTERNAL_GLOBAL(DOF_INDEX) + FORCE_INTERNAL_ELEMENTS(i, :)';
        TORQUE_GLOBAL(DOF_INDEX(1:2)) = TORQUE_GLOBAL(DOF_INDEX(1:2)) + TORQUE_INTERNAL_ELEMENTS(i, :)';
    end

    if any(~isreal(FORCE_INTERNAL_GLOBAL))
        fprintf("COMPLEX appeared: t = %g\n", STEP_TIME);
        break
    end

    
    if STEP_TIME > STEP_LOAD(1) && STEP_TIME < STEP_LOAD(2)
        FORCE_EXTERNAL_GLOBAL(ASSEMBLY_NODES_NUM + 1) = FORCE_EXTERNAL_VALUE;
    else
        FORCE_EXTERNAL_GLOBAL(ASSEMBLY_NODES_NUM + 1) = 0;
    end

    for i = 1:(2 * ASSEMBLY_NODES_NUM)
%             FORCE_INTERNAL_GLOBAL(i) = FORCE_INTERNAL_GLOBAL(i) - ...
%                MATERIAL_DAMPING * MASS_MATRIX_GLOBAL(i) * VELOCITY_CURRENT(i);
        ACCELERATION_PREVIOUS = ACCELERATION_CURRENT(i);
        ACCELERATION_CURRENT(i) = (FORCE_EXTERNAL_GLOBAL(i) - FORCE_INTERNAL_GLOBAL(i)) / MASS_MATRIX_GLOBAL(i) - ...
            MATERIAL_DAMPING * VELOCITY_CURRENT(i);
        VELOCITY_CURRENT(i) = VELOCITY_CURRENT(i) + ...
            0.5 * STEP_DT * (ACCELERATION_CURRENT(i) + ACCELERATION_PREVIOUS);
        
        COORDINATES_CURRENT(i) = COORDINATES_CURRENT(i) + ...
            VELOCITY_CURRENT(i) * STEP_DT + ...
            0.5 * ACCELERATION_PREVIOUS * STEP_DT^2;
    end

    for i = 1:ASSEMBLY_NODES_NUM
%             TORQUE_GLOBAL(i) = TORQUE_GLOBAL(i) - ...
%                 MATERIAL_DAMPING * ROTATIONAL_INERTIA_GLOBAL(i) * ANGULAR_VELOCITY_CURRENT(i);
        ANGULAR_ACCELERATION_PREVIOUS = ANGULAR_ACCELERATION_CURRENT(i);
        ANGULAR_ACCELERATION_CURRENT(i) = TORQUE_GLOBAL(i) / ROTATIONAL_INERTIA_GLOBAL(i) - ...
            MATERIAL_DAMPING * ANGULAR_VELOCITY_CURRENT(i);
        ANGULAR_VELOCITY_CURRENT(i) = ANGULAR_VELOCITY_CURRENT(i) + ...
            0.5 * STEP_DT * (ANGULAR_ACCELERATION_CURRENT(i) + ANGULAR_ACCELERATION_PREVIOUS);
        ROTATION_ANGLE = ANGULAR_VELOCITY_CURRENT(i) * STEP_DT + ...
            0.5 * ANGULAR_ACCELERATION_PREVIOUS * STEP_DT^2;
        
        PX = DIRECTORS_CURRENT(i);
        PY = DIRECTORS_CURRENT(i + ASSEMBLY_NODES_NUM);
        DIRECTORS_CURRENT(i) = ...
            cos(ROTATION_ANGLE) * PX - ...
            sin(ROTATION_ANGLE) * PY;
        DIRECTORS_CURRENT(i + ASSEMBLY_NODES_NUM) = ...
            sin(ROTATION_ANGLE) * PX + ...
            cos(ROTATION_ANGLE) * PY;
        DIRECTOR_NORMALIZATION_FACTOR = ...
            sqrt(DIRECTORS_CURRENT(i)^2 + ...
            DIRECTORS_CURRENT(i + ASSEMBLY_NODES_NUM)^2);
        DIRECTORS_CURRENT(i) = ...
            DIRECTORS_CURRENT(i) / DIRECTOR_NORMALIZATION_FACTOR;
        DIRECTORS_CURRENT(i + ASSEMBLY_NODES_NUM) = ...
            DIRECTORS_CURRENT(i + ASSEMBLY_NODES_NUM) / DIRECTOR_NORMALIZATION_FACTOR;
    end

    for i = 1:length(BOUNDARY_CONDITION)
      if BOUNDARY_CONDITION(i).type == "TRANSLATION FIX"
          COORDINATES_CURRENT(BOUNDARY_CONDITION(i).index) = BOUNDARY_CONDITION(i).value;
      elseif BOUNDARY_CONDITION(i).type == "ROTATION FIX"
          DIRECTORS_CURRENT(BOUNDARY_CONDITION(i).index) = BOUNDARY_CONDITION(i).value;
      end
    end
    if mod(STEP, STEP_PRINT) == 0
        fprintf("t = %g,\tx = %g,\ty = %g\tfx = %g\tfy = %g\n", ...
            STEP_TIME, ...
            COORDINATES_CURRENT(1), ...
            COORDINATES_CURRENT(ASSEMBLY_NODES_NUM + 1), ...
            sum(FORCE_INTERNAL_GLOBAL(1:ASSEMBLY_NODES_NUM)), ...
            sum(FORCE_INTERNAL_GLOBAL((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM)));
        fprintf("--------------------------------------------------------------\n")
        TEMP_VARS(floor(STEP / STEP_PRINT), 1) = COORDINATES_CURRENT(1);
        TEMP_VARS(floor(STEP / STEP_PRINT), 2) = COORDINATES_CURRENT(ASSEMBLY_NODES_NUM + 1);
        plot(COORDINATES_CURRENT(1:ASSEMBLY_NODES_NUM), COORDINATES_CURRENT((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM), 'r')
        hold on
        quiver(COORDINATES_CURRENT(1:ASSEMBLY_NODES_NUM), COORDINATES_CURRENT((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM), DIRECTORS_CURRENT(1:ASSEMBLY_NODES_NUM), DIRECTORS_CURRENT((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM), 'Color','b')
        quiver(COORDINATES_CURRENT(1:ASSEMBLY_NODES_NUM), ...
               COORDINATES_CURRENT((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM), ...
               FORCE_INTERNAL_GLOBAL(1:ASSEMBLY_NODES_NUM), ...
               FORCE_INTERNAL_GLOBAL((ASSEMBLY_NODES_NUM + 1):2*ASSEMBLY_NODES_NUM), 'Color','g')
        ylim([-1, 2]);
        xlim([0 2])
        daspect([1 1 1])
        grid on
        drawnow
        hold off
        writematrix([STEP_TIME, ...
                     COORDINATES_CURRENT(1), ...
                     COORDINATES_CURRENT(ASSEMBLY_NODES_NUM + 1)], ...
                     JOB_NAME + ".csv" ...
                     , 'WriteMode', 'append');
    end
    STEP_TIME = STEP_TIME + STEP_DT;
end