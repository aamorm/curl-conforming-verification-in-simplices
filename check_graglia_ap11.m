% See notes in 18.03.2019
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
% Here, a = 1; b = 2; c = 3; d = 4. I don't know which edge/face is
% but I know that they need to be mutually exclusive
xi12 = L1-L2;
xi34 = L3-L4;
chi12 = L1+L2;
chi34 = L3+L4;
chi2 = L1^2-4*L1*L2+L2^2;
chi3 = xi12*(L1^2-8*L1*L2+L2^2);
chi4 = L1^4-16*L1^3*L2+36*L1^2*L2^2-16*L1*L2^3+L2^4;
chi5 = xi12*(L1^4-24*L1^3*L2+76*L1^2*L2^2-24*L1*L2^3+L2^4);
% Legendre polynomials
legendre12_o2 = 1/2*(3*xi12^2-1);
legendre12_o3 = 1/2*(5*xi12^3-3*xi12);
legendre12_o4 = 1/8*(35*xi12^4-30*xi12^2+3);
legendre12_o5 = 1/8*(63*xi12^5-70*xi12^3+15*xi12);
legendre12_o6 = 1/16*(231*xi12^6-315*xi12^4+105*xi12^2-5);
% This belongs to Z axis (edge 12)
whitney12 = L1*gL2-L2*gL1;
% This belongs to Y axis (edge 13)
whitney13 = L1*gL3-L3*gL1;
% This belongs to X axis (edge 13)
whitney23 = L2*gL3-L3*gL2;
% This belongs to X+Y = 1 axis (edge 34)
whitney34 = L3*gL4-L4*gL3;
tolerance = 1e-10;
detailed_output = false;

disp('Edge functions')
homogeneous_order = 1;
e0_ss = whitney12;
good_function = check_space(e0_ss, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E0 function is in Nedelec space")
else
    disp("E0 function is not in Nedelec space")
end
homogeneous_order = 2;
e1_as = sqrt(3)*xi12*whitney12;
good_function = check_space(e1_as, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E1 function is in Nedelec space")
else
    disp("E1 function is not in Nedelec space")
end
homogeneous_order = 3;
e2_ss = sqrt(5)*(legendre12_o2-chi34*(chi34-2)/2)*whitney12;
good_function = check_space(e2_ss, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E2 function is in Nedelec space")
else
    disp("E2 function is not in Nedelec space")
end
homogeneous_order = 4;
e3_as = sqrt(7)*(legendre12_o3-3*chi34*(chi34-2)*xi12/2)*whitney12;
good_function = check_space(e3_as, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E3 function is in Nedelec space")
else
    disp("E3 function is not in Nedelec space")
end
homogeneous_order = 5;
e4_ss = sqrt(9)*(legendre12_o4+3*chi34*(chi34-2)*(1+40*L1*L2-9*chi12^2)/8)*whitney12;
good_function = check_space(e4_ss, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E4 function is in Nedelec space")
else
    disp("E4 function is not in Nedelec space")
end
homogeneous_order = 6;
e5_as = sqrt(11)*(legendre12_o5+5*chi34*(chi34-2)*xi12*(3+56*L1*L2-11*chi12^2)/8)*whitney12;
good_function = check_space(e5_as, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E5 function is in Nedelec space")
else
    disp("E5 function is not in Nedelec space")
end
homogeneous_order = 7;
e6_as = sqrt(13)*(legendre12_o6-5*chi34*(chi34-2)*(1+43*chi12^4+84*L1*L2*(1+12*L1*L2)-20*chi12^2*(1+21*L1*L2))/16)*whitney12;
good_function = check_space(e6_as, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("E6 function is in Nedelec space")
else
    disp("E6 function is not in Nedelec space")
end

disp('Face functions')
% The polynomial needs to be multiplied by two whitney functions in the face.
% The face chosen is the XY face, which means that L4 = 0.
% Only one edge is taken, the other is just a permutation.
% For face L4 = 0, the possible indices (edges) for
% whitney functions are 23, 34, 13.
% Here, whitney23 is taken, so c=2; d=3; a=1; b=4.
homogeneous_order = 2;
f01_s = 2*sqrt(3)*L3*whitney34;
good_function = check_space(f01_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F01 function is in Nedelec space")
else
    disp("One F01 function is not in Nedelec space")
end

homogeneous_order = 3;
f02_s = 2*sqrt(3)*L3*(5*L3-3+3*L4)*whitney34;
good_function = check_space(f02_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F02 function is in Nedelec space")
else
    disp("One F02 function is not in Nedelec space")
end
f11_a = 6*sqrt(5)*L3*xi12*whitney34;
good_function = check_space(f11_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F11 function is in Nedelec space")
else
    disp("One F11 function is not in Nedelec space")
end

homogeneous_order = 4;
f03_s = 2*sqrt(30)*L3*(7*L3^2-8*L3+2+2*L4*(4*L3-2+L4))*whitney34;
good_function = check_space(f03_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F03 function is in Nedelec space")
else
    disp("One F03 function is not in Nedelec space")
end
f12_a = 2*sqrt(30)*L3*xi12*(7*L3-3+3*L4)*whitney34;
good_function = check_space(f12_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F12 function is in Nedelec space")
else
    disp("One F12 function is not in Nedelec space")
end
f21_s = 2*sqrt(210)*L3*chi2*whitney34;
good_function = check_space(f21_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F21 function is in Nedelec space")
else
    disp("One F21 function is not in Nedelec space")
end

homogeneous_order = 5;
f04_s = 2*sqrt(15)*L3*(42*L3^3-70*L3^2+35*L3-5+5*L4*(1+chi12+chi12^2-L3*(6+5*chi12-8*L3)))*whitney34;
good_function = check_space(f04_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F04 function is in Nedelec space")
else
    disp("One F04 function is not in Nedelec space")
end
f13_a = 2*sqrt(105)*L3*xi12*(18*L3^2-16*L3+3-L4*(3+3*chi12-13*L3))*whitney34;
good_function = check_space(f13_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F13 function is in Nedelec space")
else
    disp("One F13 function is not in Nedelec space")
end
f22_s = 10*sqrt(42)*L3*chi2*(3*L3-1+L4)*whitney34;
good_function = check_space(f22_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F22 function is in Nedelec space")
else
    disp("One F22 function is not in Nedelec space")
end
f31_a = 6*sqrt(70)*L3*chi3*whitney34;
good_function = check_space(f31_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F31 function is in Nedelec space")
else
    disp("One F31 function is not in Nedelec space")
end

homogeneous_order = 6;
f05_s = 2*sqrt(105)*L3*(66*L3^4-144*L3^3+108*L3^2-32*L3+3-L4*(3*(1+chi12)*(1+chi12^2)-L3*(29+chi12*(26+23*chi12))+L3^2*(79+53*chi12)+65*L3^3))*whitney34;
good_function = check_space(f05_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F05 function is in Nedelec space")
else
    disp("One F05 function is not in Nedelec space")
end
f14_a = 6*sqrt(70)*L3*xi12*(33*L3^3-45*L3^2+18*L3-2+L4*(2*(1+chi12^2+chi12*(1-7*L3))-16*L3+29*L3^2))*whitney34;
good_function = check_space(f14_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F14 function is in Nedelec space")
else
    disp("One F14 function is not in Nedelec space")
end
f23_s = 6*sqrt(10)*L3*chi2*(55*L3^2-40*L3+6+2*L4*(3*L4+20*L3-6))*whitney34;
good_function = check_space(f23_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F23 function is in Nedelec space")
else
    disp("One F23 function is not in Nedelec space")
end
f32_a = 6*sqrt(35)*L3*chi3*(11*L3-3+3*L4)*whitney34;
good_function = check_space(f32_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F32 function is in Nedelec space")
else
    disp("One F32 function is not in Nedelec space")
end
f41_s = 6*sqrt(165)*L3*chi4*whitney34;
good_function = check_space(f41_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F41 function is in Nedelec space")
else
    disp("One F41 function is not in Nedelec space")
end

homogeneous_order = 7;
f06_s = 2*sqrt(42)*L3*(429*L3^5-1155*L3^4+1155*L3^3-525*L3^2+105*L3-7+7*L4*(1+chi12+chi12^2+chi12^3+chi12^4-(14+13*chi12+12*chi12^2+11*chi12^3)*L3+(61+48*chi12+36*chi12^2)*L3^2-8*(13+7*chi12)*L3^3+61*L3^4))*whitney34;
good_function = check_space(f06_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F06 function is in Nedelec space")
else
    disp("One F06 function is not in Nedelec space")
end
f15_a = 6*sqrt(70)*L3*xi12*(143*L3^4-264*L3^3+165*L3^2-40*L3+3+L4*(136*L3^3-3*(1+chi12+chi12^2+chi12^3)+(37+34*chi12+31*chi12^2)*L3-2*(64+47*chi12)*L3^2))*whitney34;
good_function = check_space(f15_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F15 function is in Nedelec space")
else
    disp("One F15 function is not in Nedelec space")
end
f24_s = 6*sqrt(35)*L3*chi2*(143*L3^3-165*L3^2+55*L3-5+5*L4*(1+chi12+chi12^2-(10+9*chi12)*L3+23*L3^2))*whitney34;
good_function = check_space(f24_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F24 function is in Nedelec space")
else
    disp("One F24 function is not in Nedelec space")
end
f33_a = 14*sqrt(165)*L3*chi3*(13*L3^2-8*L3+1-L4*(2-8*L3-L4))*whitney34;
good_function = check_space(f33_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F33 function is in Nedelec space")
else
    disp("One F33 function is not in Nedelec space")
end
f42_s = 6*sqrt(77)*L3*chi4*(13*L3-3+3*L4)*whitney34;
good_function = check_space(f42_s, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F42 function is in Nedelec space")
else
    disp("One F42 function is not in Nedelec space")
end
f51_a = 2*sqrt(3003)*L3*chi5*whitney34;
good_function = check_space(f51_a, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One F51 function is in Nedelec space")
else
    disp("One F51 function is not in Nedelec space")
end

disp('Volume functions')
% Same as the face. With whitney34, we have a=1, b=2, c=3, d=4.
u000_ss = 6*sqrt(35)*whitney34;

u200_ss = 12*sqrt(55)*(1+3*chi12*(5*chi12-3))*whitney34;
u110_as = 6*sqrt(2310)*xi12*(5*chi12-2)*whitney34;
u020_ss = 15*sqrt(462)*(3*xi12^2-chi12^2)*whitney34;
u101_sa = 6*sqrt(1155)*xi34*(5*chi12-1)*whitney34;
u011_aa = 60*sqrt(231)*xi12*xi34*whitney34;
u002_ss = 15*sqrt(11)*(7*xi34^2-chi34^2)*whitney34;

u300_ss = 6*sqrt(390)*(5*chi12*(3-11*chi34*chi12)-1)*whitney34;
u210_as = 6*sqrt(390)*xi12*(10+11*chi12*(6*chi12-5))*whitney34;
u120_ss = 15*sqrt(6006)*(2*chi12-1)*(3*xi12^2-chi12^2)*whitney34;
u030_as = 6*sqrt(15015)*xi12*(5*xi12^2-3*chi12^2)*whitney34;
u201_sa = 30*sqrt(91)*xi34*(1+11*chi12*(2*chi12-1))*whitney34;
u111_aa = 30*sqrt(6006)*xi12*xi34*(3*chi12-1)*whitney34;
u021_sa = 30*sqrt(3003)*xi34*(3*xi12^2-chi12)*whitney34;
u102_ss = 3*sqrt(715)*(6*chi12-1)*(7*xi34^2-chi34^2)*whitney34;
u012_as = 15*sqrt(858)*xi12*(7*xi34^2-chi34^2)*whitney34;
u003_sa = 3*sqrt(10010)*xi34*(3*xi34^2-chi34^2)*whitney34;
u100_ss = 6*sqrt(105)*(4*chi12-1)*whitney34;
u010_as = 36*sqrt(35)*xi12*whitney34;
u001_sa = 12*sqrt(105)*xi34*whitney34;

homogeneous_order = 3;
v200 = L3*L4*u000_ss;
good_function = check_space(v200, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V200 function is in Nedelec space")
else
    disp("One V200 function is not in Nedelec space")
end

homogeneous_order = 4;
v300 = L3*L4*u100_ss;
good_function = check_space(v300, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V300 function is in Nedelec space")
else
    disp("One V300 function is not in Nedelec space")
end
v210 = L3*L4*u010_as;
good_function = check_space(v210, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V210 function is in Nedelec space")
else
    disp("One V210 function is not in Nedelec space")
end
v201 = L3*L4*u001_sa;
good_function = check_space(v201, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V201 function is in Nedelec space")
else
    disp("One V201 function is not in Nedelec space")
end

homogeneous_order = 5;
v400 = L3*L4*u200_ss;
good_function = check_space(v400, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V400 function is in Nedelec space")
else
    disp("One V400 function is not in Nedelec space")
end
v310 = L3*L4*u110_as;
good_function = check_space(v310, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V310 function is in Nedelec space")
else
    disp("One V310 function is not in Nedelec space")
end
v220 = L3*L4*u020_ss;
good_function = check_space(v220, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V220 function is in Nedelec space")
else
    disp("One V220 function is not in Nedelec space")
end
v202 = L3*L4*u002_ss;
good_function = check_space(v202, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V202 function is in Nedelec space")
else
    disp("One V202 function is not in Nedelec space")
end
v301 = L3*L4*u101_sa;
good_function = check_space(v301, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V301 function is in Nedelec space")
else
    disp("One V301 function is not in Nedelec space")
end

homogeneous_order = 6;
v500 = L3*L4*u300_ss;
good_function = check_space(v500, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V500 function is in Nedelec space")
else
    disp("One V500 function is not in Nedelec space")
end
v410 = L3*L4*u210_as;
good_function = check_space(v410, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V410 function is in Nedelec space")
else
    disp("One V410 function is not in Nedelec space")
end
v401 = L3*L4*u201_sa;
good_function = check_space(v401, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V401 function is in Nedelec space")
else
    disp("One V401 function is not in Nedelec space")
end
v320 = L3*L4*u120_ss;
good_function = check_space(v320, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V320 function is in Nedelec space")
else
    disp("One V320 function is not in Nedelec space")
end
v302 = L3*L4*u102_ss;
good_function = check_space(v302, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V302 function is in Nedelec space")
else
    disp("One V302 function is not in Nedelec space")
end
v203 = L3*L4*u003_sa;
good_function = check_space(v203, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V203 function is in Nedelec space")
else
    disp("One V203 function is not in Nedelec space")
end
v230 = L3*L4*u030_as;
good_function = check_space(v230, homogeneous_order, affine_coordinates, tolerance, detailed_output);
if (good_function)
    disp("One V230 function is in Nedelec space")
else
    disp("One V230 function is not in Nedelec space")
end