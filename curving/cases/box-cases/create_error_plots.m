# this file is not a function
THIS_IS_NOT_A_FUNCTION=1;
clear THIS_IS_NOT_A_FUNCTION;

# function to append norms from a file
function [l2,h1] = append_file(l2, h1, file_name)

norms = load(file_name);
norms = norms';

l2 = [l2 norms(:,1)];
h1 = [h1 norms(:,2)];

endfunction

# define the norms
l2_quad = zeros(3,0);
h1_quad = l2_quad;
l2_cube = l2_quad;
h1_cube = l2_quad;
N = [5 10 20];

# read the norms
[l2_quad, h1_quad] = append_file(l2_quad, h1_quad, "box_n5_o2.pp");
[l2_quad, h1_quad] = append_file(l2_quad, h1_quad, "box_n10_o2.pp");
[l2_quad, h1_quad] = append_file(l2_quad, h1_quad, "box_n20_o2.pp");

[l2_cube, h1_cube] = append_file(l2_cube, h1_cube, "box_n5_o3.pp");
[l2_cube, h1_cube] = append_file(l2_cube, h1_cube, "box_n10_o3.pp");
[l2_cube, h1_cube] = append_file(l2_cube, h1_cube, "box_n20_o3.pp");


# ---------------------------------------------------------------------------
#                                      L2
# ---------------------------------------------------------------------------
figure(1), clf, hold on;

# quadratic
loglog(N, l2_quad(1,:), 'b-o',"linewidth", 2, 'markersize', 10);
loglog(N, l2_quad(2,:), 'b-s',"linewidth", 2, 'markersize', 10);
loglog(N, l2_quad(3,:), 'b-^',"linewidth", 2,  'markersize', 10);
loglog(N, N.^(-3)*100, "r--", 'linewidth', 2);

# cubic
loglog(N, l2_cube(1,:), 'k-o',"linewidth", 2,  'markersize', 10);
loglog(N, l2_cube(2,:), 'k-s',"linewidth", 2,  'markersize', 10);
loglog(N, l2_cube(3,:), 'k-^',"linewidth", 2,  'markersize', 10);
loglog(N, N.^(-4)*30, "r--", 'linewidth', 2);

# ------------------------  beautiful plot
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );

axis([3 30]);

hx=xlabel("N");
hy=ylabel("L_2 norm of error");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("img/l2.eps","-deps","-FArial", "-color");

# ---------------------------------------------------------------------------
#                                      H1
# ---------------------------------------------------------------------------
figure(1), clf, hold on;

# quadratic
loglog(N, h1_quad(1,:), 'b-o',"linewidth", 2, 'markersize', 10);
loglog(N, h1_quad(2,:), 'b-s',"linewidth", 2, 'markersize', 10);
loglog(N, h1_quad(3,:), 'b-^',"linewidth", 2,  'markersize', 10);
loglog(N, N.^(-2)*400, "r--", 'linewidth', 2);

# cubic
loglog(N, h1_cube(1,:), 'k-o',"linewidth", 2,  'markersize', 10);
loglog(N, h1_cube(2,:), 'k-s',"linewidth", 2,  'markersize', 10);
loglog(N, h1_cube(3,:), 'k-^',"linewidth", 2,  'markersize', 10);
loglog(N, N.^(-3)*200, "r--", 'linewidth', 2);

# ------------------------  beautiful plot
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );

axis([3 30]);

hx=xlabel("N");
hy=ylabel("H^1 norm of error");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("img/h1.eps","-deps","-FArial", "-color");
