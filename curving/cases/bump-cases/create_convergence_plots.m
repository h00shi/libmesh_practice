# this file is not a function
THIS_IS_NOT_A_FUNCTION=1;
clear THIS_IS_NOT_A_FUNCTION;

# function to append norms from a file
function [N,res] = append_file(file_name)

data = load(file_name);
N = data(:,1);
res = data(:,2);

endfunction

# ---------------------------------------------------------------------------
#                                     2nd order
# ---------------------------------------------------------------------------
# 1-> lagrange 2->hierarchic 3->bernstein
clf, hold on;

# lagrange
[N,res] = append_file('bump_l_o2.pp');
semilogy(N, res, 'b-;Lagrange;', 'linewidth', 2);

# hie
[N,res] = append_file('bump_h_o2.pp');
semilogy(N, res, 'r-;Hierarchic;', 'linewidth', 2);

# bern
[N,res] = append_file('bump_b_o2.pp');
semilogy(N, res, 'm-;Bernstein;', 'linewidth', 2);

# ------------------------  beautiful plot
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 13 , "location" , "northeast", "linewidth",2);

hx=xlabel("CG iteration");
hy=ylabel("LSE Residual");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("img/second-convergence.eps","-deps","-FArial", "-color");

# ---------------------------------------------------------------------------
#                                     3rd order
# ---------------------------------------------------------------------------
# 1-> lagrange 2->hierarchic 3->bernstein
clf, hold on;

# hie
[N,res] = append_file('bump_h_o3.pp');
semilogy(N, res, 'r-;Hierarchic;', 'linewidth', 2);

# bern
[N,res] = append_file('bump_b_o3.pp');
semilogy(N, res, 'm-;Bernstein;', 'linewidth', 2);

# ------------------------  beautiful plot
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );

copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 13 , "location" , "northeast", "linewidth",2);

hx=xlabel("CG iteration");
hy=ylabel("LSE Residual");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("img/third-convergence.eps","-deps","-FArial", "-color");
