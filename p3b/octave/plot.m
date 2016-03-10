N = [780, 3120, 12480];
sn = sqrt(N);

sol = load("-ascii", "out/false_dirich.txt");
lf = [sol(1,1), sol(2,1), sol(3,1)];
hf = [sol(1,2), sol(2,2), sol(3,2)];
plf = polyfit (log10(sn),  log10(lf),1);
phf = polyfit (log10(sn),  log10(hf),1);

sol = load("-ascii", "out/true_dirich.txt");
lt = [sol(1,1), sol(2,1), sol(3,1)];
ht = [sol(1,2), sol(2,2), sol(3,2)];
plt = polyfit (log10(sn),  log10(lt),1);
pht = polyfit (log10(sn),  log10(ht),1);

sol = load("-ascii", "out/true_ex.txt");
l3 = [sol(1,1), sol(2,1), sol(3,1)];
h3 = [sol(1,2), sol(2,2), sol(3,2)];
pl3 = polyfit (log10(sn),  log10(l3),1);
ph3 = polyfit (log10(sn),  log10(h3),1);

figure(1), clf, hold on;


loglog(sn, lf , ["-s;curved, slope=" ... 
	  sprintf("%.2f",-plf(1))  ";"],"linewidth", 2, "color", "k");

loglog(sn, l3 , ["--o;linear, correct cond., slope=" ... 
 	  sprintf("%.2f",-pl3(1))  ";"],"linewidth", 2, "color", "k");


loglog(sn, lt , ["-.^;linear, homogenous cond., slope=" ... 
	  sprintf("%.2f",-plt(1))  ";"],"linewidth", 2, "color", "k");

## axis
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );
axis([10 300]);

## legend
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 13 , "location" , "southwest", "linewidth",2);

##axis label and range
hx=xlabel(" sqrt(# Elems) ");
hy=ylabel("L_2 norm");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("out/l2.eps","-deps","-FArial");



figure(2), clf, hold on;

loglog(sn, hf , ["-s;curved, slope=" ... 
	  sprintf("%.2f",-phf(1))  ";"],"linewidth", 2, "color", "k");

loglog(sn, h3 , ["-.o;linear, correct cond., slope=" ... 
	  sprintf("%.2f",-ph3(1))  ";"],"linewidth", 2, "color", "k");
loglog(sn, ht , ["-.^;linear, homogenous cond., slope=" ... 
	  sprintf("%.2f",-pht(1))  ";"],"linewidth", 2, "color", "k");

## axis
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );
axis([10 300]);

## legend
copied_legend = findobj(gcf(),"type","axes","Tag","legend");
set(copied_legend, "FontSize", 13 , "location" , "southwest", "linewidth",2);

##axis label and range
hx=xlabel(" sqrt(# Elems) ");
hy=ylabel("H^1 norm");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print("out/h1.eps","-deps","-FArial");
