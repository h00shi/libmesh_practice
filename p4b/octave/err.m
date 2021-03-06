## plot error vs time
c = "kbrgmcy";

figure(1), clf, hold on;

err20=load(["../out/" pre "_err20"]);
semilogy(err20(:,1),err20(:,2), "color", c(1),"linewidth", 2);

err40=load(["../out/" pre "_err40"]);
semilogy(err40(:,1),err40(:,2), "color", c(2),"linewidth", 2);

err80=load(["../out/" pre "_err80"]);
semilogy(err80(:,1),err80(:,2), "color", c(3),"linewidth", 2);

err160=load(["../out/" pre "_err160"]);
semilogy(err160(:,1),err160(:,2), "color", c(4),"linewidth", 2);

err320=load(["../out/" pre "_err320"]);
semilogy(err320(:,1),err320(:,2), "color", c(5),"linewidth", 2);

err640=load(["../out/" pre "_err640"]);
semilogy(err640(:,1),err640(:,2), "color", c(6),"linewidth", 2);

err1280=load(["../out/" pre "_err1280"]);
semilogy(err1280(:,1),err1280(:,2), "color", c(7),"linewidth", 2);

err2560=load(["../out/" pre "_err2560"]);
semilogy(err2560(:,1),err2560(:,2), "color", c(1),"linewidth", 2);

err5120=load(["../out/" pre "_err5120"]);
semilogy(err5120(:,1),err5120(:,2), "color", c(2),"linewidth", 2);

##axis label and range
hx=xlabel("Time");
hy=ylabel("|u-u_tp|_L2");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print(["../plt/" pre "_1d_err.eps"],"-deps","-FArial", "-color");

#-------------------------------------------------------------
# refinement
#-------------------------------------------------------------

 nt = [20 40 80 160 320 640 1280 2560 5120];

 err_ = [max(err20(:,2)) max(err40(:,2)) max(err80(:,2)) ...
            max(err160(:,2)) max(err320(:,2)) ...
			max(err640(:,2))...
            max(err1280(:,2))  max(err2560(:,2)) ...
			max(err5120(:,2)) ];

p = polyfit (log10(nt),  log10(err_),1);

figure(2), clf, hold on;
loglog(nt, err_, "bo-", "linewidth", 2);
nn = [20 520];
ee = nn .^ (-2) * 100;
loglog(nn, ee, "k--", 'linewidth', 2);
printf("slope: %.3f \n", -p(1));

## axis
xt = get(gca,"XTick");
set(gca,"XTickLabel", sprintf("%.0f|",xt) );


##axis label and range
hx=xlabel("# Time Step");
hy=ylabel("maximum steady state error");
set(hx, "fontsize", 15, "linewidth", 2);
set(hy, "fontsize", 15, "linewidth", 2);
set(gca,  "fontsize", 15);

print(["../plt/" pre "_1d_ref.eps"],"-deps","-FArial", "-color");
