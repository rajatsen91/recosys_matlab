function plotDOCPerf(keyword)

%fname1 = 'perfDOCUnknown_n1000_K5_test2_p6.mat';
%fname2 = 'perfDOCUnknown_n1000_K5_test2_p7.mat';
%fname3 = 'perfDOCUnknown_n1000_K5_test2_p8.mat';
fname1 = 'perfDOCUnknown_n1000_K5_test3_p6.mat';
fname2 = 'perfDOCUnknown_n1000_K5_test3_p7.mat';

%fname4 = 'perfDOC_n1000_K5_test_p6.mat';
%fname5 = 'perfDOC_n1000_K5_test_p7.mat'; 
%fname6 = 'perfDOC_n1000_K5_test_p8.mat';

fname4 = 'perfDOC_n1000_K5_test3_p6.mat';
fname5 = 'perfDOC_n1000_K5_test3_p7.mat'; 
%fname6 = 'perfDOC_n1000_K5_test3_p8.mat';

fname7 = 'perfLinkPart_n1000_K5_p7q05.mat';
fname8 = 'perfLinkPart_n1000_K5_p7q1.mat';
fname9 = 'perfLinkPart_n1000_K5_p7q15.mat';

result_fname = strcat('perf_DOC_n1000_K5_', keyword);
n = 1000;
K = 5;

close all hidden;
f = figure;

%--------------------------------------------------------------------------
h1 = load(fname1);
res = h1.res;
clear h1;
numq = length(res.loop);
parr = zeros(1,numq);
qarr = zeros(1,numq);
err = zeros(1,numq);
for i = 1:numq
    parr(i) = res.loop{i}.p;
    qarr(i) = res.loop{i}.q;
    err(i) = mean(res.loop{i}.error);
end
%plot(qarr, 100*err, 'r-^', 'DisplayName', 'p = .6 (estimated)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%plot((parr-qarr)./sqrt(parr), 100*err, 'r-^', 'DisplayName', 'p = .6 (estimated)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
plot((parr-qarr), 100*err, 'r-s', 'DisplayName', 'DOC, p = .6 (est)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
h2 = load(fname2);
res = h2.res;
clear h2;
numq = length(res.loop);
parr = zeros(1,numq);
qarr = zeros(1,numq);
err = zeros(1,numq);
for i = 1:numq
    parr(i) = res.loop{i}.p;
    qarr(i) = res.loop{i}.q;
    err(i) = mean(res.loop{i}.error);
end
%plot(qarr, 100*err, 'b-o', 'DisplayName', 'p = .7 (estimated)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%plot((parr-qarr)./sqrt(parr), 100*err, 'r-o', 'DisplayName', 'p = .7 (estimated)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
plot((parr-qarr), 100*err, 'r-o', 'DisplayName', 'DOC, p = .7 (est)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% h3 = load(fname3);
% res = h3.res;
% clear h3;
% numq = length(res.loop);
% parr = zeros(1,numq);
% qarr = zeros(1,numq);
% err = zeros(1,numq);
% for i = 1:numq
%     parr(i) = res.loop{i}.p;
%     qarr(i) = res.loop{i}.q;
%     err(i) = mean(res.loop{i}.error);
% end
% plot(qarr, 100*err, 'r-^', 'DisplayName', 'p = .8 (estimated)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
h4 = load(fname4);
res = h4.res;
clear h4;
numq = length(res.loop);
parr = zeros(1,numq);
qarr = zeros(1,numq);
err = zeros(1,numq);
for i = 1:numq
    parr(i) = res.loop{i}.p;
    qarr(i) = res.loop{i}.q;
    err(i) = mean(res.loop{i}.error);
end
%plot(qarr, 100*err, 'r-+', 'DisplayName', 'p = .6 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%plot((parr-qarr)./sqrt(parr), 100*err, 'b-+', 'DisplayName', 'p = .6 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
plot((parr-qarr), 100*err, 'b-+', 'DisplayName', 'DOC, p = .6', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
h5 = load(fname5);
res = h5.res;
clear h5;
numq = length(res.loop);
parr = zeros(1,numq);
qarr = zeros(1,numq);
err = zeros(1,numq);
for i = 1:numq
    parr(i) = res.loop{i}.p;
    qarr(i) = res.loop{i}.q;
    err(i) = mean(res.loop{i}.error);
end
%plot(qarr, 100*err, 'b->', 'DisplayName', 'p = .7 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%plot((parr-qarr)./sqrt(parr), 100*err, 'b->', 'DisplayName', 'p = .7 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
plot((parr-qarr), 100*err, 'b->', 'DisplayName', 'DOC, p = .7', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% h6 = load(fname6);
% res = h6.res;
% clear h6;
% numq = length(res.loop);
% parr = zeros(1,numq);
% qarr = zeros(1,numq);
% err = zeros(1,numq);
% for i = 1:numq
%     parr(i) = res.loop{i}.p;
%     qarr(i) = res.loop{i}.q;
%     err(i) = mean(res.loop{i}.error);
% end
% plot(qarr, 100*err, 'b-*', 'DisplayName', 'p = .8 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
% Link Partition
h7 = load(fname7);
res = h7.res;
clear h7;
h8 = load(fname8);
res8 = h8.res;
clear h8;
h9 = load(fname9);
res9 = h9.res;
clear h9;
res.loop{2} = res8.loop{1};
res.loop{3} = res9.loop{1};

numq = length(res.loop);
parr = zeros(1,numq);
qarr = zeros(1,numq);
err = zeros(1,numq);
for i = 1:numq
    parr(i) = res.loop{i}.p;
    qarr(i) = res.loop{i}.q;
    err(i) = mean(res.loop{i}.error);
end
%plot(qarr, 100*err, 'b->', 'DisplayName', 'p = .7 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%plot((parr-qarr)./sqrt(parr), 100*err, 'b->', 'DisplayName', 'p = .7 (actual)', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
plot((parr-qarr), 100*err, 'g-*', 'DisplayName', 'LP, p = .7', 'LineWidth',2,'MarkerSize',8); grid on; hold on;
%--------------------------------------------------------------------------
%xlabel('q --->');
%xlabel('(p-q)/sqrt(p) --->');
xlabel('(p-q) --->');
ylabel('Average percentage error ->');
title_str = 'Error performance in graph with n = 1000, K = 5';
%title(title_str);
title(title_str,'FontSize',12);
%l = legend('IG, p(high) = .9', 'IG, p(high) = .85', 'LDA, p(high) = .9', 'LDA, p(high) = .85', 'Location', 'SouthEast');
%l = legend('show','Location', 'Best');
l = legend('show','Location', 'NorthWest');
set(l,'FontSize',14);
saveas(f,strcat(result_fname,'.jpg'));
saveas(f,strcat(result_fname,'.fig'));
