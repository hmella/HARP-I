% COMPARE_PRINCIPAL_VS_ROTATED
%
% demo_2D.m historically assigned the principal strains (st.p1, st.p2) to
% variables named RR and CC, while strain_curves2D.m reads the genuine
% rotated components (st.RR, st.CC) from the same saved file. The two are
% different quantities. This script quantifies the difference so the choice
% of what to display is made on evidence rather than by accident.
%
% Run demo_2D.m first (it writes motion2D_HARPI.mat).

clear; clc; close all

load('motion2D_HARPI.mat','st')
Nfr = size(st.CC,2);

gCC = mean(st.CC,1,'omitnan');
gRR = mean(st.RR,1,'omitnan');
gp1 = mean(st.p1,1,'omitnan');
gp2 = mean(st.p2,1,'omitnan');

fprintf('\n           circumferential          radial\n');
fprintf('fr     CC        p2      d       RR       p1       d\n');
for i = 1:Nfr
    fprintf('%2d  %8.4f %8.4f %7.4f  %8.4f %8.4f %7.4f\n', ...
        i, gCC(i), gp2(i), gp2(i)-gCC(i), gRR(i), gp1(i), gp1(i)-gRR(i));
end

fprintf('\nmax |p2 - CC| over frames : %.4f  (%.1f%% of peak |CC|)\n', ...
    max(abs(gp2-gCC)), 100*max(abs(gp2-gCC))/max(abs(gCC)));
fprintf('max |p1 - RR| over frames : %.4f  (%.1f%% of peak |RR|)\n', ...
    max(abs(gp1-gRR)), 100*max(abs(gp1-gRR))/max(abs(gRR)));

% pixelwise, not just the spatial mean
dCC = st.p2(:) - st.CC(:);
dRR = st.p1(:) - st.RR(:);
fprintf('\npixelwise max |p2 - CC|    : %.4f\n', max(abs(dCC),[],'omitnan'));
fprintf('pixelwise max |p1 - RR|    : %.4f\n', max(abs(dRR),[],'omitnan'));

figure('Position',[100 100 900 380]);
subplot(1,2,1); hold on; grid on
plot(1:Nfr, gCC, '-o', 'LineWidth',1.5, 'DisplayName','E_{cc} (rotated)');
plot(1:Nfr, gp2, '--s','LineWidth',1.5, 'DisplayName','2nd principal');
xlabel('frame'); ylabel('strain'); title('Circumferential'); legend('Location','best')

subplot(1,2,2); hold on; grid on
plot(1:Nfr, gRR, '-o', 'LineWidth',1.5, 'DisplayName','E_{rr} (rotated)');
plot(1:Nfr, gp1, '--s','LineWidth',1.5, 'DisplayName','1st principal');
xlabel('frame'); ylabel('strain'); title('Radial'); legend('Location','best')
