clear; clc; close all


%% LOAD DATA
% Image size
Isz = [112 112 112];
Nfr = 22;

% Load HARPI and SinMod strain
load('motion3D_HARPI.mat','st');
HRR = NaN([Isz Nfr]);
HCC = NaN([Isz Nfr]);
HLL = NaN([Isz Nfr]);
HRR(repmat(st.maskimage,[1 1 1 Nfr])) = st.RR(:);
HCC(repmat(st.maskimage,[1 1 1 Nfr])) = st.CC(:);
HLL(repmat(st.maskimage,[1 1 1 Nfr])) = st.LL(:);



%% STRAIN CURVES

% Plot settings
visible = 'off';
api = struct(...
    'AxesFontSize',  18,...
    'AxesLineWidth', 2,...
    'LegendFontSize', 17,...
    'Axis', [0 23 -0.3 0.3],...
    'XLabel', false,...
    'YLabel', false,...
    'XLabelStr', 'Frame',...
    'YLabelStr', 'Strain',...
    'YAxisTickValues', []);
plot_line_width = 2;
plot_marker_size = 3.5;

% Colors
co = [1 0 0;
      0 1 0;
      0 0 1];

% HARPI
HGLS = squeeze(mean(HLL,[1 2 3],'omitnan'));
HGRS = squeeze(mean(HRR,[1 2 3],'omitnan'));
HGCS = squeeze(mean(HCC,[1 2 3],'omitnan'));
figure,
plot(HGCS, 'LineWidth', 2, 'Color', co(1, :)); hold on;
plot(HGRS, 'LineWidth', 2, 'Color', co(2, :)); hold on;
plot(HGLS, 'LineWidth', 2, 'Color', co(3, :)); hold off;
x_range = [1 length(HGCS)];
xlim(x_range);
ylim([-0.3 0.3]); 
ax = gca;
xlabel(api.XLabelStr, 'interpreter', 'LaTeX');
ax.XAxis.TickValues = linspace(x_range(1), x_range(2), 3); % Dividir en 3 ticks
ax.YAxis.TickValues = -0.3:0.1:0.3; % Valores fijos para el eje Y
ylabel(api.YLabelStr, 'interpreter', 'LaTeX');
ax.Box = 'on';
ax.FontWeight = 'bold';
ax.FontSmoothing = 'on';
ax.FontSize = api.AxesFontSize;
ax.LineWidth = api.AxesLineWidth;
ax.TickLength = [0.025, 0.25];
ax.XAxis.TickLength = [0.025, 0.25];
ax.TickDir = 'in';
ax.TickLabelInterpreter = 'latex';


l = legend('GCS', 'GRS', 'GLS');
l.FontSize = api.LegendFontSize;
l.Interpreter = 'latex';


set(gca, 'Position', [0.1864 0.1590 0.7186 0.7660]);
set(gcf, 'Position', [680 544 440 434]);


print('-depsc', '-r600', 'curve_harpi3D');
