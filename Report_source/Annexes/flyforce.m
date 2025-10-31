% Ce script calcule et trace la masse totale de deux contrepoids 
% en fonction du parametre r1 en resolvant un systeme lineaire base 
% sur les parametres geometriques et inertiels du systeme.
% L'objectif est de determiner la valeur de r1 minimisant 
% la masse totale des contrepoids.

% Fonction de calcul des masses du contrepopids en fonction de r1
function [m1, m2] = calcul_masses(r1)
    % Parametres
    % Les longueurs sont donnees en mm, les masses en kg
    % et les inerties en kg*mm^2
    L  = 100;
    r2 = 57;
    r3 = 52;
    r4 = 45;
    Ia = 1625;
    ma = 0.217;
    Ib = 98;
    mb = 0.115;

    % Systeme d'equations lineaire a resoudre
    A = [ (L/2 - r1)/r1, -(L/2 + r1)/r1  ;
          (L/2 - r1)^2,   (L/2 + r1)^2  ];
    b = [ mb - (((r1 + r2)*r4)/(r1*r3))*ma                ;
          ((r1 + r2)/r3)*(Ia + ma*r4^2) - (Ib + mb*r1^2) ];

    % Resolution
    sol = A \ b;
    m1 = sol(1);
    m2 = sol(2);
end

% Calcule  des masses pour une plage de r
r_vals = linspace(-40, 40, 300);
M_vals = nan(size(r_vals));
for i = 1:length(r_vals)
    r = r_vals(i);
    [m1, m2] = calcul_masses(r);
    if m1 >= 0 && m2 >= 0 && isreal(m1) && isreal(m2)
        M_vals(i) = m1 + m2;
    end
end

% Trouve la valeur minimale
[~, idx] = min(M_vals);
r_opt = r_vals(idx);
M_opt = M_vals(idx);

% Trace le graphique
plot(r_vals, M_vals, 'b-', 'LineWidth', 2); hold on;

% Ligne verticale a r optimal
xline(r_opt, 'r--', 'Color', 'r');

% Point et etiquette pour M optimal
scatter(r_opt, M_opt, 60, 'r', 'filled');
text(r_opt + 1, M_opt + 0.4, sprintf('r1 = %.2f mm', r_opt), ...
     'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
text(r_opt + 1, M_opt + 0.2, sprintf('M = %.3f kg', M_opt), ...
     'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

% Etiquettes et titre
xlabel('r1 (mm)');
ylabel('M(r1) (kg)');
title('Masse totale des contrepoids en fonction de r1');
grid on;
legend('M(r1)', 'Minimum', 'Location', 'best');