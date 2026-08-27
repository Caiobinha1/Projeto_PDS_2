// =========================================================================
// PROJETO DE PDS 2 - ITEM 4.3.A: FILTRAGEM POR EQUAÇÃO DE DIFERENÇAS RECURSIVA
// =========================================================================

clc;
clear;

// 1. CARREGAR AUDIO DE VOZ
arquivo_entrada = "../teste.wav";
try
    [x_voz, Fs] = wavread(arquivo_entrada);
catch
    disp("Erro ao usar wavread. Tentando usar a funcao audioread...");
    [x_voz, Fs] = audioread(arquivo_entrada);
end

// Garante sinal mono
[linhas, colunas] = size(x_voz);
if linhas > 1 then
    sinal_voz = x_voz(1, :);
else
    sinal_voz = x_voz';
end
N_voz = length(sinal_voz);

// 2. PARÂMETROS DO FILTRO DE 2ª ORDEM
r = 0.97;

// Filtro 1: Passa-Baixas (LP) - Theta = 0
theta_LP = 0;
a_LP = -2 * r * cos(theta_LP);
b_LP = r^2;

// Filtro 2: Passa-Bandas (BP) - Theta = pi/2
theta_BP = %pi / 2;
a_BP = -2 * r * cos(theta_BP);
b_BP = r^2;

// Filtro 3: Passa-Altas (HP) - Theta = pi
theta_HP = %pi;
a_HP = -2 * r * cos(theta_HP);
b_HP = r^2;

// Vetores de coeficientes para o loop
as = [a_LP, a_BP, a_HP];
bs = [b_LP, b_BP, b_HP];
nomes_filtros = ["Passa-Baixas", "Passa-Bandas", "Passa-Altas"];
prefixos = ["LP", "BP", "HP"];

// Loop para aplicar cada um dos 3 filtros no sinal de voz
for f_idx = 1:3
    a = as(f_idx);
    b = bs(f_idx);
    nome = nomes_filtros(f_idx);
    prefixo = prefixos(f_idx);
    
    disp("Aplicando filtro " + nome + " usando equacao de diferencas recursiva...");
    
    // Vetor de saida
    y = zeros(1, N_voz);
    
    // Condicoes iniciais manuais (para n=1 e n=2)
    y(1) = sinal_voz(1);
    y(2) = sinal_voz(2) - a * y(1);
    
    // Equacao recursiva: y[n] = x[n] - a*y[n-1] - b*y[n-2]
    for n = 3:N_voz
        y(n) = sinal_voz(n) - a * y(n-1) - b * y(n-2);
    end
    
    // Normalizacao
    normalix = max(abs(min(y)), abs(max(y)));
    y_norm = y / normalix;
    
    // Salva o audio correspondente
    nome_saida = "../Audios_Processados/output_4_3_A_" + prefixo + ".wav";
    try
        wavwrite(y_norm, Fs, nome_saida);
    catch
        audiowrite(nome_saida, y_norm', Fs);
    end
    
    disp("Audio filtrado salvo: " + nome_saida);
end

disp("Fim do script 4.3.A!");
