// =========================================================================
// PROJETO DE PDS 2 - ITEM 4.3.B: FILTRAGEM POR CONVOLUÇÃO TEMPORAL TRUNCADA
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

[linhas, colunas] = size(x_voz);
if linhas > 1 then
    sinal_voz = x_voz(1, :);
else
    sinal_voz = x_voz';
end
N_voz = length(sinal_voz);

// Parâmetro do filtro de 2ª ordem
r = 0.97;

// Comprimentos de truncamento (Nc = 64 e Nc = 512)
Ncs = [64, 512];

// Loop pelos comprimentos de truncamento
for i = 1:length(Ncs)
    Nc = Ncs(i);
    disp("Gerando respostas ao impulso truncadas com Nc = " + string(Nc) + "...");
    
    // Vetor de indices de tempo de 0 ate Nc-1 para o filtro
    n_filtro = 0:Nc-1;
    
    // Resposta ao impulso teorica para cada filtro:
    // Passa-Baixas (Theta = 0): h[n] = (n+1) * r^n
    h_LP = (n_filtro + 1) .* (r .^ n_filtro);
    
    // Passa-Bandas (Theta = pi/2): h[n] = r^n * sin((n+1)*pi/2)
    h_BP = (r .^ n_filtro) .* sin((n_filtro + 1) * %pi / 2);
    
    // Passa-Altas (Theta = pi): h[n] = (n+1) * (-r)^n
    h_HP = (n_filtro + 1) .* ((-r) .^ n_filtro);
    
    // Lista de filtros para o loop
    filtros = [h_LP; h_BP; h_HP];
    prefixos = ["LP", "BP", "HP"];
    nomes = ["Passa-Baixas", "Passa-Bandas", "Passa-Altas"];
    
    // B.1) Janela Retangular (direto da equacao, sem multiplicar por outra janela)
    for f_idx = 1:3
        h_atual = filtros(f_idx, :);
        prefixo = prefixos(f_idx);
        nome = nomes(f_idx);
        
        // Faz a convolucao do sinal de voz com o filtro
        y_completo = convol(h_atual, sinal_voz);
        y = y_completo(1:N_voz); // Mantem o tamanho
        
        // Normaliza e salva
        normalix = max(abs(min(y)), abs(max(y)));
        y_norm = y / normalix;
        nome_saida = "../Audios_Processados/output_4_3_B1_ret_" + prefixo + "_Nc_" + string(Nc) + ".wav";
        try
            wavwrite(y_norm, Fs, nome_saida);
        catch
            audiowrite(nome_saida, y_norm', Fs);
        end
        disp("Salvo (Janela Retangular): " + nome_saida);
    end
    
    // B.2) Janela Hamming (multiplicando a resposta h pela janela 'hm' do Scilab)
    w_ham = window('hm', Nc); // Janela de Hamming de tamanho Nc
    for f_idx = 1:3
        h_atual = filtros(f_idx, :);
        h_windowed = h_atual .* w_ham'; // Multiplica amostra a amostra
        prefixo = prefixos(f_idx);
        nome = nomes(f_idx);
        
        // Convolucao
        y_completo = convol(h_windowed, sinal_voz);
        y = y_completo(1:N_voz);
        
        // Normaliza e salva
        normalix = max(abs(min(y)), abs(max(y)));
        y_norm = y / normalix;
        nome_saida = "../Audios_Processados/output_4_3_B2_ham_" + prefixo + "_Nc_" + string(Nc) + ".wav";
        try
            wavwrite(y_norm, Fs, nome_saida);
        catch
            audiowrite(nome_saida, y_norm', Fs);
        end
        disp("Salvo (Janela Hamming): " + nome_saida);
    end
end

disp("Fim do script 4.3.B!");
