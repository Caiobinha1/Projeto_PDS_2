// =========================================================================
// PROJETO DE PDS 2 - ITEM 4.3.C E D: CONVOLUÇÃO NO DOMÍNIO DA FREQUÊNCIA (FFT)
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

// Comprimento do sinal de entrada (P)
P = length(sinal_voz);
disp("Comprimento do sinal de voz (P): " + string(P) + " amostras.");

// Comprimento do filtro truncado (Nc = 512)
Nc = 512;
disp("Comprimento da resposta ao impulso (Nc): " + string(Nc) + " amostras.");

// 2. DETERMINAR COMPRIMENTO M DA FFT
// Para evitar aliasing temporal na convolucao circular, precisamos de:
// M >= P + Nc - 1
comprimento_minimo = P + Nc - 1;
disp("Comprimento minimo para evitar aliasing: " + string(comprimento_minimo));

// Escolhemos a proxima potencia de 2 para otimizacao da FFT
potencia = ceil(log2(comprimento_minimo));
M = 2^potencia;
disp("Comprimento M escolhido para a FFT (potencia de 2): " + string(M));

// 3. PROJETO DOS FILTROS E CONVOLUÇÃO NO DOMÍNIO DA FREQUÊNCIA
r = 0.97;
n_filtro = 0:Nc-1;

// Gerando a resposta ao impulso do filtro Passa-Baixas (LP) como exemplo para demonstrar
// h[n] = (n+1) * r^n
h_LP = (n_filtro + 1) .* (r .^ n_filtro);

// Aplica a janela de Hanning ao filtro para reduzir ripple na resposta em frequencia
w_hn = window('hn', Nc); // Janela Hanning de tamanho Nc
h_LP_hn = h_LP .* w_hn'; // Resposta ao impulso truncada com janela Hanning

// PADDING (Preenchimento com zeros) ate o comprimento M
sinal_padded = [sinal_voz, zeros(1, M - P)];
h_padded_ret = [h_LP, zeros(1, M - Nc)];
h_padded_hn = [h_LP_hn, zeros(1, M - Nc)];

// 4. CÁLCULO DAS FFTS DE COMPRIMENTO M
X = fft(sinal_padded);
H_ret = fft(h_padded_ret);
H_hn = fft(h_padded_hn);

// Multiplicacao no dominio da frequencia
Y_ret = X .* H_ret;
Y_hn = X .* H_hn;

// Inversa da FFT para voltar ao dominio do tempo
y_ret_completo = real(ifft(Y_ret));
y_hn_completo = real(ifft(Y_hn));

// Cortamos para o tamanho original P do sinal de entrada
y_ret = y_ret_completo(1:P);
y_hn = y_hn_completo(1:P);

// 5. NORMALIZAÇÃO E EXPORTAÇÃO DOS ÁUDIOS (ITEM D)
normalix_ret = max(abs(min(y_ret)), abs(max(y_ret)));
y_ret_norm = y_ret / normalix_ret;

normalix_hn = max(abs(min(y_hn)), abs(max(y_hn)));
y_hn_norm = y_hn / normalix_hn;

// Salva os audios na pasta organizada
nome_saida_ret = "../Audios_Processados/output_4_3_C_LP_ret_Nc_512.wav";
nome_saida_hn = "../Audios_Processados/output_4_3_C_LP_hn_Nc_512.wav";

try
    wavwrite(y_ret_norm, Fs, nome_saida_ret);
    wavwrite(y_hn_norm, Fs, nome_saida_hn);
catch
    audiowrite(nome_saida_ret, y_ret_norm', Fs);
    audiowrite(nome_saida_hn, y_hn_norm', Fs);
end

disp("Audios filtrados via FFT salvos com sucesso!");
disp("Fim do script 4.3.C e D!");
