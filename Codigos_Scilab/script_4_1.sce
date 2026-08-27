// =========================================================================
// PROJETO DE PDS 2 - ITEM 4.1 E 4.2: ANÁLISE ESPECTRAL DE SINAIS COSSENOIDAIS
// =========================================================================

clc;
clear;

// Tamanhos de FFT sugeridos pelo professor
Ns = [512, 8192];

// Loop para passar pelos comprimentos N = 512 e N = 8192
for i = 1:length(Ns)
    N = Ns(i);
    disp("Processando para N = " + string(N) + "...");
    
    // Vetor de indices de tempo de 0 ate N-1
    n = 0:N-1;
    
    // Sinal 1: cos(2*pi/3 * n)
    x1 = cos((2 * %pi / 3) * n);
    
    // Sinal 2: cos(pi/4 * n)
    x2 = cos((%pi / 4) * n);
    
    // =========================================================================
    // ITEM 4.1: JANELA RETANGULAR (SINAL DIRETO, SEM MULTIPLICAR POR OUTRA JANELA)
    // =========================================================================
    X1_ret = abs(fft(x1));
    X2_ret = abs(fft(x2));
    
    // =========================================================================
    // ITEM 4.2: JANELA DE HANNING (APLICADA ANTES DA FFT)
    // Multiplicacao amostra a amostra ( .* ) pela janela hn do Scilab
    // =========================================================================
    w_hn = window('hn', N); // Gera a janela de Hanning de tamanho N
    x1_hn = x1 .* w_hn';     // w_hn e retornado como coluna, entao transpomos para linha
    x2_hn = x2 .* w_hn';
    
    X1_hn = abs(fft(x1_hn));
    X2_hn = abs(fft(x2_hn));
    
    // Eixo de frequencias discretas k (bins de 0 ate N-1)
    k = 0:N-1;
    
    // Gerar graficos comparativos para este comprimento N
    // Sinal 1: cos(2*pi/3 * n)
    scf(2*i - 1); // Abre ou seleciona a janela grafica
    clf();
    subplot(2, 1, 1);
    plot(k, X1_ret, "b");
    xtitle("Espectro cos(2*pi/3*n) - N = " + string(N) + " - Janela Retangular", "Bin de Frequencia k", "Magnitude");
    
    subplot(2, 1, 2);
    plot(k, X1_hn, "r");
    xtitle("Espectro cos(2*pi/3*n) - N = " + string(N) + " - Janela Hanning", "Bin de Frequencia k", "Magnitude");
    
    // Sinal 2: cos(pi/4 * n)
    scf(2*i);
    clf();
    subplot(2, 1, 1);
    plot(k, X2_ret, "b");
    xtitle("Espectro cos(pi/4*n) - N = " + string(N) + " - Janela Retangular", "Bin de Frequencia k", "Magnitude");
    
    subplot(2, 1, 2);
    plot(k, X2_hn, "r");
    xtitle("Espectro cos(pi/4*n) - N = " + string(N) + " - Janela Hanning", "Bin de Frequencia k", "Magnitude");
    
    disp("Graficos para N = " + string(N) + " gerados com sucesso!");
end

disp("Fim do script 4.1 e 4.2!");
