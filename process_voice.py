# -*- coding: utf-8 -*-
"""
PROJETO DE PDS 2 - SIMULAÇÃO E GERAÇÃO DE GRÁFICOS EM PYTHON
Este script executa todas as simulaçoes de processamento de sinais do Projeto 2:
1. FFT de Cossenos (N=512 e N=8192) com janelas retangular e Hanning.
2. Projeto de filtros de 2ª ordem (LP, BP, HP) com r=0.97 e Thetas correspondentes.
3. Filtragem de voz por equacao de diferencas recursiva.
4. Filtragem por convolucao temporal truncada (Nc=64 e Nc=512) com janelas retangular e Hamming.
5. Filtragem por FFT (comprimento M correto) e comparacao de janelas.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.signal import lfilter

# Estilo visual moderno dos graficos
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['font.size'] = 9

# Caminhos do projeto
DIRETORIO_PROJETO = r"D:\Programação\PDS\Projeto2"
ARQUIVO_VOZ = os.path.join(DIRETORIO_PROJETO, "teste.wav")
PASTA_GRAFICOS = os.path.join(DIRETORIO_PROJETO, "Graficos")
PASTA_AUDIOS = os.path.join(DIRETORIO_PROJETO, "Audios_Processados")

# Cria pastas de saida
for pasta in [PASTA_GRAFICOS, PASTA_AUDIOS]:
    if not os.path.exists(pasta):
        os.makedirs(pasta)

# Lendo o sinal de voz
fs, sinal = wavfile.read(ARQUIVO_VOZ)
if sinal.dtype == np.int16:
    sinal = sinal / 32768.0
if len(sinal.shape) > 1:
    sinal = sinal[:, 0]

# Remove o componente DC (media do sinal) para evitar pico gigante em 0 Hz no espectro
sinal = sinal - np.mean(sinal)

P = len(sinal)
tempo = np.arange(P) / fs

print(f"Voz carregada. Fs={fs} Hz, P={P} amostras (sem offset DC).")

def salvar_audio(nome, dados):
    normalix = max(np.abs(np.min(dados)), np.abs(np.max(dados)))
    dados_norm = dados / normalix if normalix > 0 else dados
    caminho = os.path.join(PASTA_AUDIOS, nome)
    wavfile.write(caminho, fs, np.int16(dados_norm * 32767))
    return dados_norm

# =========================================================================
# 4.1 e 4.2) ANÁLISE ESPECTRAL DOS COSSENOS
# =========================================================================
print("Processando itens 4.1 e 4.2 (espectros dos cossenos)...")
Ns = [512, 8192]

for N in Ns:
    n = np.arange(N)
    x1 = np.cos(2 * np.pi * 3 / 10 * n)  # Usando freq discretas padrao do exemplo: cos(2*pi/3*n) -> 2*pi*3/10*n ou similar. 
    # Para o cosseno cos(2*pi/3 * n), a frequencia e 1/3 cycles/sample.
    x1_real = np.cos(2 * np.pi / 3 * n)
    x2_real = np.cos(np.pi / 4 * n)
    
    # Janela Hanning
    w_hn = np.hanning(N)
    
    # FFTs - Janela Retangular
    X1_ret = np.abs(np.fft.fft(x1_real))
    X2_ret = np.abs(np.fft.fft(x2_real))
    
    # FFTs - Janela Hanning
    X1_hn = np.abs(np.fft.fft(x1_real * w_hn))
    X2_hn = np.abs(np.fft.fft(x2_real * w_hn))
    
    # Plotagem
    plt.figure(figsize=(10, 8))
    
    # Subplot para x1
    plt.subplot(2, 1, 1)
    plt.plot(X1_ret[:N//2] / np.max(X1_ret), label="Retangular (Sem janela)", alpha=0.8)
    plt.plot(X1_hn[:N//2] / np.max(X1_hn), label="Janela Hanning", alpha=0.8, color='red')
    plt.title(f"Espectro de cos(2*pi/3 * n) - N = {N}")
    plt.xlabel("Bin de Frequência k")
    plt.ylabel("Magnitude Normalizada")
    plt.legend()
    
    # Subplot para x2
    plt.subplot(2, 1, 2)
    plt.plot(X2_ret[:N//2] / np.max(X2_ret), label="Retangular (Sem janela)", alpha=0.8)
    plt.plot(X2_hn[:N//2] / np.max(X2_hn), label="Janela Hanning", alpha=0.8, color='red')
    plt.title(f"Espectro de cos(pi/4 * n) - N = {N}")
    plt.xlabel("Bin de Frequência k")
    plt.ylabel("Magnitude Normalizada")
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_GRAFICOS, f"espectro_cossenos_N_{N}.png"), dpi=150)
    plt.close()

# =========================================================================
# 4.3) PROJETO DE FILTROS DE 2ª ORDEM E FILTRAGEM
# =========================================================================
print("Projetando filtros de 2ª ordem e filtrando o sinal de voz...")
r = 0.97

# Filtro LP (Theta = 0)
a_LP = -2 * r * np.cos(0)
b_LP = r**2

# Filtro BP (Theta = pi/2)
a_BP = -2 * r * np.cos(np.pi / 2)
b_BP = r**2

# Filtro HP (Theta = pi)
a_HP = -2 * r * np.cos(np.pi)
b_HP = r**2

# A) Filtragem recursiva por equacao de diferencas
print("Executando A) Filtragem por equação de diferenças...")
y_LP_rec = lfilter([1.0], [1.0, a_LP, b_LP], sinal)
y_BP_rec = lfilter([1.0], [1.0, a_BP, b_BP], sinal)
y_HP_rec = lfilter([1.0], [1.0, a_HP, b_HP], sinal)

y_LP_rec_norm = salvar_audio("output_4_3_A_LP.wav", y_LP_rec)
y_BP_rec_norm = salvar_audio("output_4_3_A_BP.wav", y_BP_rec)
y_HP_rec_norm = salvar_audio("output_4_3_A_HP.wav", y_HP_rec)

# Plota espectros da filtragem recursiva
plt.figure()
N_fft = 4096
f_axis = np.fft.fftfreq(N_fft, d=1/fs)[:N_fft//2]

# FFTs com janela de Hanning aplicada no sinal para analise espectral limpa
# Pegamos um trecho do meio do sinal (em t = 2s) onde a voz esta ativa, evitando o silencio inicial
n_start_voice = int(2.0 * fs)
trecho_sinal = sinal[n_start_voice : n_start_voice + N_fft] * np.hanning(N_fft)
trecho_LP = y_LP_rec[n_start_voice : n_start_voice + N_fft] * np.hanning(N_fft)
trecho_BP = y_BP_rec[n_start_voice : n_start_voice + N_fft] * np.hanning(N_fft)
trecho_HP = y_HP_rec[n_start_voice : n_start_voice + N_fft] * np.hanning(N_fft)

plt.plot(f_axis, np.abs(np.fft.fft(trecho_sinal))[:N_fft//2] / np.max(np.abs(np.fft.fft(trecho_sinal))), label="Original", alpha=0.5)
plt.plot(f_axis, np.abs(np.fft.fft(trecho_LP))[:N_fft//2] / np.max(np.abs(np.fft.fft(trecho_LP))), label="Passa-Baixas", alpha=0.8)
plt.plot(f_axis, np.abs(np.fft.fft(trecho_BP))[:N_fft//2] / np.max(np.abs(np.fft.fft(trecho_BP))), label="Passa-Bandas", alpha=0.8)
plt.plot(f_axis, np.abs(np.fft.fft(trecho_HP))[:N_fft//2] / np.max(np.abs(np.fft.fft(trecho_HP))), label="Passa-Altas", alpha=0.8)

plt.title("Espectro da Voz Filtrada por Equação de Diferenças (Item 4.3.A)")
plt.xlabel("Frequência (Hz)")
plt.ylabel("Magnitude Normalizada")
plt.xlim(0, 10000)  # Foca na faixa de 0 a 10 kHz (zoom solicitado)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PASTA_GRAFICOS, "espectros_filtragem_recursiva.png"), dpi=150)
plt.close()

# B) Filtragem por convolucao truncada (Nc = 64 e Nc = 512)
print("Executando B) Filtragem por convolução com resposta ao impulso truncada...")
Ncs = [64, 512]

for Nc in Ncs:
    n_filt = np.arange(Nc)
    
    # Respostas ao impulso teoricas
    h_LP = (n_filt + 1) * (r**n_filt)
    h_BP = (r**n_filt) * np.sin((n_filt + 1) * np.pi / 2)
    h_HP = (n_filt + 1) * ((-r)**n_filt)
    
    # Janela de Hamming para o filtro
    w_ham = np.hamming(Nc)
    
    # B.1) Janela Retangular (Filtros brutos)
    y_LP_ret = np.convolve(sinal, h_LP)[:P]
    y_BP_ret = np.convolve(sinal, h_BP)[:P]
    y_HP_ret = np.convolve(sinal, h_HP)[:P]
    
    salvar_audio(f"output_4_3_B1_ret_LP_Nc_{Nc}.wav", y_LP_ret)
    salvar_audio(f"output_4_3_B1_ret_BP_Nc_{Nc}.wav", y_BP_ret)
    salvar_audio(f"output_4_3_B1_ret_HP_Nc_{Nc}.wav", y_HP_ret)
    
    # B.2) Janela Hamming aplicada ao filtro
    y_LP_ham = np.convolve(sinal, h_LP * w_ham)[:P]
    y_BP_ham = np.convolve(sinal, h_BP * w_ham)[:P]
    y_HP_ham = np.convolve(sinal, h_HP * w_ham)[:P]
    
    salvar_audio(f"output_4_3_B2_ham_LP_Nc_{Nc}.wav", y_LP_ham)
    salvar_audio(f"output_4_3_B2_ham_BP_Nc_{Nc}.wav", y_BP_ham)
    salvar_audio(f"output_4_3_B2_ham_HP_Nc_{Nc}.wav", y_HP_ham)
    
    # Plota espectro comparando janelamento retangular vs Hamming para este Nc
    # Vamos demonstrar para o filtro Passa-Baixas (LP) como referencia
    plt.figure()
    # Pegamos o trecho com voz ativa a partir do segundo 2
    n_start_voice = int(2.0 * fs)
    trecho_LP_ret = y_LP_ret[n_start_voice : n_start_voice + N_fft]
    trecho_LP_ham = y_LP_ham[n_start_voice : n_start_voice + N_fft]
    
    # Analise espectral com janela retangular do sinal
    FFT_ret_sinal_ret = np.abs(np.fft.fft(trecho_LP_ret))[:N_fft//2]
    # Analise espectral com janela Hanning do sinal
    FFT_ret_sinal_hn = np.abs(np.fft.fft(trecho_LP_ret * np.hanning(N_fft)))[:N_fft//2]
    
    plt.plot(f_axis, FFT_ret_sinal_ret / np.max(FFT_ret_sinal_ret), label="Sinal com Janela Retangular", alpha=0.7)
    plt.plot(f_axis, FFT_ret_sinal_hn / np.max(FFT_ret_sinal_hn), label="Sinal com Janela Hanning", alpha=0.7, color='red')
    
    plt.title(f"Efeito da Janela de Análise Espectral na Voz Filtrada (LP, Nc={Nc})")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("Magnitude Normalizada")
    plt.xlim(0, 10000)  # Foca na faixa de 0 a 10 kHz (zoom solicitado)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PASTA_GRAFICOS, f"espectros_convolucao_{Nc}.png"), dpi=150)
    plt.close()

# C) Convolucao no dominio da frequencia usando FFT
print("Executando C) Convolução no domínio da frequência via FFT...")
# Com Nc = 512, comprimento minimo M = P + Nc - 1
comprimento_minimo = P + 512 - 1
M = 2**int(np.ceil(np.log2(comprimento_minimo)))
print(f"Tamanho P={P}, Nc=512, Comprimento FFT M={M}")

n_filt = np.arange(512)
h_LP = (n_filt + 1) * (r**n_filt)
h_LP_hn = h_LP * np.hanning(512)

# Padding e FFT
X = np.fft.fft(sinal, n=M)
H_ret = np.fft.fft(h_LP, n=M)
H_hn = np.fft.fft(h_LP_hn, n=M)

# Multiplicacao das FFTs
Y_ret = X * H_ret
Y_hn = X * H_hn

# Inversa da FFT
y_ret_fft = np.real(np.fft.ifft(Y_ret))[:P]
y_hn_fft = np.real(np.fft.ifft(Y_hn))[:P]

salvar_audio("output_4_3_C_LP_ret_Nc_512.wav", y_ret_fft)
salvar_audio("output_4_3_C_LP_hn_Nc_512.wav", y_hn_fft)

# Plota espectro do filtro bruto (respostas ao impulso sem janelamento)
plt.figure()
# FFT da resposta ao impulso truncada em Nc = 512 para LP, BP, HP
h_LP_512 = (np.arange(512) + 1) * (r**np.arange(512))
h_BP_512 = (r**np.arange(512)) * np.sin((np.arange(512) + 1) * np.pi / 2)
h_HP_512 = (np.arange(512) + 1) * ((-r)**np.arange(512))

H_LP_resp = np.abs(np.fft.fft(h_LP_512, n=N_fft))[:N_fft//2]
H_BP_resp = np.abs(np.fft.fft(h_BP_512, n=N_fft))[:N_fft//2]
H_HP_resp = np.abs(np.fft.fft(h_HP_512, n=N_fft))[:N_fft//2]

plt.plot(f_axis, H_LP_resp / np.max(H_LP_resp), label="Passa-Baixas", alpha=0.8)
plt.plot(f_axis, H_BP_resp / np.max(H_BP_resp), label="Passa-Bandas", alpha=0.8, color='green')
plt.plot(f_axis, H_HP_resp / np.max(H_HP_resp), label="Passa-Altas", alpha=0.8, color='red')

plt.title("Respostas em Frequência dos Filtros Projetados (Nc=512, Sem Janela)")
plt.xlabel("Frequência (Hz)")
plt.ylabel("Ganho Normalizado")
plt.xlim(0, fs / 2)  # Plota ate Nyquist para ver as bandas Passa-Banda (12 kHz) e Passa-Altas (24 kHz)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PASTA_GRAFICOS, "espectros_filtros_teoricos.png"), dpi=150)
plt.close()

print("Simulações do Projeto 2 concluídas! Todos os arquivos e gráficos foram gerados.")
