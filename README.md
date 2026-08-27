# Projeto 2: Análise Espectral e Filtragem de Sinais

Este repositório contém o código-fonte, dados e resultados do Projeto 2 da disciplina de Processamento Digital de Sinais (EEL7522) da Universidade Federal de Santa Catarina (UFSC).

## Integrantes
- Caio Missfeld Carlos (Matrícula: 22202674)
- Miguel Sória da Luz (Matrícula: 22100860)

## Estrutura do Projeto

O projeto está organizado na seguinte estrutura:

- **input_voice.wav**: Áudio original gravado (5.9 segundos, 48 kHz, mono) usado para processamento.
- **process_voice.py**: Script Python para processamento de sinais, que calcula os espectros dos cossenos, projeta os filtros de 2ª ordem, realiza a convolução temporal/frequencial na voz e salva os gráficos.
- **relatorio_pds.pdf**: Relatório final do projeto contendo as análises dos resultados e gráficos de espectro.
- **Codigos_Scilab/**: Pasta com os scripts Scilab (.sce) correspondentes aos itens solicitados:
  - **script_4_1.sce**: FFT de cossenos com janelas retangular e Hanning para N=512 e N=8192.
  - **script_4_3_A.sce**: Filtragem recursiva de voz por equação de diferenças.
  - **script_4_3_B.sce**: Filtragem por convolução temporal com resposta ao impulso truncada em Nc=64 e Nc=512 com janelas retangular e Hamming.
  - **script_4_3_C_D.sce**: Convolução no domínio da frequência via FFT e IFFT.
- **Audios_Processados/**: Pasta contendo todos os áudios resultantes da filtragem no formato WAV.
- **Graficos/**: Pasta com as imagens dos espectros de frequências gerados pelo script.

## Funcionamento

O projeto aborda dois temas principais em processamento digital de sinais:

1. **Análise Espectral e Vazamento:** Avaliação da transformada de Fourier (FFT) em cossenos. O cosseno em $\pi/4$ cai exatamente em uma frequência múltipla da amostragem (bin inteiro), não sofrendo vazamento. O cosseno em $2\pi/3$ não coincide com bins inteiros, gerando espalhamento espectral (vazamento). A janela de Hanning reduz significativamente os lobos secundários do vazamento ao custo de um leve alargamento da banda de pico.
2. **Projeto de Ressoadores de 2ª Ordem:** Projetamos filtros passa-baixas ($\Theta=0$), passa-bandas ($\Theta=\pi/2$) e passa-altas ($\Theta=\pi$) com raio $r=0.97$.
   - **Tempo (Recursivo):** Filtro IIR implementado diretamente pela equação de diferenças recursiva.
   - **Tempo (Filtro FIR Truncado):** Resposta ao impulso truncada em 64 e 512 amostras com janela retangular e Hamming. A janela de Hamming reduz ondulações e o truncamento longo (512) provê maior seletividade.
   - **Frequência (FFT/IFFT):** Multiplicação das FFTs das sequências. O comprimento M necessário para a FFT é determinado por $M \ge P + Nc - 1$, arredondado para a próxima potência de 2 para otimização rápida de processamento, garantindo que não ocorra aliasing temporal na convolução linear resultante.
