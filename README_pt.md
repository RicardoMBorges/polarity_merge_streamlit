# Polarity Merge

**Integração de redes moleculares entre os modos de ionização positivo e negativo por LC–MS/MS**  
**LAABio — Laboratório de Análise e Avaliação da Biodiversidade | IPPN–UFRJ**  
**Desenvolvido por Ricardo M. Borges**

O **Polarity Merge** é uma aplicação [Streamlit](https://streamlit.io/) que conecta redes moleculares adquiridas nos modos de ionização **positivo (POS)** e **negativo (NEG)**, utilizando concordância de **massa neutra** e evidências de **tempo de retenção (RT)**. O programa preserva as redes originais e acrescenta conexões entre polaridades, exportáveis para o [Cytoscape](https://cytoscape.org/).

> **Versão documentada: v0.7.** A hipótese de pareamento atualmente implementada é **`[M+H]+ ↔ [M−H]−`**. As associações produzidas são **candidatas**, não identificações químicas definitivas.

## Sumário

1. [Objetivo e princípio](#objetivo-e-princípio)
2. [Principais funcionalidades](#principais-funcionalidades)
3. [Instalação e publicação](#instalação-e-publicação)
4. [Arquivos de entrada](#arquivos-de-entrada)
5. [Tutorial passo a passo](#tutorial-passo-a-passo)
6. [Fundamentos do pareamento por massa](#fundamentos-do-pareamento-por-massa)
7. [Distribuição de ΔRT e filtragem por MAD](#distribuição-de-δrt-e-filtragem-por-mad)
8. [Correlação de abundâncias](#correlação-de-abundâncias)
9. [Visualização dos espectros MS/MS](#visualização-dos-espectros-msms)
10. [Exportação e Cytoscape](#exportação-e-cytoscape)
11. [Correspondências múltiplas e ambiguidades](#correspondências-múltiplas-e-ambiguidades)
12. [Exemplo de interpretação](#exemplo-de-interpretação)
13. [Limitações e controle de qualidade](#limitações-e-controle-de-qualidade)
14. [Solução de problemas](#solução-de-problemas)
15. [Reprodutibilidade e citação](#reprodutibilidade-e-citação)

## Objetivo e princípio

A ionização por electrospray (ESI) nos modos positivo e negativo pode revelar conjuntos complementares de metabólitos. Como as redes moleculares são frequentemente construídas separadamente, seus identificadores de nós e escores de similaridade não constituem, por si só, uma correspondência entre polaridades.

O Polarity Merge procura responder:

> **Um feature detectado em POS e outro detectado em NEG podem representar a mesma molécula neutra?**

A hipótese é avaliada a partir da massa neutra reconstruída e da compatibilidade cromatográfica. O aplicativo **não funde os nós**, **não impõe correspondência um-para-um** e **não interpreta a similaridade de fragmentação entre polaridades como equivalente ao cosseno intrapolaridade**.

## Principais funcionalidades

- Importação de duas redes moleculares independentes em **GraphML**.
- Busca de pares POS↔NEG com base em `[M+H]+ / [M−H]−`.
- Tolerância de massa configurável em **ppm** ou **Da**.
- Janela inicial de RT para geração de candidatos.
- Cálculo de `ΔRT = RT_NEG − RT_POS`, com sinal, em minutos e segundos.
- Estimativa robusta do deslocamento cromatográfico por **mediana** e **MAD**.
- Filtragem por **mediana ± k × MAD**, com `k` ajustável (padrão: 3).
- Gráficos de concordância de RT e distribuição de ΔRT, com referências em ±2 e ±3 MAD.
- Tabela editável de candidatos.
- Correlação opcional de abundâncias entre amostras correspondentes.
- Visualização lado a lado dos espectros MS/MS de POS e NEG.
- Exportação de **GraphML**, CSV de conexões aceitas e CSV de auditoria de RT.
- Preservação dos nós, das conexões intrapolaridade e de seus atributos.
- Identidade visual do LAABio na barra lateral.
- Execução explícita pelo botão **Run Polarity Merge**.

## Instalação e publicação

### Requisitos

Recomenda-se Python 3.10 ou superior. As dependências estão em `requirements.txt`: `streamlit`, `pandas`, `numpy`, `networkx`, `plotly` e `openpyxl`.

### Execução local

```bash
git clone https://github.com/RicardoMBorges/polarity_merge_streamlit.git
cd polarity_merge_streamlit

python -m venv .venv

# Windows PowerShell
.venv\Scripts\Activate.ps1

# Linux / macOS
# source .venv/bin/activate

pip install -r requirements.txt
streamlit run app.py
```

O Streamlit informará o endereço local da aplicação, normalmente `http://localhost:8501`.

### Publicação no Streamlit Community Cloud

1. Envie ao repositório GitHub os arquivos `app.py`, `requirements.txt`, `README.md`, `README_pt.md` e a pasta `assets/` com o logotipo.
2. Acesse [share.streamlit.io](https://share.streamlit.io/).
3. Crie uma aplicação e conecte o repositório.
4. Informe **`app.py`** como arquivo principal.
5. Publique e abra a URL disponibilizada.

A versão atual não exige banco de dados nem credenciais próprias. Para dados confidenciais ou não publicados, considere as políticas de acesso e armazenamento da infraestrutura utilizada.

### Estrutura recomendada

```text
polarity_merge_streamlit/
├── app.py
├── requirements.txt
├── README.md
├── README_pt.md
├── .gitignore
└── assets/
    └── LAABio_logo.png
```

Os botões de tutorial da aplicação podem apontar diretamente para:

- **English tutorial:** https://github.com/RicardoMBorges/polarity_merge_streamlit/blob/main/README.md
- **Tutorial em português:** https://github.com/RicardoMBorges/polarity_merge_streamlit/blob/main/README_pt.md

## Arquivos de entrada

### Obrigatórios: duas redes GraphML

| Arquivo | Função |
|---|---|
| **Positive GraphML** | Rede molecular do modo positivo |
| **Negative GraphML** | Rede molecular do modo negativo |

Para participar do pareamento, cada nó precisa fornecer:

- **m/z:** atributo `mz` ou equivalente reconhecido pelo aplicativo.
- **RT:** atributo `rt_min` ou `rt`, em **minutos**.
- **Identificador:** ID do nó no GraphML.

O aplicativo mantém os atributos originais, incluindo, quando disponíveis, `score`, `EdgeScore`, `matched_peaks`, `EdgeType`, `deltamz` e anotações.

**Atenção às unidades:** arquivos MGF frequentemente registram `RTINSECONDS`, enquanto o pareamento utiliza o RT do GraphML em **minutos**.

### Opcionais: espectros MGF

| Arquivo | Função |
|---|---|
| **Positive MGF** | Espectros MS/MS em POS |
| **Negative MGF** | Espectros MS/MS em NEG |

O visualizador procura os espectros pelos identificadores presentes nos metadados MGF, incluindo `FEATURE_ID` e `SCANS`. Para que o espectro seja recuperado, o ID do nó deve corresponder ao identificador indexado no MGF.

### Opcionais: tabelas quantitativas

São necessárias duas tabelas feature × amostra para calcular correlação entre polaridades. Cada tabela deve conter uma coluna de ID e colunas numéricas de abundância por amostra.

**Exemplo POS (valores ilustrativos):**

| FeatureID | Amostra_A_POS | Amostra_B_POS | Amostra_C_POS |
|---|---:|---:|---:|
| 1721 | 150000 | 220000 | 180000 |
| 5427 | 80000 | 65000 | 99000 |

**Exemplo NEG (valores ilustrativos):**

| FeatureID | Amostra_A_NEG | Amostra_B_NEG | Amostra_C_NEG |
|---|---:|---:|---:|
| 85 | 140000 | 210000 | 170000 |
| 186 | 76000 | 62000 | 95000 |

**Importante:** GraphML e MGF, isoladamente, não fornecem perfis completos de intensidade feature × amostra. O campo `FEATURE_MS1_HEIGHT` de um espectro MGF não substitui uma matriz quantitativa.

## Tutorial passo a passo

### Etapa 1 — Carregar as redes

Na barra lateral, carregue os arquivos **Positive GraphML** e **Negative GraphML**. Se desejar examinar espectros MS/MS, carregue também os arquivos MGF correspondentes.

### Etapa 2 — Definir as tolerâncias

Valores iniciais sugeridos:

| Parâmetro | Valor inicial | Interpretação |
|---|---:|---|
| **RT tolerance (min)** | `0.100` | Janela inicial de ±6 segundos |
| **MS1 tolerance (ppm)** | `5` | Erro máximo permitido de massa neutra |
| **Filter by observed ΔRT distribution** | Ativado | Aplica o filtro baseado no deslocamento observado |
| **Acceptance band (± k × MAD)** | `3.0` | Intervalo final centrado na mediana |

Esses valores são **pontos de partida**, não critérios universais. Se o deslocamento entre as aquisições for superior a 6 segundos, amplie a janela inicial de RT antes de executar.

### Etapa 3 — Executar o processamento

Clique em **Run Polarity Merge**. O aplicativo não recalcula automaticamente quando arquivos ou parâmetros são modificados.

Serão exibidas as contagens de nós POS, nós NEG, conexões originais e conexões candidatas entre polaridades.

### Etapa 4 — Revisar os candidatos

Na aba **Polarity matches**, examine:

- Identificadores dos nós POS e NEG.
- Valores de m/z dos precursores.
- Massas neutras reconstruídas.
- Erro de massa em ppm e Da.
- RT de cada polaridade.
- ΔRT com sinal e em valor absoluto.
- Identificador do par (`pair_id`).
- Campo editável **Keep edge** (`accepted`).

Desmarque pares que não devam ser considerados na análise. Essa decisão manual é independente da classificação automática de RT.

### Etapa 5 — Explorar o deslocamento cromatográfico

Na seção **Retention-time agreement**, observe:

- **Median ΔRT:** deslocamento central entre NEG e POS.
- **MAD ΔRT:** dispersão robusta ao redor da mediana.
- **Mean |ΔRT|:** média dos deslocamentos absolutos.
- **95th percentile |ΔRT|:** percentil 95 dos deslocamentos absolutos.
- Gráfico **RT POS × RT NEG**.
- Histograma de **ΔRT**, com mediana e linhas de ±2/±3 MAD.
- Classificação de cada par como **Inlier** ou **Outlier**.

Para `ΔRT = RT_NEG − RT_POS`:

- **ΔRT negativo:** o feature NEG eluiu antes do POS.
- **ΔRT positivo:** o feature NEG eluiu depois do POS.

O critério final é centrado no **deslocamento observado**, não obrigatoriamente em zero.

### Etapa 6 — Examinar a correlação de abundâncias (opcional)

Ative **Use cross-polarity intensity correlation** e carregue as duas tabelas quantitativas.

Na aba **Sample mapping & correlation**:

1. Escolha a coluna identificadora de features em cada tabela.
2. Revise a proposta automática de correspondência entre nomes de amostras.
3. Corrija manualmente as associações POS↔NEG.
4. Examine `intensity_pearson_r` e `n_paired_samples`.

A correspondência deve representar **as mesmas amostras biológicas**. Não associe injeções diferentes apenas porque os nomes são parecidos.

### Etapa 7 — Visualizar os espectros MS/MS (opcional)

Na aba **Spectral Pair Viewer**, selecione um `pair_id` e compare os espectros POS e NEG **lado a lado**.

Os controles incluem:

- **Peaks to display:** All, Top 10, Top 20, Top 50 ou Top 100.
- **Minimum relative intensity (%):** 0, 1, 2, 5 ou 10%.
- **Label most intense peaks:** 0, 3, 5, 10 ou 20.

A intensidade é normalizada separadamente para o pico base de cada espectro. O aplicativo **não utiliza representação espelhada nem calcula cosseno entre polaridades**.

### Etapa 8 — Exportar os resultados

Na aba **Export**, baixe:

| Arquivo | Conteúdo |
|---|---|
| `polarity_merge_network.graphml` | Redes POS e NEG originais com novas conexões cromatograficamente aceitas |
| `polarity_merge_edges_accepted.csv` | Pares que efetivamente serão incluídos como novas conexões |
| `polarity_merge_edges_rt_audit.csv` | Pares revisados e suas classificações de RT |

**Atenção:** a tabela de auditoria reflete o conjunto de candidatos mantidos na revisão manual. Para preservar também os candidatos desmarcados, exporte previamente **Download candidate edges (CSV)** na aba **Polarity matches**.

## Fundamentos do pareamento por massa

### Reconstrução da massa neutra

O aplicativo considera o par:

\[
[M+H]^+ \longleftrightarrow [M-H]^-
\]

Utilizando a massa do próton:

\[
m_p=1{,}007276466621\ \mathrm{Da}
\]

a massa neutra é estimada por:

\[
M_{POS}=(m/z)_{POS}-m_p
\]

\[
M_{NEG}=(m/z)_{NEG}+m_p
\]

A diferença entre as massas reconstruídas é:

\[
\Delta M=M_{NEG}-M_{POS}
\]

O erro assinado em ppm é:

\[
\mathrm{erro}_{ppm}=\frac{\Delta M}{|M_{POS}|}\times10^6
\]

A aceitação por massa exige que o **valor absoluto** do erro esteja dentro da tolerância selecionada em ppm ou Da.

Para esse par de íons, a diferença esperada entre os m/z dos precursores é:

\[
(m/z)_{POS}-(m/z)_{NEG}\approx2m_p
\]

\[
2m_p=2{,}014552933242\ \mathrm{Da}
\]

Entretanto, o aplicativo realiza a comparação por **massa neutra reconstruída**, e não apenas pela diferença direta entre os precursores.

### Ressalva sobre adutos

O algoritmo atual **não confirma automaticamente** que os íons observados sejam realmente `[M+H]+` e `[M−H]−`. Adutos de sódio, amônio, formiato ou cloreto, multímeros, fragmentos gerados na fonte, isotopólogos e atribuições incorretas podem produzir associações enganosas.

## Distribuição de ΔRT e filtragem por MAD

### Primeira etapa: geração de candidatos

A janela inicial exige:

\[
|RT_{NEG}-RT_{POS}|\leq RT_{\mathrm{tolerância}}
\]

Essa janela utiliza a diferença absoluta de RT e serve apenas para **selecionar candidatos**.

### Segunda etapa: estimativa do deslocamento observado

Para cada candidato:

\[
\Delta RT_i=RT_{NEG,i}-RT_{POS,i}
\]

O centro robusto da distribuição é:

\[
c=\mathrm{mediana}(\Delta RT_i)
\]

O desvio absoluto mediano é:

\[
MAD=\mathrm{mediana}(|\Delta RT_i-c|)
\]

Com multiplicador `k`, um candidato é classificado como **Inlier** quando:

\[
|\Delta RT_i-c|\leq k\times MAD
\]

ou, equivalentemente:

\[
c-k\times MAD\leq\Delta RT_i\leq c+k\times MAD
\]

A condição de fronteira é **inclusiva**. A versão v0.7 aplica uma pequena tolerância numérica para evitar classificações incorretas devido à representação de ponto flutuante.

### Por que mediana e MAD?

A mediana e o MAD são menos sensíveis a valores extremos do que média e desvio-padrão. São úteis quando a maioria dos pares representa uma tendência cromatográfica comum, mas há associações potencialmente incorretas.

**Detalhe metodológico:** o programa utiliza o **MAD bruto**, sem multiplicação pelo fator de consistência normal `1,4826`. Portanto, **±3 MAD não corresponde automaticamente a ±3 desvios-padrão** e não representa um intervalo de confiança ou um controle formal de taxa de falsos positivos.

### O que acontece com os outliers?

Os candidatos classificados como **Outlier**:

- **Não** são exportados como novas conexões `PolarityMerge` no GraphML.
- **Permanecem** na tabela de auditoria de RT.
- **Não provocam exclusão** dos nós originais.
- **Não removem** conexões intrapolaridade das redes originais.

### Considerações metodológicas

1. A janela inicial deve ser suficientemente ampla para capturar o deslocamento real. Um par excluído na primeira etapa não pode ser recuperado pelo filtro MAD.
2. A distribuição é calculada a partir de **conexões candidatas**, e não de identidades químicas previamente confirmadas.
3. Um mesmo nó pode participar de vários pares e, portanto, contribuir várias vezes para a distribuição.
4. A exclusão manual de candidatos pode modificar a mediana e o MAD utilizados na revisão e exportação.
5. O aplicativo não realiza sucessivas rodadas automáticas de exclusão de outliers.
6. Se houver menos de três valores finitos de ΔRT, a classificação informa **Insufficient data**. Com o filtro habilitado, esses pares não são exportados como conexões apoiadas por RT.
7. Se `MAD = 0`, são aceitos os valores efetivamente iguais à mediana; isso exige atenção em dados com RT arredondado ou discretizado.

**Recomendação científica:** sempre que possível, estime ou valide o deslocamento usando padrões, amostras QC ou pares POS↔NEG identificados independentemente. Inspecione também a estabilidade do deslocamento ao longo do gradiente cromatográfico.

## Correlação de abundâncias

Quando existem matrizes quantitativas e mapeamento correto das amostras, o programa pode calcular:

\[
r_{POS,NEG}=\mathrm{corr}(\mathbf{x}_{POS},\mathbf{x}_{NEG})
\]

onde cada vetor contém as intensidades do feature nas **mesmas amostras** analisadas em polaridades diferentes.

A correlação é **evidência complementar**, não um critério obrigatório para adicionar a conexão.

Uma correlação elevada não comprova identidade química; uma correlação baixa também não a exclui necessariamente. Efeitos de lote, supressão iônica, saturação, valores ausentes e diferenças de eficiência de ionização podem afetar `r`.

## Visualização dos espectros MS/MS

O visualizador foi desenvolvido para comparação **qualitativa** de espectros adquiridos em polaridades distintas.

- Os espectros POS e NEG são apresentados lado a lado.
- Cada espectro é normalizado independentemente para intensidade relativa.
- A seleção dos picos mais intensos é feita separadamente para cada polaridade.
- Os principais valores de m/z podem ser anotados.
- **Não é calculado cosseno espectral entre POS e NEG.**

As vias de fragmentação podem diferir substancialmente entre os modos de ionização. Os escores de similaridade existentes nas redes originais permanecem associados **apenas às respectivas conexões intrapolaridade**.

## Exportação e Cytoscape

### Organização dos identificadores

Para evitar colisões entre IDs iguais nas duas redes, o GraphML exportado utiliza prefixos:

```text
POS::1721
NEG::85
```

O ID original é preservado no atributo `original_node_id`, e a polaridade em `polarity`.

As conexões originais recebem `edge_type = MolecularNetworking` quando esse atributo não existia. As novas conexões recebem `edge_type = PolarityMerge`.

### Principais atributos das conexões

| Atributo | Significado |
|---|---|
| `pair_id` | Identificador do par candidato |
| `pos_node`, `neg_node` | IDs originais dos features |
| `pos_mz`, `neg_mz` | m/z dos precursores |
| `neutral_mass_pos`, `neutral_mass_neg` | Massas neutras reconstruídas |
| `delta_neutral_mass_da` | Diferença assinada de massa neutra |
| `mass_error_ppm` | Erro assinado em ppm |
| `delta_rt_sec` | Diferença assinada de RT |
| `abs_delta_rt_sec` | Diferença absoluta de RT |
| `rt_center_sec`, `rt_mad_sec` | Parâmetros da distribuição |
| `rt_deviation_sec`, `rt_robust_z` | Desvio em relação ao centro e razão desvio/MAD |
| `rt_status` | Classificação cromatográfica |
| `edge_type` | Tipo de conexão |

Quando disponíveis, também podem constar `intensity_pearson_r` e `n_paired_samples`.

### Importação no Cytoscape

1. Abra o Cytoscape.
2. Selecione **File → Import → Network from File**.
3. Carregue `polarity_merge_network.graphml`.
4. No painel **Style**, configure:
   - **Node Fill Color** por `polarity`.
   - **Node Label** por `original_node_id`.
   - **Edge Stroke Color** e **Line Type** por `edge_type`.
5. Diferencie visualmente `MolecularNetworking` e `PolarityMerge`.
6. Examine componentes, conexões entre polaridades e correspondências múltiplas.

Um estilo externo `Polarity_Merge_Cytoscape_Style.xml` pode ser importado separadamente. Esse arquivo **não é automaticamente incluído** na exportação da aplicação v0.7.

**Limitação de estrutura:** a exportação atual utiliza `networkx.Graph`, não `MultiGraph`. Redes que contenham arestas paralelas entre o mesmo par de nós podem ter essas arestas colapsadas. Verifique as contagens de conexões ao trabalhar com multigrafos.

## Correspondências múltiplas e ambiguidades

É possível que um nó POS corresponda a vários nós NEG, ou vice-versa. O Polarity Merge **preserva essas relações 1:N e N:1**.

Entre as explicações possíveis estão:

- **Peak splitting:** divisão artificial de um pico cromatográfico durante o processamento.
- Diferenças de integração e desconvolução entre polaridades.
- Isômeros próximos ou parcialmente coeluídos.
- Redundância de features ou atribuição incorreta de adutos.
- Concordância acidental de massa neutra e RT.

**Não se deve impor automaticamente um único par com base apenas no menor erro de massa.** É recomendável examinar cromatogramas de íons extraídos (EIC/XIC), formatos dos picos, espectros MS/MS e perfis quantitativos.

## Exemplo de interpretação

**Exemplo ilustrativo:** considere uma distribuição com:

- Mediana de ΔRT: **−4,5 s**
- MAD: **2,7 s**
- Multiplicador `k`: **3**

O intervalo de aceitação será:

\[
-4{,}5-3(2{,}7)\leq\Delta RT\leq-4{,}5+3(2{,}7)
\]

\[
\boxed{-12{,}6\ \mathrm{s}\leq\Delta RT\leq+3{,}6\ \mathrm{s}}
\]

Nesse cenário:

- **ΔRT = −5,0 s:** Inlier, próximo ao deslocamento observado.
- **ΔRT = +3,6 s:** Inlier, exatamente na fronteira inclusiva.
- **ΔRT = +7,2 s:** Outlier, fora do intervalo.

**Importante:** com a tolerância inicial padrão de **0,100 min (6 s)**, o par de +7,2 s sequer seria gerado como candidato. Para examiná-lo, é necessário aumentar a janela inicial (por exemplo, para 0,15 min) antes de executar.

## Limitações e controle de qualidade

1. **Associação não é identificação.** Concordância de massa e RT não confirma uma estrutura.
2. **Hipótese específica de adutos.** A versão atual modela `[M+H]+ ↔ [M−H]−`.
3. **Sem alinhamento cromatográfico automático.** O deslocamento é estimado para filtragem, mas os RT originais não são transformados.
4. **Possível circularidade.** O centro é estimado a partir dos próprios candidatos avaliados.
5. **Multiplicidade influencia a distribuição.** Um mesmo feature pode participar de vários pares.
6. **MAD não fornece p-valor ou FDR.** O limiar é exploratório.
7. **Correlação não é obrigatória.** Quando disponível, funciona como evidência complementar.
8. **Recuperação de MGF depende de IDs.** Identificadores inconsistentes impedem localizar espectros.
9. **Similaridade intrapolaridade não é similaridade entre polaridades.**
10. **Multigrafos podem exigir tratamento adicional** para preservar arestas paralelas.
11. **Mudanças de parâmetros exigem nova execução** do botão Run Polarity Merge.
12. **A auditoria depende da revisão manual.** Salve a tabela inicial de candidatos para documentar exclusões.

## Solução de problemas

| Problema | Verificação / solução |
|---|---|
| Nenhum candidato encontrado | Confira unidades de RT, atributos `mz`, hipótese de adutos e tolerâncias. |
| Distribuição de ΔRT truncada | Amplie a janela inicial de RT. |
| Muitos pares 1:N | Verifique peak splitting, coeluição, integração e redundância. |
| Espectro não encontrado | Confira a correspondência entre ID GraphML e `FEATURE_ID` no MGF. |
| Correlação indisponível | Forneça matrizes quantitativas completas e revise o mapeamento das amostras. |
| `Insufficient data` | Há menos de três candidatos com ΔRT finito para estimar a distribuição. |
| Mediana/MAD mudaram após revisão | A exclusão manual altera o conjunto usado no cálculo. |
| Estilo não aparece no Cytoscape | Confira os atributos `polarity` e `edge_type`; importe o XML de estilo separadamente. |
| Alterações nos controles não afetam os resultados | Clique novamente em **Run Polarity Merge**. |
| Falha na exportação GraphML | Verifique atributos incompatíveis, valores não escalares e tipo de grafo. |

## Reprodutibilidade e citação

Para documentar uma análise, registre:

- Versão do aplicativo ou commit Git utilizado.
- Arquivos de entrada e, preferencialmente, seus hashes.
- Hipótese de par iônico.
- Tolerância de massa em ppm ou Da.
- Janela inicial de RT.
- Ativação do filtro de distribuição, valor de `k`, mediana e MAD.
- Número de candidatos, inliers, outliers e exclusões manuais.
- Mapeamento de amostras e pré-processamento, quando houver correlação.
- Arquivos GraphML e CSV finais.

### Sugestão de texto para Materiais e Métodos

> As redes moleculares obtidas nos modos de ionização positivo e negativo foram integradas utilizando o Polarity Merge (LAABio, IPPN–UFRJ). As conexões candidatas entre polaridades foram geradas pela comparação de massas neutras reconstruídas sob a hipótese dos íons `[M+H]+` e `[M−H]−`, respeitando tolerâncias predefinidas de erro de massa e de tempo de retenção. O deslocamento cromatográfico assinado (`RT_NEG − RT_POS`) foi caracterizado pela mediana e pelo desvio absoluto mediano (MAD). Quando habilitado, o filtro de distribuição excluiu da exportação as conexões com ΔRT fora do intervalo mediana ± k × MAD. Todos os nós e as conexões intrapolaridade originais, incluindo seus atributos de similaridade espectral, foram preservados. Correspondências múltiplas não foram automaticamente reduzidas a relações um-para-um, e as associações foram interpretadas como putativas, não como identificações confirmadas.

### Como citar o software

Enquanto não houver DOI ou publicação específica, informe o autor, o repositório GitHub e a versão/commit utilizado:

**Ricardo M. Borges. Polarity Merge. LAABio, IPPN–UFRJ.**  
https://github.com/RicardoMBorges/polarity_merge_streamlit

### Licença

Este documento não estabelece uma licença de software. Adicione um arquivo `LICENSE` ao repositório caso deseje definir explicitamente as condições de uso e redistribuição.

---

**LAABio — Laboratório de Análise e Avaliação da Biodiversidade**  
**Instituto de Pesquisas de Produtos Naturais (IPPN)**  
**Universidade Federal do Rio de Janeiro (UFRJ)**  
**Polarity Merge | v0.7**
