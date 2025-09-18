# PhytoIn 0.2.0

- Nova função `summary.param()` (S3) com saída concisa e robusta.
- `plot.param()` reescrito com seleção de tema via string/objeto **ggplot2**.
- `collector.curve()` e `rarefaction()` com intervalos de confiança (ribbons) e tema configurável.
- `AGB()` sem dependência de `rappdirs`; cache via `tools::R_user_dir()` e `BIOMASS::createCache()`.
- `BAplot()` corrigido (bindings globais) e compatível com parcelas retangulares e coordenadas individuais.
- Novos datasets: `quadrat2_plot.df`, `quadrat2_tree.df`, `quadrat3_rect.df`.
- Documentação, exemplos e NAMESPACE/DESCRIPTION revisados.

## First public release on CRAN

- `phytoparam()`, `summary.param()`, `plot.param()` para parâmetros fitossociológicos e índices de diversidade.
- `AGB()` para biomassa acima do solo, carbono e CO₂ equivalente (wrapper do **BIOMASS**).
- `stratvol()` para volume estratificado por classes de DAP.
- `collector.curve()` para curvas de acumulação de espécies.
- `rarefaction()` para rarefação individual com intervalos de confiança.
- `BAplot()` para visualização de áreas basais (parcelas retangulares e coordenadas individuais).
- Conjuntos de dados incluídos: `quadrat.df`, `point.df`, `quadrat2_plot.df`, `quadrat2_tree.df`, `quadrat3_rect.df`.

# PhytoIn 0.1.0 (development version)

- Versão interna inicial.
- Protótipos das funções centrais de análise fitossociológica.
- Não submetido ao CRAN.




