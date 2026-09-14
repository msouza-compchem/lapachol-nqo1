# Roteiro de execução (uso interno — não versionar no README público)

Este arquivo é para você, não para o comitê de admissão. O repositório público é
todo em inglês; este é o seu mapa de trabalho. Se preferir, deixe-o fora do git.

---

## Cronograma realista

| Semana | O que fazer | Tempo de máquina |
|--------|-------------|------------------|
| 1 | Instalar ORCA, Open Babel, GROMACS e o ambiente conda. Rodar `00_get_structure.sh`. Busca conformacional com CREST. | 2 h |
| 2 | Otimização + frequências (`01_opt_freq.inp`). Conferir ausência de frequência imaginária. | 2–3 h (deixe rodando à noite) |
| 3 | Single point com def2-TZVP e CPCM; ânion radicalar. Rodar `orca_descriptors.py`. Gerar os cubes e a figura no VMD. | 3–4 h |
| 4 | Preparar receptor 1D4A, validar o protocolo redockando duroquinona, docar o lapachol. | 30 min |
| 5–6 | Parametrizar o ligante no CGenFF, montar o sistema, minimizar e equilibrar. **Esta é a etapa que trava todo mundo.** | 1 dia de trabalho |
| 7 | Produção de 100 ns no Colab, com checkpoint no Drive. | 2–3 dias de relógio |
| 8 | `md_analysis.py`, preencher a seção 6 do README, revisar tudo em inglês. | 1 dia |

Dois meses trabalhando à noite e nos fins de semana. É factível.

---

## Os três pontos onde o projeto pode empacar

**1. Confôrmero errado.** O lapachol tem cadeia prenila flexível. Se você otimizar
a partir de um confôrmero ruim, os orbitais de fronteira saem sem sentido e você não
vai perceber. Rode o CREST antes. Custa 20 minutos e salva o projeto.

**2. Memória no ORCA.** `%maxcore` é por core. Com 8 GB, `nprocs 4` e
`%maxcore 1200` é o teto seguro. Se a máquina começar a usar swap, o cálculo não
trava com erro — ele só fica lento a ponto de parecer travado.

**3. Parametrização do FAD.** É o obstáculo real. Se consumir mais de uma semana,
tome a decisão pragmática: rode a MD sem o FAD e escreva isso na seção de limitações
do README. Reconhecer uma limitação por escrito é sinal de maturidade científica;
esconder é o contrário.

---

## Checklist antes de tornar o repositório público

- [ ] Nenhum arquivo acima de 50 MB (`git count-objects -vH` para conferir)
- [ ] README com as figuras aparecendo corretamente no GitHub
- [ ] Tudo em inglês, inclusive comentários de código e nomes de arquivo
- [ ] Seção 6 (Results) preenchida — repositório com resultado em branco é pior que
      repositório nenhum
- [ ] Uma figura forte no topo do README (a sobreposição HOMO/LUMO, ou o sítio ativo
      com o lapachol encaixado)
- [ ] `docs/REFERENCES.md` conferido contra os registros das editoras
- [ ] Descrição curta e tópicos no GitHub: `computational-chemistry`, `dft`,
      `molecular-dynamics`, `natural-products`, `drug-discovery`
- [ ] Link do repositório no currículo, na carta de intenção e na assinatura do e-mail

---

## Como transformar isso em conversa com orientador

Quando escrever para Gustavo Seabra (UF/CNPD3) ou Andrés Cisneros (UNT/CASCaM), o
projeto vira o segundo parágrafo do e-mail — o parágrafo que mostra que você já tem
pergunta própria. Algo assim:

> I recently put together a small reproducible pipeline combining DFT descriptors,
> docking and MD for lapachol, a naphthoquinone from a Brazilian tree, against NQO1
> [link]. It was run almost entirely on a laptop. What interests me about your group's
> work on [assunto específico deles] is [motivo específico], and I would like to ask
> whether you are considering students for the next cycle.

Curto, específico, com prova anexada. É isso que faz o professor responder.
