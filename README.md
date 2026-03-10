# LA_HFpEF_RNAseq

<div align="center">

```text
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ LA_HFpEF_RNAseq                                                              ┃
┃ 左心房重构転写体档案 / レフトアトリウム・トランスクリプトーム解析                  ┃
┃ GHOST-SHELL // BIO-NET // HFpEF TWO-HIT MODEL                                ┃
┃ CTRL × 2HM :: DESeq2 · GSEA · StringTie · WGCNA · rMATS                      ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
```

**MULTI-LAYER TRANSCRIPTOMIC ATLAS OF LEFT ATRIAL REMODELLING**  
**IN A MURINE HFD + L-NAME TWO-HIT HFpEF MODEL**

[![STATUS: ONLINE](https://img.shields.io/badge/STATUS-ONLINE-00E5FF?style=for-the-badge)]()
[![PIPELINE: ACTIVE](https://img.shields.io/badge/PIPELINE-ACTIVE-FF2BD6?style=for-the-badge)]()
[![MODEL: CTRL_vs_2HM](https://img.shields.io/badge/MODEL-CTRL%20vs%202HM-FF8A00?style=for-the-badge)]()
[![SPECIES: MOUSE](https://img.shields.io/badge/SPECIES-MOUSE-A6FF00?style=for-the-badge)]()

[![Container: Apptainer SIF](https://img.shields.io/badge/CONTAINER-Apptainer%20SIF-1F6FEB?style=flat-square)](container/LA_HFpEF_baked.def)
[![Genome: GRCm39 / GENCODE M38](https://img.shields.io/badge/GENOME-GRCm39%20%2F%20GENCODE%20M38-00C2A8?style=flat-square)](https://www.gencodegenes.org/mouse/)
[![R ≥ 4.3](https://img.shields.io/badge/R-%E2%89%A5%204.3-7A5CFF?style=flat-square&logo=r)](https://www.r-project.org/)
[![Platform: Ubuntu 22.04](https://img.shields.io/badge/PLATFORM-Ubuntu%2022.04-E95420?style=flat-square&logo=ubuntu)](https://ubuntu.com/)
[![License: MIT](https://img.shields.io/badge/LICENSE-MIT-F7DF1E?style=flat-square)](LICENSE)
[![Output: Publication Ready](https://img.shields.io/badge/OUTPUT-Publication%20Ready-19C37D?style=flat-square)]()

</div>

---

<a id="toc"></a>
## 导航索引 // ナビゲーション // Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Pipeline Architecture](#pipeline-architecture)
- [Directory Structure](#directory-structure)
- [Quick Start](#quick-start)
  - [Native Execution](#native-execution)
  - [Container Execution](#container-execution)
- [Script Reference](#script-reference)
- [Input Requirements](#input-requirements)
- [Output Summary](#output-summary)
- [Computational Environment](#computational-environment)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)

---

<a id="overview"></a>
## 概览 // オーバービュー // Overview

> [!NOTE]
> **SYSTEM PROFILE**  
> This repository is a fully engineered, end-to-end bulk RNA-seq pipeline for characterising the **multi-layer transcriptomic landscape of the left atrium (LA)** in a murine **Heart Failure with preserved Ejection Fraction (HFpEF)** model.

### Experimental Frame // 実验配置

| Group | Label | n | Model |
|---|---|---:|---|
| Control | `CTRL` | 5 | Standard chow + water |
| Two-hit HFpEF | `2HM` | 5 | High-fat diet (HFD) + L-NAME |

### Analysis Layers // 解析层

| Layer | Mission | Core Outputs |
|---|---|---|
| Gene layer | Differential expression + systems response | DE tables, VST matrix, PCA, volcano, heatmap |
| Pathway layer | Hallmark / GO:BP / TF activity | GSEA curves, NES barplots, enrichment maps |
| Transcript layer | Novel transcript discovery | merged GTF, gffcompare, CPAT, transcript DE |
| Splicing layer | Differential alternative splicing | rMATS event tables, sashimi plots, enrichment |

This pipeline spans **gene-level DE**, **GSEA/pathway analysis**, **novel transcript discovery**, and **alternative splicing**, and is designed to generate publication-quality tables and figures suitable for high-impact cardiovascular genomics work.

All analyses are fully reproducible via a baked **Apptainer SIF** container. A single master script, `00_MASTER_RUN.sh`, orchestrates the full workflow with coloured logging, per-step timing, and automated I/O preflight checks.

---

<a id="scientific-background"></a>
## 科学背景 // バックグラウンド // Scientific Background

**HFpEF** now accounts for more than half of all heart failure cases and still lacks disease-modifying therapies. **Left atrial (LA) remodelling**—including fibrosis, electrical dysfunction, and metabolic reprogramming—is increasingly recognised as a major driver of symptoms and atrial fibrillation risk in HFpEF. The transcriptomic programs underlying LA remodelling, however, remain incompletely resolved.

This project uses the **two-hit mouse model**:

- **60 kcal% HFD**
- **L-NAME 0.5 g/L in drinking water**
- **15-week exposure**

The model recapitulates key human HFpEF phenotypes including obesity, hypertension, diastolic dysfunction, and exercise intolerance. LA tissue RNA-seq (**N = 5 per group**, polyA-enriched bulk sequencing) is analysed at three major biological resolution layers:

1. **Gene-level** — differential expression, pathway enrichment, TF activity, co-expression networks (WGCNA)
2. **Transcript-level** — novel transcript discovery (StringTie + gffcompare + CPAT), isoform-level differential expression
3. **Splicing-level** — differential alternative splicing events (rMATS-turbo), sashimi visualisation, functional enrichment

---

<a id="pipeline-architecture"></a>
## 线路架构 // パイプライン構造 // Pipeline Architecture

```text
FASTQ (raw, N=10)
    │
    ▼
[00_MASTER_RUN.sh] ────────────────────────────────────────────────────────────┐
    │                                                                          │
    ├── BRANCH 0 :: 共有上游 / Shared Upstream                                 │
    │   └── 01_upstream_pipeline.sh                                            │
    │        fastp           :: QC + trim                                      │
    │        STAR            :: splice-aware alignment (GRCm39)                │
    │        featureCounts   :: gene-level counts                              │
    │        MultiQC         :: aggregate QC report                            │
    │                                                                          │
    ├── BRANCH A :: 基因层 / Gene-Level Analysis                               │
    │   ├── 02_deseq2_gene.R           DESeq2 + apeglm shrinkage               │
    │   │                              VST matrix, volcano, PCA, heatmap, MA   │
    │   ├── 03_gsea_pathway.R          Hallmark + GO:BP via fgsea              │
    │   │                              NES barplots, GSEA curves               │
    │   ├── 04_gsea_advanced.R         enrichment map, TF inference            │
    │   │                              leading-edge heatmap, network plots      │
    │   └── 05_trait_wgcna.R           limma trait models + gene WGCNA         │
    │                                                                          │
    ├── BRANCH B :: 转录本发现 / Transcript Discovery                          │
    │   ├── 06_stringtie_pipeline.sh  assembly → merge → gffcompare           │
    │   │                              → re-quantification → prepDE            │
    │   ├── 07_novel_characterisation.sh  CPAT + cis-target annotation         │
    │   ├── 08_deseq2_novel_transcripts.R  isoform-level DESeq2               │
    │   └── 09_wgcna_coexpression.R   transcript-level WGCNA                   │
    │                                                                          │
    └── BRANCH C :: 可变剪接 / Alternative Splicing                            │
        ├── 10_rmats.sh               rMATS-turbo CTRL vs 2HM                 │
        │                             SE/MXE/A5SS/A3SS/RI filtering            │
        ├── 11_batch_sashimi.sh       rmats2sashimiplot batch                  │
        ├── 12_splicing_plots.R       event summary, PSI distributions         │
        └── 13_splicing_enrichment.R  GO/KEGG enrichment on spliced genes      │
                                                                             (logs)
```

---

<a id="directory-structure"></a>
## 目录结构 // ディレクトリ構成 // Directory Structure

```text
LA_HFpEF_RNAseq/
├── scripts/                         # All analysis scripts (numbered)
│   ├── 00_setup_env.sh
│   ├── 00_MASTER_RUN.sh             ← One-shot controller
│   ├── 01_upstream_pipeline.sh
│   ├── 02_deseq2_gene.R
│   ├── 03_gsea_pathway.R
│   ├── 04_gsea_advanced.R
│   ├── 05_trait_wgcna.R
│   ├── 06_stringtie_pipeline.sh
│   ├── 07_novel_characterisation.sh
│   ├── 08_deseq2_novel_transcripts.R
│   ├── 09_wgcna_coexpression.R
│   ├── 10_rmats.sh
│   ├── 11_batch_sashimi.sh
│   ├── 12_splicing_plots.R
│   └── 13_splicing_enrichment.R
│
├── container/
│   └── LA_HFpEF_baked.def           ← Apptainer definition file
│
├── envs/
│   ├── LA_HFpEF_all.yml             ← Primary conda env (r_rnaseq)
│   └── LA_HFpEF_rmats.yml           ← rMATS-dedicated env (rmats_env)
│
├── configs/
│   └── cardiac_sashimi_whitelist.txt ← Gene whitelist for sashimi
│
├── data/                            ← [NOT committed — user-provided]
│   ├── raw_fastq/                   # FASTQ files: {SampleID}_R1/R2.fastq.gz
│   └── trimmed_fastq/               # fastp output (auto-generated)
│
├── refs/                            ← [NOT committed — user-provided]
│   ├── genome/                      # GRCm39 primary assembly FASTA
│   ├── annotation/                  # GENCODE vM38 GTF
│   └── star_index/                  # Pre-built STAR genome index
│
├── results/                         ← All pipeline outputs
│   ├── 00_upstream/                 # fastp / STAR / featureCounts / MultiQC
│   ├── 01_gene_level/               # DESeq2 / GSEA / WGCNA outputs
│   ├── 02_transcript_level/         # StringTie / CPAT / novel DE outputs
│   └── 03_splicing/                 # rMATS / sashimi / enrichment outputs
│
├── logs/                            ← Timestamped run logs + version manifests
└── docs/
    ├── pipeline_io_map.tsv
    ├── FINAL_PIPELINE_AUDIT_REPORT.md
    └── SCRIPT_RENUMBERING.md
```

---

<a id="quick-start"></a>
## 快速启动 // クイックスタート // Quick Start

### Runtime Profile // 运行条件

| Item | Requirement |
|---|---|
| OS | **Ubuntu 22.04 LTS** (tested; other Linux distros should work) |
| RAM | **≥ 32 GB** (64 GB recommended for STAR indexing) |
| CPU | Modern multi-core processor (tuned for **AMD Ryzen 9 7900X**, 24 threads) |
| Storage | **≥ 200 GB** free |
| Execution mode | **Apptainer ≥ 1.1** or **Conda/Mamba** |

> [!IMPORTANT]
> **Reference data must be provided by the user.** See [Input Requirements](#input-requirements).

<a id="native-execution"></a>
### 本地执行 // ネイティブ実行 // Native Execution

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/LA_HFpEF_RNAseq.git
cd LA_HFpEF_RNAseq

# 2. Bootstrap the conda environment (r_rnaseq + rmats_env)
bash scripts/00_setup_env.sh

# 3. Place your FASTQs and reference files (see Input Requirements)

# 4. Run the full pipeline
bash scripts/00_MASTER_RUN.sh --cores 20

# --- OR selectively skip branches ---
bash scripts/00_MASTER_RUN.sh --cores 20 --skip-upstream
bash scripts/00_MASTER_RUN.sh --cores 20 --skip-splicing
bash scripts/00_MASTER_RUN.sh --cores 20 --skip-transcript --skip-splicing
```

**Available flags for `00_MASTER_RUN.sh`:**

| Flag | Description |
|---|---|
| `--cores N` | CPU threads for R and Bash tools (default: `20`) |
| `--project-root PATH` | Override project root (default: inferred from script path) |
| `--skip-setup` | Skip `00_setup_env.sh` (automatically set in containers) |
| `--skip-upstream` | Skip fastp / STAR / featureCounts / MultiQC |
| `--skip-gene` | Skip Branch A: DESeq2, GSEA, WGCNA (02–05) |
| `--skip-transcript` | Skip Branch B: StringTie, novel characterisation (06–09) |
| `--skip-splicing` | Skip Branch C: rMATS, sashimi, enrichment (10–13) |

<a id="container-execution"></a>
### 容器执行 // コンテナ実行 // Container Execution (Recommended)

Container execution guarantees exact software reproducibility across systems. All tools are baked into a single Apptainer SIF — **no internet access or runtime installs required**.

```bash
# 1. Build the SIF from the definition file (one-time; requires root or fakeroot)
apptainer build container/LA_HFpEF_baked.sif container/LA_HFpEF_baked.def

# 2. Run the full pipeline inside the container
REPO="$(pwd)"
apptainer exec \
  --cleanenv \
  --bind "${REPO}:/work" \
  container/LA_HFpEF_baked.sif \
  bash -lc 'cd /work && bash scripts/00_MASTER_RUN.sh --project-root /work --cores 20'

# 3. Test the container environment (runs %test block)
apptainer test container/LA_HFpEF_baked.sif
```

> [!TIP]
> **Container design principles**  
> - Both `r_rnaseq` and `rmats_env` conda environments are on `PATH` simultaneously  
> - Plotting is gracefully suppressed inside the container; computation and table export still run normally  
> - Script 11 (`11_batch_sashimi.sh`) exits cleanly in container mode  
> - Runtime package installs are strictly prohibited inside the SIF

---

<a id="script-reference"></a>
## 脚本索引 // スクリプト一覧 // Script Reference

| # | Script | Branch | Description | Key Tools |
|---|---|---|---|---|
| 00 | `00_setup_env.sh` | Setup | Bootstrap directory tree + conda environments | Conda/Mamba |
| 00 | `00_MASTER_RUN.sh` | Controller | One-shot orchestrator with coloured UI, timing, preflight | bash |
| 01 | `01_upstream_pipeline.sh` | Shared | fastp → STAR → featureCounts → MultiQC | fastp, STAR, subread, MultiQC |
| 02 | `02_deseq2_gene.R` | A | DESeq2 CTRL vs 2HM, apeglm shrinkage, VST, QC plots | DESeq2, apeglm |
| 03 | `03_gsea_pathway.R` | A | GSEA: MSigDB Hallmark + GO:BP via fgsea | fgsea, msigdbr |
| 04 | `04_gsea_advanced.R` | A | Enrichment map, TF activity (DoRothEA), leading-edge heatmap | decoupleR, igraph |
| 05 | `05_trait_wgcna.R` | A | limma trait models + gene-level WGCNA module-trait associations | limma, WGCNA |
| 06 | `06_stringtie_pipeline.sh` | B | StringTie assembly → merge → gffcompare → prepDE count matrices | StringTie, gffcompare |
| 07 | `07_novel_characterisation.sh` | B | CPAT coding potential + cis-target annotation for novel transcripts | CPAT, bedtools |
| 08 | `08_deseq2_novel_transcripts.R` | B | DESeq2 on novel (non-reference) transcripts | DESeq2 |
| 09 | `09_wgcna_coexpression.R` | B | Transcript-level WGCNA + guilt-by-association | WGCNA |
| 10 | `10_rmats.sh` | C | rMATS-turbo: 5 event types, CTRL vs 2HM, FDR/ΔPSI filtering | rMATS-turbo |
| 11 | `11_batch_sashimi.sh` | C | Batch sashimi plots for cardiac splicing hits | rmats2sashimiplot |
| 12 | `12_splicing_plots.R` | C | Event summary, PSI volcano, per-type barplots, results tables | ggplot2, patchwork |
| 13 | `13_splicing_enrichment.R` | C | GO/KEGG enrichment on differentially spliced gene sets | clusterProfiler |

---

<a id="input-requirements"></a>
## 输入要求 // 入力条件 // Input Requirements

### 1. Raw FASTQ Files

Place paired-end FASTQ files in `data/raw_fastq/`:

```text
data/raw_fastq/
├── CTRL_1_R1.fastq.gz    CTRL_1_R2.fastq.gz
├── CTRL_2_R1.fastq.gz    CTRL_2_R2.fastq.gz
├── CTRL_3_R1.fastq.gz    CTRL_3_R2.fastq.gz
├── CTRL_4_R1.fastq.gz    CTRL_4_R2.fastq.gz
├── CTRL_5_R1.fastq.gz    CTRL_5_R2.fastq.gz
├── 2HM_1_R1.fastq.gz     2HM_1_R2.fastq.gz
...
└── 2HM_5_R2.fastq.gz
```

### 2. Reference Genome & Annotation

```bash
# Download GENCODE M38 (GRCm39)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz

# Decompress and place in refs/
gunzip -c GRCm39.primary_assembly.genome.fa.gz > refs/genome/GRCm39.primary_assembly.genome.fa
gunzip -c gencode.vM38.primary_assembly.annotation.gtf.gz > refs/annotation/gencode.vM38.primary_assembly.annotation.gtf
```

### 3. STAR Index

Build the STAR genome index before running the upstream pipeline (requires ~30 GB RAM):

```bash
STAR --runMode genomeGenerate \
     --genomeDir refs/star_index \
     --genomeFastaFiles refs/genome/GRCm39.primary_assembly.genome.fa \
     --sjdbGTFfile refs/annotation/gencode.vM38.primary_assembly.annotation.gtf \
     --sjdbOverhang 149 \
     --runThreadN 20
```

### 4. Sample Metadata

Ensure `data/sample_metadata.csv` exists with columns:

```text
SampleID,Group,Condition,MouseID
CTRL_1,CTRL,Control,M01
...
2HM_5,2HM,HFpEF,M10
```

---

<a id="output-summary"></a>
## 输出总览 // 出力一覧 // Output Summary

| Branch | Key Deliverables | Location |
|---|---|---|
| Upstream | MultiQC HTML report, featureCounts gene matrix, alignment stats | `results/00_upstream/` |
| Gene DE | DESeq2 results table, VST matrix, volcano plot, PCA, heatmap | `results/01_gene_level/01_deseq2/` |
| GSEA | Hallmark + GO:BP enrichment tables, NES barplots, GSEA curves | `results/01_gene_level/02_gsea/` |
| Advanced GSEA | Enrichment map, TF activity scores, leading-edge heatmap | `results/01_gene_level/03_gsea_advanced/` |
| WGCNA (gene) | Module-trait heatmap, hub gene table, dendrogram | `results/01_gene_level/04_wgcna/` |
| StringTie | Merged GTF, gffcompare stats, transcript/gene count matrices | `results/02_transcript_level/01_stringtie/` |
| Novel transcripts | CPAT coding scores, novel transcript characterisation table | `results/02_transcript_level/02_cpat/` |
| Novel DE | Novel transcript DESeq2 results, volcano | `results/02_transcript_level/03_deseq2_novel/` |
| WGCNA (transcript) | Transcript co-expression modules, guilt-by-association | `results/02_transcript_level/04_wgcna/` |
| rMATS | Significant splicing events (SE/MXE/A5SS/A3SS/RI), ΔPSI tables | `results/03_splicing/01_rmats/` |
| Sashimi | Per-gene sashimi plots (PDF/PNG) for cardiac prioritised hits | `results/03_splicing/03_sashimi/` |
| Splicing plots | PSI volcano, event type barplots, summary tables | `results/03_splicing/04_plots/` |
| Splicing enrichment | GO/KEGG enrichment on spliced genes (dotplots, tables) | `results/03_splicing/05_enrichment/` |
| Logs | Per-step timestamped logs + software version manifest | `logs/` |

---

<a id="computational-environment"></a>
## 计算环境 // 計算環境 // Computational Environment

### Software Versions

| Tool | Version | Role |
|---|---|---|
| R | ≥ 4.3.0 | All R analyses |
| DESeq2 | ≥ 1.42 | Differential expression |
| apeglm | ≥ 1.24 | Log-fold-change shrinkage |
| WGCNA | ≥ 1.73 | Co-expression networks |
| fgsea | ≥ 1.28 | Fast GSEA |
| clusterProfiler | ≥ 4.10 | GO/KEGG enrichment |
| decoupleR | ≥ 2.8 | TF activity inference |
| fastp | ≥ 0.23 | FASTQ QC and trimming |
| STAR | ≥ 2.7.10 | Splice-aware alignment |
| featureCounts | ≥ 2.0 | Gene-level quantification |
| MultiQC | ≥ 1.21 | QC aggregation |
| StringTie | ≥ 2.2 | Transcript assembly |
| gffcompare | ≥ 0.12 | GTF annotation comparison |
| CPAT | ≥ 3.0 | Coding potential assessment |
| rMATS-turbo | ≥ 4.3 | Differential splicing |
| rmats2sashimiplot | ≥ 3.0 | Sashimi visualisation |
| samtools | ≥ 1.19 | BAM handling |

### Hardware Configuration

| Resource | Specification |
|---|---|
| CPU | AMD Ryzen 9 7900X (12 cores / 24 threads) |
| RAM | 64 GB DDR5 |
| OS | Ubuntu 22.04 LTS |
| Default threads | `N_CORES=20` |

### Conda Environments

Two isolated environments handle all dependencies:

- **`r_rnaseq`** — Primary environment: R ≥ 4.3, all Bioconductor/CRAN packages, fastp, STAR, featureCounts, StringTie, gffcompare, CPAT, samtools, MultiQC
- **`rmats_env`** — Dedicated rMATS environment: Python 3.8, rMATS-turbo, rmats2sashimiplot (isolated to prevent Python dependency conflicts)

```bash
# Recreate environments from lockfiles
conda env create -f envs/LA_HFpEF_all.yml      # r_rnaseq
conda env create -f envs/LA_HFpEF_rmats.yml    # rmats_env
```

---

<a id="reproducibility"></a>
## 可复现性 // 再現性 // Reproducibility

This pipeline is designed for full computational reproducibility:

1. **Baked Apptainer SIF** — all software pinned at build time; no network access required at runtime
2. **`set -euo pipefail`** — every Bash script exits immediately on error, unset variable, or pipe failure
3. **`set.seed(42)`** — all stochastic R operations (WGCNA `blockwiseModules`, permutation tests) use a fixed seed
4. **Timestamped logs** — every run generates a unique log in `logs/` with full stdout/stderr capture
5. **Software version manifest** — `logs/software_versions_<timestamp>.txt` records exact tool versions at each run
6. **Dplyr namespace prefixing** — all tidyverse calls use explicit namespaces (for example `dplyr::filter()`) to prevent masking conflicts
7. **Portable `PROJECT_ROOT` inference** — all scripts resolve paths relative to their own location; no hardcoded absolute paths

> [!WARNING]
> Reproducibility depends on correct placement of **FASTQ**, **reference genome**, **annotation**, and **sample metadata** files before execution.

---

<a id="citation"></a>
## 引用 // 引用情報 // Citation

If you use this pipeline or the results it produces, please cite:

> **[Manuscript in preparation]**  
> Haobo *et al.* — *Multi-layer transcriptomic landscape of left atrial remodelling in a murine HFpEF model.*  
> Loughrey Lab, University of Glasgow.

Please also cite the key underlying tools:

- Love *et al.* (2014) — DESeq2. *Genome Biology* 15:550
- Langfelder & Horvath (2008) — WGCNA. *BMC Bioinformatics* 9:559
- Subramanian *et al.* (2005) — GSEA. *PNAS* 102:15545
- Shen *et al.* (2014) — rMATS. *PNAS* 111:E5359
- Pertea *et al.* (2015) — StringTie. *Nature Biotechnology* 33:290
- Dobin *et al.* (2013) — STAR. *Bioinformatics* 29:15

---

<a id="acknowledgements"></a>
## 致谢 // 謝辞 // Acknowledgements

This work was performed in the **Loughrey Lab**, Institute of Cardiovascular & Medical Sciences, University of Glasgow. The two-hit HFpEF model and physiological phenotyping were conducted by the Loughrey Lab experimental team.

Computational resources were provided by the University of Glasgow computing infrastructure.

<div align="center">

```text
[ Loughrey Lab ] · [ University of Glasgow ] · [ Cardiovascular & Medical Sciences ]
SIGNAL STABLE :: ARCHIVE SEALED :: END OF FILE
```

</div>
