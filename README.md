# Pipeline Borrelia

> **Extraction et analyse phylogénétique des séquences de *Borrelia crocidurae* à partir de données de séquençage métagénomique**

---

## Structure du projet

```
projet_borrelia/
├── data/
│   ├── raw/                                        # Sequences netagenomiques
│   │      
│   │
│   ├── clean/                                       # Fastq apres trimming
│   │       ├── *_R1_clean.fastq.gz
│   │       └── *_R2_clean.fastq.gz
│   │
│   ├── kraken_db/                                   # Base de donnee kraken personaliser
│   │   
│   │   
│   │       
│   │
│   └── reference/
│       ├── Borrelia_crocidurae_Achema                 # Génomes de référence
|       |
│       ├── annotation/                                # Base de donnee d'annotation SnpEff
|       |        
|       ├── souches_refference/                         # Refference du Base de donnee kraken        
|       |
|       ├── phylogeny_genome/                           # Sequences complete de Borrelia
|       |
|       └── phylogeny_locus/                            # Sequences de Borrelia du locus detecter 
|
|
│
├── qc/                                                   # Controle qualites des sequences 
│   ├── fastqc/
│   └── multiqc/
│
├── kraken/                                                 # Classification, estimation et Borrelia crocidurae   |      ├── *.kraken                                           extraite 
│      ├── *.report
│      └── ...                                                  
│       
│
├── mapping/                                                  # Resultats du mapping
│   ├── bam/
│   │   └── *.bam
│   └── coverage/
│       └── *.txt
│
├── variants/                                                   # Resultats du variant calling et consensus
|
│
├── phylogeny/
|     ├── genome
|     |                                                # Resultats de la pylogenie par genome cpmplet
|     └── locus                                        # Resultats de la phylogenie par locus              
|
|
│
├── scripts/                                            # Scripts du pipeline
│   ├── install.sh                                     
│   ├── run_pipeline.sh                                  
│   ├── 01_qc.sh                                         
|   ├── 02_kraken.sh                                        
│   ├── 03_qc_borrelia.sh
|   ├── 04_mapping.sh
│   ├── 05_variants.sh
|   ├── 06_annotation.sh
│   ├── 07_genome_phylogeny.sh
│   ├── 08_locus_detection.sh
│   └── 09_locus_phylogeny.sh
│
├── results/
|     ├── logs/                                           # Journeaux d'execution
|     ├── annotation/                                     # Resultats de l'annotation
|     ├── locus_detection                                 # Resultats du locus detection
|     └── locus+phylogeny                                 # Resultas du locus phylogeny
|
├── environment.yml                                        # Outils du pipeline
├── README.md                                              # Presentation et explication du pipeline
└── .gitignore                                             # dossiers et fichiers exclus du depot
```

---

## Pipeline

```
Données métagénomiques (FASTQ — 2 lanes)
        │
        ▼
  Fusion des lanes L001 + L002
        │
        ▼
  [Étape 1] Contrôle qualité     → FastQC, Fastp, MultiQC
        │
        ▼
  [Étape 2] Extraction Borrelia  → Kraken2, Bracken, KrakenTools
        │
        ▼
  [Étape 3] Alignement           → BWA mem, Samtools
        │
        ▼
  [Étape 4] Variants & Consensus → bcftools
        │
        ▼
  [Étape 5] Phylogénie           → MAFFT, IQ-TREE2, FigTree/iTOL
        │
        ▼
  Identification du clade Borrelia
```

---

## Installation

### 1. Cloner le dépôt

```bash
git clone https://github.com/<username>/projet_borrelia.git
cd projet_borrelia
```

### 2. Installer les outils

```bash
bash scripts/install.sh
```

### 3. Activer l'environnement

```bash
conda activate borrelia_pipeline
```

---

## Configuration

Ouvrir `scripts/pipeline_borrelia.sh` et modifier :

```bash
THREADS=8                              # nombre de CPU disponibles
KRAKEN_DB="/path/to/kraken2_db"       # chemin vers la base Kraken2
```

---

## Utilisation

### Format des données attendu

```
data/raw/
├── sample1_L001_R1.fastq.gz
├── sample1_L001_R2.fastq.gz
├── sample1_L002_R1.fastq.gz
├── sample1_L002_R2.fastq.gz
├── sample2_L001_R1.fastq.gz
└── ...
```

### Un seul échantillon

```bash
conda activate borrelia_pipeline

bash scripts/pipeline_borrelia.sh \
    data/raw/sample_L001_R1.fastq.gz \
    data/raw/sample_L001_R2.fastq.gz \
    data/raw/sample_L002_R1.fastq.gz \
    data/raw/sample_L002_R2.fastq.gz
```

### Tous les échantillons d'un coup

```bash
conda activate borrelia_pipeline
bash scripts/run_all_samples.sh
```

---

## Résultats

| Dossier | Fichiers produits |
|---|---|
| `qc/` | `*_multiqc_report.html`, `*_fastp_report.html` |
| `kraken/` | `*_kraken_report.txt`, `*_borrelia_R1/R2.fastq` |
| `mapping/` | `*_aligned.bam`, `*_flagstat.txt`, `*_coverage.txt` |
| `variants/` | `*_variants_filtered.vcf.gz`, `*_consensus.fasta` |
| `phylogeny/` | `*_tree.treefile`, `*_aligned.fasta` |
| `results/` | `*_pipeline_log.txt` |

---

## Visualisation de l'arbre

```bash
# Ouvrir dans FigTree
figtree phylogeny/sample_tree.treefile
```

Ou en ligne : [https://itol.embl.de](https://itol.embl.de)

---

## Outils

| Outil | Version | Rôle |
|---|---|---|
| FastQC | 0.12.1 | Qualité des reads |
| Fastp | 0.23.4 | Trimming |
| MultiQC | 1.21 | Rapport global |
| Kraken2 | 2.1.3 | Classification taxonomique |
| Bracken | 3.1 | Abondance par espèce |
| KrakenTools | 1.2 | Extraction reads Borrelia |
| BWA | 0.7.18 | Alignement |
| Samtools | 1.21 | Manipulation BAM |
| bcftools | 1.21 | Variants et consensus |
| MAFFT | 7.525 | Alignement multiple |
| IQ-TREE2 | 2.3.6 | Phylogénie |

---

## Auteurs

- **Stagiaire** : Étudiante L3 Bioinformatique
- **Encadrant** : Ndongo — Master 2 Bioinformatique
