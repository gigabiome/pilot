# AMRBench wiki

**AMRBench** is an open benchmark for predicting antimicrobial resistance (MIC and R/S) from bacterial whole-genome sequences across the seven ESKAPE-E pathogens. It is built by Gigabiome (lead), UCL / Royal Free (clinical isolates and MIC testing) and causaLens (causal AI), and funded by Innovate UK, *Frontier AI Benchmarking Datasets*.

This wiki documents the **pilot**: a complete, working version of the benchmark that runs on synthetic data shaped like the real deliverables. When the clinical genomes arrive, the same code builds the real benchmark.

> **Status (October 2026, project month M2):** the pilot runs end to end on mock data, from mock genomes through splits, three baselines and harness scoring to the leaderboard. Real sequencing batch 1 is due at M3 (November 2026).

## Start here

| If you want to… | Read |
|---|---|
| understand what AMRBench is and how it fits the grant | [Project overview](Project-Overview.md) |
| run the pilot on your machine | [Getting started](Getting-Started.md) |
| understand how the pieces fit together | [Architecture](Architecture.md) |
| deliver data from the sequencing/QC pipeline | [Data contract](Data-Contract.md) |
| submit a model as an external team | [Submission format](Submission-Format.md) |
| see what the pilot has taught us so far | [Pilot findings](Pilot-Findings.md) |
| prepare for the real data | [Real-data migration](Real-Data-Migration.md) |

## All pages

**Concepts and design**
- [Project overview](Project-Overview.md): aims, consortium, tasks, what has changed since the proposal
- [Architecture](Architecture.md): pipeline stages, the two seams, directory layout
- [Labels and breakpoints](Labels-and-Breakpoints.md): MIC to log2 conversion, censoring, R/S calls
- [Splits and leakage](Splits-and-Leakage.md): clone clusters, leak-free splits, audits, manifests
- [Tasks and metrics](Tasks-and-Metrics.md): T1–T7 definitions, formulas, targets

**Components**
- [Data contract](Data-Contract.md): every input file, column by column
- [Mock data](Mock-Data.md): the simulated world and its planted ground truth
- [Features](Features.md): k-mers, pan-genome matrix, gene graphs
- [Baseline models](Baseline-Models.md): the model API, the three baselines, adding a model
- [Submission format](Submission-Format.md): the harness contract for any model
- [Configuration reference](Configuration-Reference.md): every key in `benchmark.yaml`

**Project**
- [Pilot findings](Pilot-Findings.md): results and lessons from the mock run
- [Real-data migration](Real-Data-Migration.md): checklist for switching to the clinical data
- [Decision log](Decision-Log.md): design decisions, made and open
- [Roadmap](Roadmap.md): pilot work mapped to work packages, deliverables and gates
- [Glossary](Glossary.md)
- [Troubleshooting & FAQ](Troubleshooting-FAQ.md)

## About this wiki
The pages live in `pilot/wiki/` next to the code, so documentation changes go through the same review as code changes. They follow GitHub-wiki conventions (`Home.md`, `_Sidebar.md`), so they can be pushed unchanged to the GitHub wiki once the repository is public.
