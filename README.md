# Mathematical modelling of SDH-b loss in chromaffin cells

This repository contains MATLAB code, parameter files, experimental data, and supporting model components used to investigate the metabolic consequences of succinate dehydrogenase subunit B (SDH-b) loss in chromaffin cells.

SDH links the tricarboxylic acid cycle to mitochondrial electron transport as Complex II. Loss of SDH-b disrupts this connection and produces the pseudohypoxic metabolic phenotype associated with SDH-deficient phaeochromocytoma. The model represented here integrates central carbon metabolism with mitochondrial energetics, electron transport, ion and volume regulation, and the proton motive force to examine how chromaffin cells adapt to this disruption.

## Scientific question

Why can SDH-b-deficient chromaffin cells retain mitochondrial Complex I function and metabolic fitness despite severe disruption of Complex II and the TCA cycle?

The modelling framework was developed alongside experimental measurements and 13C metabolic-flux analysis. It was used to investigate the consequences of SDH-b loss for metabolic fluxes, electron-transport-chain activity, mitochondrial membrane potential, ATP synthesis, proton balance, and mitochondrial volume.

## Main findings

The simulations indicate that retention of Complex I is associated with cofactor oxidation and helps SDH-b-deficient cells manage mitochondrial swelling and limit reversal of ATP synthase. The analysis also identifies mitochondrial proton leakage and control of the proton gradient across the inner mitochondrial membrane as important determinants of metabolic fitness following SDH-b loss.

## Associated publication

**Vera-Sigüenza E, Rana H, Nashebi R, Cloete I, Kl’učková K, Spill F, Tennant DA. _A Mathematical Exploration of SDH-b Loss in Chromaffin Cells._ Bulletin of Mathematical Biology 87, 53 (2025).**

[Read the published article](https://doi.org/10.1007/s11538-025-01427-z)

The article describes the biological rationale, experimental measurements, model construction, parameterisation, benchmarking, and SDH-b knockout simulations associated with this research code.

## Repository note

The main SDH modelling material is retained in its original research structure. The `archive/dynamical_systems_2022/` directory contains earlier Birmingham-era exploratory modelling material preserved for provenance and is not presented as part of the final published SDH-b model.

## Licence

The repository is distributed under the MIT License.