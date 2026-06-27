## TODO For Manuscript

### Missing Data/Planned Analyses

- [x] n for MGI-BioVU (add to scripts/manuscript/Figure 1A)
- [x] population details
- [ ] Use UKBB for vitamin D and update scripts
- [ ] Test HDL-C and triglyceride associations with serum cholesterol/BMD
- [ ] Pick methodology to account for UKBB overlapping
- [ ] **PCSK9 pQTL-based MR on heel BMD (UKB-PPP, ~1–2 days)** — Replace current 3-SNP cis PCSK9 instruments with protein-level instruments from the UKB Pharma Proteomics Project. Stronger null, higher power, cleaner interpretation. Makes the PCSK9-null leg of the HMGCR-specific / PCSK9-null argument bulletproof. Add result to drug-target MR BMD script (`drug_target_mr_bmd.qmd`).
- [ ] **Colocalization of HMGCR signals (LDL-C / BMD / HMGCR eQTL, ~1–2 days)** — Test whether the HMGCR locus signals for LDL-C, heel BMD, and HMGCR expression share a single causal variant (coloc H4) or are distinct (H3 = pleiotropy). Use `coloc` or `SuSiE-coloc` R packages. Shared variant strengthens causal MR interpretation and pre-empts reviewer concern about cis-instrument validity.
- [ ] **HMGCR eQTL-based MR across tissues (eQTLGen + GTEx, ~3–5 days)** — Use cis-eQTLs for HMGCR expression in liver/whole blood (eQTLGen) and bone-relevant peripheral tissues (GTEx: adipose, muscle, artery) as instruments for heel BMD and serum calcium. Tests whether the bone/calcium signal is driven by hepatic (systemic LDL-C lowering) or peripheral (local GGPP depletion in bone cells) HMGCR expression. A significant non-hepatic signal directly supports the cell-autonomous osteoclast hypothesis.