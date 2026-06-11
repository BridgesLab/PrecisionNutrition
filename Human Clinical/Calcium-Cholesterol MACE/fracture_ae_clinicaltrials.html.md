---
title: "Fracture Adverse Events in Lipid-Lowering Outcome Trials (ClinicalTrials.gov)"
author: "Dave Bridges"
date: today
format:
  html:
    toc: true
    toc-location: right
    keep-md: true
    code-fold: true
    code-summary: "Show the code"
  gfm:
    html-math-method: webtex
theme: journal
execute:
  echo: true
  warning: false
knitr:
  opts_chunk:
    fig-path: "figures-fracture-ae/"
    dev: ["png", "pdf"]
    dpi: 300
    dev.args:
      png:
        type: "cairo-png"
      pdf:
        family: "sans"
---


::: {.cell}

```{.r .cell-code}
library(httr2)
library(tidyverse)   # dplyr, purrr, stringr, tidyr, ggplot2
library(broom)
library(metafor)     # random-effects meta-analysis
library(knitr)
```
:::


## Motivation

This is a **preliminary, hypothesis-generating** scan of fracture-related adverse
events (AEs) reported on ClinicalTrials.gov for large lipid-lowering cardiovascular
outcome trials. The motivating question is whether lowering LDL-C — and by extension
the cholesterol/calcium axis we have been studying in the EHR cohort — is associated
with a signal in fracture risk. We pull the structured `adverseEventsModule` from the
v2 API, filter to fracture terms, and compare active vs. control arms within each
trial before pooling.

Two drug classes are included:

- **Non-statin LDL-lowering** (PCSK9 inhibitors, ezetimibe) — the original target set.
- **Statins** — added here for context, since they are the dominant LDL-lowering class
  and carry their own (mixed) literature on bone/fracture effects.

::: {.callout-warning}
## Caveats up front
- AE modules report **all-cause** fracture counts from safety surveillance, not
  adjudicated efficacy endpoints. Ascertainment differs across trials.
- Denominators (`numAtRisk`) are the safety-population arm sizes, not person-time, so
  these are **crude risk ratios**, not incidence-rate ratios. Differential follow-up
  between arms is *not* accounted for.
- Many landmark statin trials (4S, WOSCOPS, HPS, LIPID, CARDS) predate mandatory
  results posting and have **no machine-readable AE module** — they will silently
  return empty. Treat the statin set as a convenience sample of registry-era trials.
- A wrong NCT fails silently: it returns either an empty result *or — worse — a real
  but unrelated trial*. In an earlier draft `NCT00831441` (APPRAISE-2, apixaban) and
  `NCT00153231` (a vaginal-prolapse surgery study) had been mislabeled as statin
  trials and one slipped into the pool. Every NCT and **every arm title** below has now
  been confirmed against a live registry pull, and arm roles are assigned by an explicit
  hand-checked map rather than parsed from free-text arm titles.
:::

## Trial registry


::: {.cell}

```{.r .cell-code}
# class:      drug class for grouping
# comparator: "placebo" or "active" (more- vs less-intensive LDL lowering)
# has_ae:     TRUE = confirmed posted AE module on a live pull (2026-06); FALSE = no module
trials <- tribble(
  ~nct,          ~trial,             ~drug,          ~class,        ~comparator, ~has_ae,
  # --- Non-statin LDL-lowering (original set) ---
  "NCT01764633", "FOURIER",          "evolocumab",   "PCSK9i",      "placebo",   TRUE,
  "NCT01663402", "ODYSSEY OUTCOMES", "alirocumab",   "PCSK9i",      "placebo",   TRUE,
  "NCT00202878", "IMPROVE-IT",       "ezetimibe",    "ezetimibe",   "active",    TRUE,   # eze+simva vs simva
  "NCT00125593", "SHARP",            "eze+simva",    "ezetimibe",   "placebo",   TRUE,   # vs placebo (confounded)
  "NCT02993406", "CLEAR Outcomes",   "bempedoic",    "ATP-citrate", "placebo",   TRUE,
  # --- Statins ---
  "NCT00239681", "JUPITER",          "rosuvastatin", "statin",      "placebo",   TRUE,   # rosuva 20 vs placebo
  # Registry-era statin trials kept for the coverage table but with NO posted AE
  # module (confirmed empty on a live pull), so they contribute nothing downstream:
  "NCT00147602", "SPARCL",           "atorvastatin", "statin",      "placebo",   FALSE,  # atorva 80 vs placebo
  "NCT00327691", "SATURN",           "rosuvastatin", "statin",      "active",    FALSE,  # rosuva vs atorva
  "NCT00468923", "HOPE-3",           "rosuvastatin", "statin",      "placebo",   FALSE   # rosuva vs placebo (corrected NCT)
  # Dropped: NCT00831441 = APPRAISE-2 (apixaban), NCT00153231 = prolapse surgery.
  # Both were wrong NCTs masquerading as statin trials in an earlier draft.
)

trials |> arrange(class, trial) |> kable()
```

::: {.cell-output-display}


|nct         |trial            |drug         |class       |comparator |has_ae |
|:-----------|:----------------|:------------|:-----------|:----------|:------|
|NCT02993406 |CLEAR Outcomes   |bempedoic    |ATP-citrate |placebo    |TRUE   |
|NCT01764633 |FOURIER          |evolocumab   |PCSK9i      |placebo    |TRUE   |
|NCT01663402 |ODYSSEY OUTCOMES |alirocumab   |PCSK9i      |placebo    |TRUE   |
|NCT00202878 |IMPROVE-IT       |ezetimibe    |ezetimibe   |active     |TRUE   |
|NCT00125593 |SHARP            |eze+simva    |ezetimibe   |placebo    |TRUE   |
|NCT00468923 |HOPE-3           |rosuvastatin |statin      |placebo    |FALSE  |
|NCT00239681 |JUPITER          |rosuvastatin |statin      |placebo    |TRUE   |
|NCT00327691 |SATURN           |rosuvastatin |statin      |active     |FALSE  |
|NCT00147602 |SPARCL           |atorvastatin |statin      |placebo    |FALSE  |


:::
:::


Arm roles are assigned by an **explicit map** keyed on the registry's stable event-group
id (`EG000` / `EG001`), confirmed against a live pull of each trial's `eventGroups`. This
replaces the old free-text title heuristic, which silently mislabeled the lone-`Simvastatin`
control arm in IMPROVE-IT and an apixaban arm as "active". Only trials with a posted AE
module appear here.


::: {.cell}

```{.r .cell-code}
# role: "treatment" (LDL-lowering / more-intensive arm) vs "control"
arm_roles <- tribble(
  ~nct,          ~group_id, ~role,
  "NCT01764633", "EG000",   "control",     # FOURIER     Placebo
  "NCT01764633", "EG001",   "treatment",   # FOURIER     Evolocumab
  "NCT01663402", "EG000",   "control",     # ODYSSEY     Placebo
  "NCT01663402", "EG001",   "treatment",   # ODYSSEY     Alirocumab
  "NCT00202878", "EG000",   "treatment",   # IMPROVE-IT  Ezetimibe/Simvastatin
  "NCT00202878", "EG001",   "control",     # IMPROVE-IT  Simvastatin (active comparator)
  "NCT00125593", "EG000",   "treatment",   # SHARP       Simvastatin Plus Ezetimibe
  "NCT00125593", "EG001",   "control",     # SHARP       Placebo
  "NCT02993406", "EG000",   "treatment",   # CLEAR       Bempedoic Acid 180 mg
  "NCT02993406", "EG001",   "control",     # CLEAR       Placebo Comparator
  "NCT00239681", "EG000",   "control",     # JUPITER     Placebo
  "NCT00239681", "EG001",   "treatment"    # JUPITER     Rosuvastatin 20 mg
)
```
:::


## Pull adverse-event modules


::: {.cell}

```{.r .cell-code}
fetch_ae <- function(nct) {
  aem <- request(paste0("https://clinicaltrials.gov/api/v2/studies/", nct)) |>
    req_url_query(format = "json",
                  fields = "resultsSection.adverseEventsModule") |>
    req_user_agent("fracture-AE-pull (research; davebrid@umich.edu)") |>
    req_retry(max_tries = 3) |>
    req_perform() |>
    resp_body_json() |>
    pluck("resultsSection", "adverseEventsModule")

  if (is.null(aem)) return(tibble(nct = nct))

  groups <- map_dfr(aem$eventGroups,
                    \(g) tibble(group_id = g$id, arm = g$title))

  parse_events <- function(events, category) {
    map_dfr(events %||% list(), \(e) {
      map_dfr(e$stats, \(s) tibble(
        category     = category,
        term         = e$term,
        organ        = e$organSystem  %||% NA_character_,
        group_id     = s$groupId,
        num_affected = s$numAffected  %||% NA_integer_,
        num_at_risk  = s$numAtRisk    %||% NA_integer_
      ))
    })
  }

  bind_rows(parse_events(aem$seriousEvents, "serious"),
            parse_events(aem$otherEvents,   "other")) |>
    left_join(groups, by = "group_id") |>
    mutate(nct = nct, .before = 1)
}

# Be polite to the API; a small pause between calls.
ae_raw <- map_dfr(trials$nct, \(x) {
  out <- fetch_ae(x)
  Sys.sleep(0.5)
  out
})

# Did the live pull match what the registry table claims (has_ae)?
ae_coverage <- ae_raw |>
  group_by(nct) |>
  summarise(pulled_ae = any(!is.na(category)), .groups = "drop") |>
  right_join(trials, by = "nct") |>
  transmute(trial, drug, class, comparator,
            expected_ae = has_ae,
            pulled_ae = coalesce(pulled_ae, FALSE))

ae_coverage |> arrange(desc(pulled_ae), class, trial) |> kable()
```

::: {.cell-output-display}


|trial            |drug         |class       |comparator |expected_ae |pulled_ae |
|:----------------|:------------|:-----------|:----------|:-----------|:---------|
|CLEAR Outcomes   |bempedoic    |ATP-citrate |placebo    |TRUE        |TRUE      |
|FOURIER          |evolocumab   |PCSK9i      |placebo    |TRUE        |TRUE      |
|ODYSSEY OUTCOMES |alirocumab   |PCSK9i      |placebo    |TRUE        |TRUE      |
|IMPROVE-IT       |ezetimibe    |ezetimibe   |active     |TRUE        |TRUE      |
|SHARP            |eze+simva    |ezetimibe   |placebo    |TRUE        |TRUE      |
|JUPITER          |rosuvastatin |statin      |placebo    |TRUE        |TRUE      |
|HOPE-3           |rosuvastatin |statin      |placebo    |FALSE       |FALSE     |
|SATURN           |rosuvastatin |statin      |active     |FALSE       |FALSE     |
|SPARCL           |atorvastatin |statin      |placebo    |FALSE       |FALSE     |


:::
:::


## Filter to fracture terms


::: {.cell}

```{.r .cell-code}
# Broad fracture regex. We exclude "fracture" appearing only as part of unrelated
# MedDRA terms is unlikely, but inspect the distinct terms to be safe.
fracture <- ae_raw |>
  filter(str_detect(term, regex("fracture", ignore_case = TRUE))) |>
  left_join(trials, by = "nct")

# Sanity-check exactly which MedDRA terms were captured.
fracture |>
  distinct(term, organ) |>
  arrange(term) |>
  kable(caption = "Distinct fracture MedDRA terms captured")
```

::: {.cell-output-display}


Table: Distinct fracture MedDRA terms captured

|term                        |organ                                           |
|:---------------------------|:-----------------------------------------------|
|Acetabulum Fracture         |Injury, poisoning and procedural complications  |
|Acetabulum fracture         |Injury, poisoning and procedural complications  |
|Ankle Fracture              |Injury, poisoning and procedural complications  |
|Ankle fracture              |Injury, poisoning and procedural complications  |
|Avulsion fracture           |Injury, poisoning and procedural complications  |
|Cervical Vertebral Fracture |Injury, poisoning and procedural complications  |
|Cervical vertebral fracture |Injury, poisoning and procedural complications  |
|Chance Fracture             |Injury, poisoning and procedural complications  |
|Clavicle Fracture           |Injury, poisoning and procedural complications  |
|Clavicle fracture           |Injury, poisoning and procedural complications  |
|Comminuted fracture         |Injury, poisoning and procedural complications  |
|Compression Fracture        |Injury, poisoning and procedural complications  |
|Facial Bones Fracture       |Injury, poisoning and procedural complications  |
|Facial bones fracture       |Injury, poisoning and procedural complications  |
|Femoral Neck Fracture       |Injury, poisoning and procedural complications  |
|Femoral neck fracture       |Injury, poisoning and procedural complications  |
|Femur Fracture              |Injury, poisoning and procedural complications  |
|Femur fracture              |Injury, poisoning and procedural complications  |
|Fibula Fracture             |Injury, poisoning and procedural complications  |
|Fibula fracture             |Injury, poisoning and procedural complications  |
|Foot Fracture               |Injury, poisoning and procedural complications  |
|Foot fracture               |Injury, poisoning and procedural complications  |
|Forearm Fracture            |Injury, poisoning and procedural complications  |
|Forearm fracture            |Injury, poisoning and procedural complications  |
|Fracture                    |Injury, poisoning and procedural complications  |
|Fracture Displacement       |Injury, poisoning and procedural complications  |
|Fracture Malunion           |Musculoskeletal and connective tissue disorders |
|Fracture Nonunion           |Musculoskeletal and connective tissue disorders |
|Fracture displacement       |Injury, poisoning and procedural complications  |
|Fracture nonunion           |Musculoskeletal and connective tissue disorders |
|Fracture pain               |Musculoskeletal and connective tissue disorders |
|Fracture treatment          |Surgical and medical procedures                 |
|Fractured Ischium           |Injury, poisoning and procedural complications  |
|Fractured Sacrum            |Injury, poisoning and procedural complications  |
|Fractured Skull Depressed   |Injury, poisoning and procedural complications  |
|Fractured ischium           |Injury, poisoning and procedural complications  |
|Fractured sacrum            |Injury, poisoning and procedural complications  |
|Hand Fracture               |Injury, poisoning and procedural complications  |
|Hand fracture               |Injury, poisoning and procedural complications  |
|Hip Fracture                |Injury, poisoning and procedural complications  |
|Hip fracture                |Injury, poisoning and procedural complications  |
|Humerus Fracture            |Injury, poisoning and procedural complications  |
|Humerus fracture            |Injury, poisoning and procedural complications  |
|Ilium fracture              |Injury, poisoning and procedural complications  |
|Jaw Fracture                |Injury, poisoning and procedural complications  |
|Jaw fracture                |Injury, poisoning and procedural complications  |
|Limb fracture               |Injury, poisoning and procedural complications  |
|Lisfranc fracture           |Injury, poisoning and procedural complications  |
|Lower Limb Fracture         |Injury, poisoning and procedural complications  |
|Lower limb fracture         |Injury, poisoning and procedural complications  |
|Lumbar Vertebral Fracture   |Injury, poisoning and procedural complications  |
|Lumbar vertebral fracture   |Injury, poisoning and procedural complications  |
|Multiple Fractures          |Injury, poisoning and procedural complications  |
|Multiple fractures          |Injury, poisoning and procedural complications  |
|Open Fracture               |Injury, poisoning and procedural complications  |
|Open reduction of fracture  |Surgical and medical procedures                 |
|Osteoporotic Fracture       |Musculoskeletal and connective tissue disorders |
|Osteoporotic fracture       |Musculoskeletal and connective tissue disorders |
|Patella Fracture            |Injury, poisoning and procedural complications  |
|Patella fracture            |Injury, poisoning and procedural complications  |
|Pathological Fracture       |Musculoskeletal and connective tissue disorders |
|Pathological fracture       |Musculoskeletal and connective tissue disorders |
|Pelvic Fracture             |Injury, poisoning and procedural complications  |
|Pelvic fracture             |Injury, poisoning and procedural complications  |
|Periprosthetic Fracture     |Injury, poisoning and procedural complications  |
|Periprosthetic fracture     |Injury, poisoning and procedural complications  |
|Pubic rami fracture         |Injury, poisoning and procedural complications  |
|Pubis Fracture              |Injury, poisoning and procedural complications  |
|Pubis fracture              |Injury, poisoning and procedural complications  |
|Radius Fracture             |Injury, poisoning and procedural complications  |
|Radius fracture             |Injury, poisoning and procedural complications  |
|Rib Fracture                |Injury, poisoning and procedural complications  |
|Rib fracture                |Injury, poisoning and procedural complications  |
|Sacroiliac fracture         |Injury, poisoning and procedural complications  |
|Scapula Fracture            |Injury, poisoning and procedural complications  |
|Scapula fracture            |Injury, poisoning and procedural complications  |
|Skull Fracture              |Injury, poisoning and procedural complications  |
|Skull fracture              |Injury, poisoning and procedural complications  |
|Skull fractured base        |Injury, poisoning and procedural complications  |
|Spinal Compression Fracture |Injury, poisoning and procedural complications  |
|Spinal Fracture             |Injury, poisoning and procedural complications  |
|Spinal compression fracture |Injury, poisoning and procedural complications  |
|Spinal fracture             |Injury, poisoning and procedural complications  |
|Sternal Fracture            |Injury, poisoning and procedural complications  |
|Sternal fracture            |Injury, poisoning and procedural complications  |
|Thoracic Vertebral Fracture |Injury, poisoning and procedural complications  |
|Thoracic vertebral fracture |Injury, poisoning and procedural complications  |
|Tibia Fracture              |Injury, poisoning and procedural complications  |
|Tibia fracture              |Injury, poisoning and procedural complications  |
|Tooth Fracture              |Injury, poisoning and procedural complications  |
|Traumatic Fracture          |Injury, poisoning and procedural complications  |
|Traumatic fracture          |Injury, poisoning and procedural complications  |
|Ulna Fracture               |Injury, poisoning and procedural complications  |
|Ulna fracture               |Injury, poisoning and procedural complications  |
|Upper Limb Fracture         |Injury, poisoning and procedural complications  |
|Upper limb fracture         |Injury, poisoning and procedural complications  |
|Wrist Fracture              |Injury, poisoning and procedural complications  |
|Wrist fracture              |Injury, poisoning and procedural complications  |


:::
:::



::: {.cell}

```{.r .cell-code}
# Collapse to trial x arm x (serious/other). num_at_risk is constant within an
# event group, so max() recovers the denominator. Sum across the (possibly several)
# distinct fracture terms within an arm.
# Keep nct + group_id so we can join the explicit arm-role map (group_id is the
# registry's stable key; arm titles are free text and not trustworthy as a key).
fracture_by_arm <- fracture |>
  group_by(nct, trial, drug, class, comparator, group_id, arm, category) |>
  summarise(affected  = sum(num_affected, na.rm = TRUE),
            n_at_risk = max(num_at_risk,  na.rm = TRUE),
            .groups = "drop") |>
  mutate(pct = 100 * affected / n_at_risk)

fracture_by_arm |>
  arrange(class, trial, category, desc(arm)) |>
  select(trial, drug, class, arm, category, affected, n_at_risk, pct) |>
  kable(digits = 2, caption = "Fracture AE counts by trial, arm, and seriousness")
```

::: {.cell-output-display}


Table: Fracture AE counts by trial, arm, and seriousness

|trial            |drug         |class       |arm                   |category | affected| n_at_risk|  pct|
|:----------------|:------------|:-----------|:---------------------|:--------|--------:|---------:|----:|
|CLEAR Outcomes   |bempedoic    |ATP-citrate |Placebo Comparator    |serious  |       95|      6964| 1.36|
|CLEAR Outcomes   |bempedoic    |ATP-citrate |Bempedoic Acid 180 mg |serious  |       78|      7001| 1.11|
|FOURIER          |evolocumab   |PCSK9i      |Placebo               |serious  |      121|     13756| 0.88|
|FOURIER          |evolocumab   |PCSK9i      |Evolocumab            |serious  |      138|     13769| 1.00|
|ODYSSEY OUTCOMES |alirocumab   |PCSK9i      |Placebo               |serious  |      101|      9443| 1.07|
|ODYSSEY OUTCOMES |alirocumab   |PCSK9i      |Alirocumab            |serious  |       92|      9451| 0.97|
|IMPROVE-IT       |ezetimibe    |ezetimibe   |Simvastatin           |serious  |      198|      9077| 2.18|
|IMPROVE-IT       |ezetimibe    |ezetimibe   |Ezetimide/Simvastatin |serious  |      179|      9067| 1.97|
|JUPITER          |rosuvastatin |statin      |ROSUVASTATIN 20 MG    |serious  |       87|      8869| 0.98|
|JUPITER          |rosuvastatin |statin      |PLACEBO               |serious  |       92|      8864| 1.04|


:::
:::


## Assign arm roles (treatment vs. control)

Each event group is labeled by joining the explicit `arm_roles` map (keyed on the stable
`group_id`). The table below pairs each role with the **actual registry arm title** so you
can eyeball that "treatment" and "control" are not flipped — particularly for the
active-comparator trial (IMPROVE-IT) and the confounded SHARP comparison.


::: {.cell}

```{.r .cell-code}
fracture_roles <- fracture_by_arm |>
  inner_join(arm_roles, by = c("nct", "group_id"))

# Show the actual arm title behind each role so flips are easy to catch.
fracture_roles |>
  distinct(trial, class, comparator, group_id, arm, role) |>
  arrange(class, trial, role) |>
  kable(caption = "Arm-role assignment with registry arm titles — verify nothing is flipped")
```

::: {.cell-output-display}


Table: Arm-role assignment with registry arm titles — verify nothing is flipped

|trial            |class       |comparator |group_id |arm                   |role      |
|:----------------|:-----------|:----------|:--------|:---------------------|:---------|
|CLEAR Outcomes   |ATP-citrate |placebo    |EG001    |Placebo Comparator    |control   |
|CLEAR Outcomes   |ATP-citrate |placebo    |EG000    |Bempedoic Acid 180 mg |treatment |
|FOURIER          |PCSK9i      |placebo    |EG000    |Placebo               |control   |
|FOURIER          |PCSK9i      |placebo    |EG001    |Evolocumab            |treatment |
|ODYSSEY OUTCOMES |PCSK9i      |placebo    |EG000    |Placebo               |control   |
|ODYSSEY OUTCOMES |PCSK9i      |placebo    |EG001    |Alirocumab            |treatment |
|IMPROVE-IT       |ezetimibe   |active     |EG001    |Simvastatin           |control   |
|IMPROVE-IT       |ezetimibe   |active     |EG000    |Ezetimide/Simvastatin |treatment |
|JUPITER          |statin      |placebo    |EG000    |PLACEBO               |control   |
|JUPITER          |statin      |placebo    |EG001    |ROSUVASTATIN 20 MG    |treatment |


:::
:::


::: {.callout-note}
**IMPROVE-IT** is an *active-comparator* trial (ezetimibe+simvastatin vs simvastatin
alone), so its "treatment vs control" contrast is *more- vs less-intensive* LDL lowering,
not drug-vs-placebo. **SHARP** randomized eze+simva vs placebo, so its effect cannot be
attributed to ezetimibe alone. Both are flagged `comparator` in the registry and IMPROVE-IT
is excluded from the placebo-controlled pool below.
:::

## Per-trial risk ratios


::: {.cell}

```{.r .cell-code}
# 2x2 per trial x seriousness: treatment vs control. Keep the arm titles so the
# RR table is self-documenting. Drop rows lacking a clean two-arm contrast.
rr_input <- fracture_roles |>
  group_by(trial, drug, class, comparator, category, role) |>
  summarise(affected = sum(affected), n_at_risk = max(n_at_risk),
            arm = first(arm), .groups = "drop") |>
  pivot_wider(names_from = role,
              values_from = c(affected, n_at_risk, arm)) |>
  filter(!is.na(affected_treatment), !is.na(affected_control),
         n_at_risk_treatment > 0, n_at_risk_control > 0)

# Risk ratio with Wald CI on the log scale; Haldane–Anscombe 0.5 continuity
# correction when any cell is zero.
rr_tab <- rr_input |>
  mutate(
    a = affected_treatment, n1 = n_at_risk_treatment,
    c = affected_control,   n0 = n_at_risk_control,
    cc = if_else(a == 0 | c == 0, 0.5, 0),
    rr      = ((a + cc) / (n1 + 2 * cc)) / ((c + cc) / (n0 + 2 * cc)),
    log_rr  = log(rr),
    se_log  = sqrt(1 / (a + cc) - 1 / (n1 + 2 * cc) +
                   1 / (c + cc) - 1 / (n0 + 2 * cc)),
    lci = exp(log_rr - 1.96 * se_log),
    uci = exp(log_rr + 1.96 * se_log)
  )

rr_tab |>
  transmute(trial, drug, class, comparator, category,
            `treatment arm` = str_glue("{arm_treatment} ({a}/{n1})"),
            `control arm`   = str_glue("{arm_control} ({c}/{n0})"),
            rr = round(rr, 2),
            `95% CI` = str_glue("{round(lci,2)}–{round(uci,2)}")) |>
  arrange(class, category, trial) |>
  kable(caption = "Crude fracture risk ratio (treatment vs control) by trial, with arm labels")
```

::: {.cell-output-display}


Table: Crude fracture risk ratio (treatment vs control) by trial, with arm labels

|trial            |drug         |class       |comparator |category |treatment arm                    |control arm                  |   rr|95% CI    |
|:----------------|:------------|:-----------|:----------|:--------|:--------------------------------|:----------------------------|----:|:---------|
|CLEAR Outcomes   |bempedoic    |ATP-citrate |placebo    |serious  |Bempedoic Acid 180 mg (78/7001)  |Placebo Comparator (95/6964) | 0.82|0.61–1.1  |
|FOURIER          |evolocumab   |PCSK9i      |placebo    |serious  |Evolocumab (138/13769)           |Placebo (121/13756)          | 1.14|0.89–1.45 |
|ODYSSEY OUTCOMES |alirocumab   |PCSK9i      |placebo    |serious  |Alirocumab (92/9451)             |Placebo (101/9443)           | 0.91|0.69–1.21 |
|IMPROVE-IT       |ezetimibe    |ezetimibe   |active     |serious  |Ezetimide/Simvastatin (179/9067) |Simvastatin (198/9077)       | 0.91|0.74–1.11 |
|JUPITER          |rosuvastatin |statin      |placebo    |serious  |ROSUVASTATIN 20 MG (87/8869)     |PLACEBO (92/8864)            | 0.95|0.71–1.27 |


:::
:::


## Forest plot and random-effects pooling

We pool **serious** fracture events across placebo-controlled trials with a
random-effects model (REML). Active-comparator trials are excluded from the pooled
estimate but shown in the per-trial table above.


::: {.cell}

```{.r .cell-code}
# Restrict to serious events and to placebo-controlled trials (per the registry
# `comparator` flag — no hardcoded trial list). Active-comparator trials (IMPROVE-IT)
# are excluded automatically. Empty-AE statin trials never reach here.
meta_dat <- rr_tab |>
  filter(category == "serious", comparator == "placebo")

message("Pooling ", nrow(meta_dat), " placebo-controlled trials: ",
        paste(meta_dat$trial, collapse = ", "))

if (nrow(meta_dat) >= 2) {
  es <- escalc(measure = "RR",
               ai = a, n1i = n1,    # active arm: events, n
               ci = c, n2i = n0,    # control arm: events, n
               data = meta_dat,
               slab = trial)

  res <- rma(yi, vi, data = es, method = "REML")

  forest(res, atransf = exp,
         at = log(c(0.25, 0.5, 1, 2, 4)),
         xlab = "Risk ratio (serious fracture, active vs control)",
         mlab = "RE model (REML)",
         header = c("Trial", "RR [95% CI]"))

  res
} else {
  message("Not enough trials with serious fracture data to pool — ",
          "verify the unverified NCTs and re-run.")
}
```

::: {.cell-output-display}
![](figures-fracture-ae/meta-1.png){width=2400}
:::

::: {.cell-output .cell-output-stdout}

```

Random-Effects Model (k = 4; tau^2 estimator: REML)

tau^2 (estimated amount of total heterogeneity): 0.0030 (SE = 0.0187)
tau (square root of estimated tau^2 value):      0.0551
I^2 (total heterogeneity / total variability):   13.21%
H^2 (total variability / sampling variability):  1.15

Test for Heterogeneity:
Q(df = 3) = 3.1897, p-val = 0.3633

Model Results:

estimate      se     zval    pval    ci.lb   ci.ub    
 -0.0395  0.0756  -0.5219  0.6018  -0.1877  0.1088    

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::
:::



::: {.cell}

```{.r .cell-code}
# Subgroup pooling by drug class, if enough trials returned data.
if (exists("res") && nrow(meta_dat) >= 2) {
  es <- escalc(measure = "RR", ai = a, n1i = n1, ci = c, n2i = n0,
               data = meta_dat, slab = trial)
  res_class <- rma(yi, vi, mods = ~ class, data = es, method = "REML")
  res_class
}
```

::: {.cell-output .cell-output-stdout}

```

Mixed-Effects Model (k = 4; tau^2 estimator: REML)

tau^2 (estimated amount of residual heterogeneity):     0.0073 (SE = 0.0357)
tau (square root of estimated tau^2 value):             0.0853
I^2 (residual heterogeneity / unaccounted variability): 28.85%
H^2 (unaccounted variability / sampling variability):   1.41
R^2 (amount of heterogeneity accounted for):            0.00%

Test for Residual Heterogeneity:
QE(df = 1) = 1.4054, p-val = 0.2358

Test of Moderators (coefficients 2:3):
QM(df = 2) = 1.2677, p-val = 0.5306

Model Results:

             estimate      se     zval    pval    ci.lb   ci.ub    
intrcpt       -0.2025  0.1742  -1.1623  0.2451  -0.5439  0.1389    
classPCSK9i    0.2322  0.2070   1.1219  0.2619  -0.1734  0.6378    
classstatin    0.1460  0.2445   0.5973  0.5503  -0.3331  0.6252    

---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::
:::


## Session info


::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

```
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.5

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Detroit
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] knitr_1.51          metafor_5.0-1       numDeriv_2016.8-1.1
 [4] metadat_1.6-0       Matrix_1.7-5        broom_1.0.13       
 [7] lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0      
[10] dplyr_1.2.1         purrr_1.2.2         readr_2.2.0        
[13] tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.3      
[16] tidyverse_2.0.0     httr2_1.2.2        

loaded via a namespace (and not attached):
 [1] rappdirs_0.3.4     generics_0.1.4     stringi_1.8.7      lattice_0.22-9    
 [5] hms_1.1.4          digest_0.6.39      magrittr_2.0.5     evaluate_1.0.5    
 [9] grid_4.6.0         timechange_0.4.0   RColorBrewer_1.1-3 fastmap_1.2.0     
[13] jsonlite_2.0.0     backports_1.5.1    scales_1.4.0       codetools_0.2-20  
[17] cli_3.6.6          rlang_1.2.0        withr_3.0.2        yaml_2.3.12       
[21] otel_0.2.0         tools_4.6.0        tzdb_0.5.0         mathjaxr_2.0-0    
[25] curl_7.1.0         vctrs_0.7.3        R6_2.6.1           lifecycle_1.0.5   
[29] htmlwidgets_1.6.4  pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
[33] glue_1.8.1         xfun_0.57          tidyselect_1.2.1   rstudioapi_0.18.0 
[37] farver_2.1.2       htmltools_0.5.9    nlme_3.1-169       rmarkdown_2.31    
[41] compiler_4.6.0     S7_0.2.2          
```


:::
:::


## Next steps

- **Verify the `FALSE` NCTs** in the registry table — several statin entries are
  best-guesses and may lack a posted AE module (they'll show `has_ae = FALSE` above).
- Swap crude risk ratios for **incidence-rate ratios** where person-time is recoverable
  (some modules expose `frequencyThreshold` / time-frame metadata).
- Hand-correct arm roles for active-comparator trials and decide whether to keep them.
- If a real signal emerges, cross-reference against the EHR-derived osteoporosis /
  fracture competing-risks analyses in this project for triangulation.
```

