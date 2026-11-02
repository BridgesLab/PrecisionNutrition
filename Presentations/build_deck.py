"""Build the Surface Lab Update slide deck — clean version."""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from PIL import Image
import os

# ── Colors ──────────────────────────────────────────────────────────────────
BLUE      = RGBColor(0x00, 0x27, 0x4C)
MAIZE     = RGBColor(0xFF, 0xCB, 0x05)
WHITE     = RGBColor(0xFF, 0xFF, 0xFF)
LGRAY     = RGBColor(0xF5, 0xF5, 0xF5)
DGRAY     = RGBColor(0x40, 0x40, 0x40)
MBLU      = RGBColor(0x1A, 0x5C, 0x96)

SLIDE_W   = Inches(13.333)
SLIDE_H   = Inches(7.5)

FIG_BASE  = "/Users/davebrid/Documents/GitHub/PrecisionNutrition"
SCRATCH   = "/private/tmp/claude-503/-Users-davebrid-Documents-GitHub-PrecisionNutrition/09589a09-42a9-48c8-ad83-9849e1ca0e49/scratchpad"
GRANT     = "/Users/davebrid/Documents/GitHub/GrantApplications/NIH/R01/2026-07 Calcium Cholesterol/figures"

FIGURES = {
    "mouse_scatter":    f"{FIG_BASE}/Mouse Genetics/Cholesterol/figures/calcium-stratification-1.png",
    "do_tree":          f"{FIG_BASE}/Mouse Genetics/Cholesterol/figures/ggparty-bootstrapped-tree-1.png",
    "meta_forest":      f"{FIG_BASE}/Meta Analysis/figures/calcium-bayesian-forest-1.png",
    "mr_summary":       f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/mr-results-summary-2.png",
    "bayesian_dir":     f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/cholesterol-calcium-bayesian-inference-2.png",
    "bmd_forest":       f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/bone-forest-primary-1.png",
    "drug_target_both": f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/drug-target-mr-bmd-ca-1.png",
    "bayes_bar":        f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/bayes-bar-plot-1.png",
    "mechanisms":       f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/summary-mechanisms-1.png",
    "mech_primary":     f"{FIG_BASE}/Human Genetics/calcium-cholesterol/figures/summary-primary-mechanisms-1.png",
    "statin_ldlc":      f"{FIG_BASE}/Human Clinical/Calcium-Cholesterol MACE/figures-osteo-cr/statin_primary_model_only.png",
    "biochem":          f"{SCRATCH}/biochemistry-framework.png",
}


# ── Helpers ──────────────────────────────────────────────────────────────────
prs = Presentation()
prs.slide_width  = SLIDE_W
prs.slide_height = SLIDE_H
blank = prs.slide_layouts[6]


def rect(sl, l, t, w, h, fill=None, line=None):
    s = sl.shapes.add_shape(1, l, t, w, h)
    if fill:
        s.fill.solid(); s.fill.fore_color.rgb = fill
    else:
        s.fill.background()
    if line:
        s.line.color.rgb = line
    else:
        s.line.fill.background()
    return s


def txt(sl, text, l, t, w, h, sz=18, bold=False, italic=False,
        color=None, align=PP_ALIGN.LEFT, wrap=True):
    tb = sl.shapes.add_textbox(l, t, w, h)
    tf = tb.text_frame; tf.word_wrap = wrap
    p = tf.paragraphs[0]; p.alignment = align
    r = p.add_run(); r.text = text
    r.font.size = Pt(sz); r.font.bold = bold; r.font.italic = italic
    if color: r.font.color.rgb = color
    return tb


def bullets(sl, items, l, t, w, h, sz=14, color=None):
    tb = sl.shapes.add_textbox(l, t, w, h)
    tf = tb.text_frame; tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        indent = 1 if item.startswith("  ") else 0
        p.level = indent
        r = p.add_run(); r.text = item.strip()
        r.font.size = Pt(sz - indent * 2)
        if color: r.font.color.rgb = color
    return tb


def fig(sl, path, l, t, max_w, max_h):
    """Add image preserving aspect ratio within (max_w × max_h) box, centred."""
    if not os.path.exists(path):
        return False
    with Image.open(path) as im:
        iw, ih = im.size
    scale = min(max_w / iw, max_h / ih)
    w = iw * scale; h = ih * scale
    offset_l = l + (max_w - w) / 2
    offset_t = t + (max_h - h) / 2
    sl.shapes.add_picture(path, offset_l, offset_t, w, h)
    return True


def header(sl, title, sub=None):
    rect(sl, 0, 0, SLIDE_W, Inches(1.1), fill=BLUE)
    txt(sl, title, Inches(0.3), Inches(0.08), Inches(12.7), Inches(0.65),
        sz=28, bold=True, color=WHITE)
    if sub:
        txt(sl, sub, Inches(0.3), Inches(0.72), Inches(12.7), Inches(0.34),
            sz=15, color=MAIZE)
    rect(sl, 0, Inches(1.1), SLIDE_W, Inches(0.04), fill=MAIZE)


def footer(sl, label="Bridges Lab · University of Michigan"):
    rect(sl, 0, Inches(7.2), SLIDE_W, Inches(0.3), fill=BLUE)
    txt(sl, label, Inches(0.2), Inches(7.21), Inches(12), Inches(0.25),
        sz=9, color=RGBColor(0xCC, 0xCC, 0xCC))


# ─── SLIDE 1: TITLE ──────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=BLUE)
rect(sl, 0, Inches(5.9), SLIDE_W, Inches(0.12), fill=MAIZE)

txt(sl, "Cholesterol and Bone Demineralization:",
    Inches(0.8), Inches(1.4), Inches(11.7), Inches(1.1),
    sz=36, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
txt(sl, "From Observational Associations to Causal Mechanisms",
    Inches(0.8), Inches(2.5), Inches(11.7), Inches(0.8),
    sz=26, color=MAIZE, align=PP_ALIGN.CENTER)
txt(sl, "Dave Bridges Lab  |  University of Michigan",
    Inches(0.8), Inches(3.5), Inches(11.7), Inches(0.5),
    sz=18, color=WHITE, align=PP_ALIGN.CENTER)
txt(sl, "Surface Lab Update  —  November 13, 2026",
    Inches(0.8), Inches(4.1), Inches(11.7), Inches(0.4),
    sz=16, color=RGBColor(0xAA, 0xCC, 0xFF), align=PP_ALIGN.CENTER)
txt(sl, "Seven converging lines of evidence: DO mice · meta-analysis · bidirectional MR · "
       "HMGCR drug-target MR · BMD mechanism MR · robust overlap-corrected MR · EHR statin analysis",
    Inches(1.5), Inches(5.1), Inches(10.3), Inches(1.0),
    sz=12, italic=True, color=RGBColor(0xCC, 0xCC, 0xFF), align=PP_ALIGN.CENTER)


# ─── SLIDE 2: DO MICE ────────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "The Discovery: Calcium–Cholesterol in Diversity Outbred Mice",
       "Cousineau et al., Nutrients 2024  |  Svenson cohort, n = 850  |  NCD + HFHS diets")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(4.4), Inches(5.85), fill=WHITE)
txt(sl, "Serendipitous Discovery", Inches(0.25), Inches(1.3), Inches(4.2), Inches(0.38),
    sz=15, bold=True, color=BLUE)
bullets(sl, [
    "Hypothesis-naïve bootstrapped regression tree",
    "Serum calcium: top predictor of cholesterol",
    "  (after diet and triglycerides)",
    "1 mg/dL ↑ Ca → 12.7 ± 0.8 mg/dL ↑ Chol",
    "  p = 3.0 × 10⁻⁴³;  ω²ₚ = 22%",
    "Spearman ρ ≈ 0.40–0.50 across conditions",
    "Replicated in both sexes and both diets",
    "",
    "→ Top panel: tree (Ca = first split)",
    "→ Bottom panel: Ca–Chol scatter by strata",
],
    Inches(0.25), Inches(1.72), Inches(4.2), Inches(5.1), sz=13, color=DGRAY)

# Regression tree (top right)
fig(sl, FIGURES["do_tree"],      Inches(4.7), Inches(1.2),  Inches(8.45), Inches(2.9))
# Calcium-cholesterol scatter (bottom right)
fig(sl, FIGURES["mouse_scatter"],Inches(4.7), Inches(4.15), Inches(8.45), Inches(2.9))


# ─── SLIDE 3: META-ANALYSIS ──────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Meta-Analysis of Calcium–Cholesterol Associations in Humans",
       "10 studies · 88,939 participants · 17 subgroups")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(4.0), Inches(5.85), fill=WHITE)
txt(sl, "Summary Statistics", Inches(0.25), Inches(1.3), Inches(3.8), Inches(0.4),
    sz=15, bold=True, color=BLUE)
bullets(sl, [
    "Frequentist meta-analysis:",
    "  r = 0.17 ± 0.018,  p = 4.2 × 10⁻⁸",
    "",
    "Bayesian hierarchical model:",
    "  Prior from DO mouse data (ρ ≈ 0.45)",
    "  Posterior r = 0.176 (95% CI: 0.140–0.214)",
    "  P(r > 0) = 1.000 (Bayes Factor > 100)",
    "  P(r < 0.39) = 1.000 — attenuated vs mice",
    "",
    "Between-study variance:  τ = 0.071",
    "",
    "Global reach: Europe, Asia, Middle East",
    "Both sexes; ages 20–75+",
],
    Inches(0.25), Inches(1.75), Inches(3.8), Inches(5.0), sz=13, color=DGRAY)

fig(sl, FIGURES["meta_forest"], Inches(4.3), Inches(1.2), Inches(8.8), Inches(5.85))


# ─── SLIDE 4: BIDIRECTIONAL MR ───────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Bidirectional MR: Chol → Ca, Not the Reverse",
       "Two-sample MR, UK Biobank instruments → Michigan Genomics Initiative (n = 102,441)")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(5.8), Inches(5.85), fill=WHITE)
txt(sl, "Results", Inches(0.25), Inches(1.3), Inches(5.6), Inches(0.4),
    sz=15, bold=True, color=BLUE)

rows = [
    ("FORWARD  (Cholesterol → Calcium)", True, BLUE),
    ("Total Chol → Ca:   β = 0.065,   p = 5.8×10⁻⁴   ✓", False, DGRAY),
    ("LDL-C → Ca:   β = 0.053,   p = 0.006   ✓", False, DGRAY),
    ("Robust: MR-Egger, WM, PRESSO, MR-RAPS", False, RGBColor(0x55, 0x55, 0x55)),
    ("", False, DGRAY),
    ("REVERSE  (Calcium → Cholesterol)", True, RGBColor(0x99, 0x00, 0x00)),
    ("Ca → Total Chol:   β = 0.032,   p = 0.39   ✗  NULL", False, DGRAY),
    ("Ca → LDL-C:   β = −0.001,   p = 0.98   ✗  NULL", False, DGRAY),
    ("", False, DGRAY),
    ("Instruments:", True, BLUE),
    ("LDL-C: 236 SNPs, R² = 0.089, F = 173", False, DGRAY),
    ("Total Chol: 285 SNPs, R² = 0.097, F = 159", False, DGRAY),
    ("Calcium: 277 SNPs, R² = 0.064, F = 95", False, DGRAY),
    ("", False, DGRAY),
    ("Bayesian model comparison:", True, BLUE),
    ("96.5% posterior for Chol → Ca  (unidirectional)", False, RGBColor(0x00, 0x66, 0x00)),
]
y = Inches(1.75)
for text, bold, color in rows:
    if text:
        txt(sl, text, Inches(0.3), y, Inches(5.5), Inches(0.3),
            sz=12.5, bold=bold, color=color)
    y += Inches(0.285)

fig(sl, FIGURES["mr_summary"],   Inches(6.1), Inches(1.2), Inches(3.7), Inches(5.85))
fig(sl, FIGURES["bayesian_dir"], Inches(9.9), Inches(1.2), Inches(3.2), Inches(5.85))


# ─── SLIDE 5: DRUG-TARGET MR ─────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Drug-Target MR: HMGCR/Mevalonate Drives Both BMD Loss and Calcium Rise",
       "cis-SNPs for HMGCR, PCSK9, NPC1L1 as instruments; heel eBMD (Morris 2019, n=426k) + serum Ca (MGI, UKB)")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

# Key annotation strip above figures
ann = [
    ("HMGCR → BMD: β=−0.115, p<10⁻¹²  ✓", RGBColor(0x00, 0x66, 0x00)),
    ("PCSK9 → BMD: β=−0.014, p=0.61  ✗", RGBColor(0x99, 0x00, 0x00)),
    ("HMGCR → Ca:  β=+0.185, p<10⁻²⁰  ✓", RGBColor(0x00, 0x66, 0x00)),
    ("PCSK9 → Ca:  β=+0.004, p=0.87  ✗", RGBColor(0x99, 0x00, 0x00)),
]
for ai, (atxt, acol) in enumerate(ann):
    ax = Inches(0.3 + ai * 3.2)
    rect(sl, ax, Inches(1.28), Inches(3.05), Inches(0.32), fill=acol)
    txt(sl, atxt, ax + Inches(0.05), Inches(1.3), Inches(3.0), Inches(0.26),
        sz=11.5, bold=True, color=WHITE)

txt(sl,
    "Same mevalonate fingerprint for both outcomes:  HMGCR (statin target) drives BMD ↓ and Ca ↑; "
    "PCSK9i (mevalonate-independent) is null for both.  Effect sizes 3–5× larger than genome-wide MR.",
    Inches(0.3), Inches(1.65), Inches(12.7), Inches(0.42),
    sz=12, bold=True, color=BLUE, align=PP_ALIGN.CENTER)

# Two-panel HMGCR/PCSK9/NPC1L1 → Ca and BMD figure (full width)
fig(sl, FIGURES["drug_target_both"], Inches(0.2), Inches(2.12), Inches(12.9), Inches(4.92))


# ─── SLIDE 6: BMD MECHANISM ──────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Cholesterol Lowers BMD — But Not Through Vitamin D",
       "Drug-target & downstream MR; heel BMD (Morris 2019, n=426,824) and 25-OH Vit D (MGI-BioVU)")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(5.3), Inches(5.85), fill=WHITE)

# Hypothesis 1 — Vitamin D route (REFUTED)
rect(sl, Inches(0.25), Inches(1.28), Inches(5.1), Inches(0.34), fill=RGBColor(0xAA, 0x00, 0x00))
txt(sl, "Hypothesis 1: Chol → ↑ Vitamin D → bone effect  ✗  REFUTED",
    Inches(0.3), Inches(1.3), Inches(4.95), Inches(0.28),
    sz=11, bold=True, color=WHITE)
rect(sl, Inches(0.25), Inches(1.62), Inches(5.1), Inches(1.48),
     fill=RGBColor(0xFF, 0xEE, 0xEE), line=RGBColor(0xCC, 0x88, 0x88))
bullets(sl, [
    "Cholesterol is a Vit D precursor → hypothesis: higher",
    "  Chol → higher 25-OH Vit D → better bone",
    "MR: Total Chol → 25-OH Vit D (Revez, UKB n ≈ 418k):",
    "  β = −0.137,  p = 10⁻²⁹  ← INVERSE sign",
    "  Higher cholesterol → lower vitamin D",
    "Bayesian model weight: 2.7% for VitD route",
],
    Inches(0.3), Inches(1.66), Inches(4.9), Inches(1.4), sz=11,
    color=RGBColor(0x77, 0x00, 0x00))

# Hypothesis 2 — BMD route (SUPPORTED)
rect(sl, Inches(0.25), Inches(3.18), Inches(5.1), Inches(0.34), fill=RGBColor(0x00, 0x66, 0x00))
txt(sl, "Hypothesis 2: Chol → ↓ BMD  ✓  WELL SUPPORTED",
    Inches(0.3), Inches(3.2), Inches(4.95), Inches(0.28),
    sz=11, bold=True, color=WHITE)
rect(sl, Inches(0.25), Inches(3.52), Inches(5.1), Inches(2.52),
     fill=RGBColor(0xEE, 0xFF, 0xEE), line=RGBColor(0x88, 0xCC, 0x88))
bullets(sl, [
    "Genome-wide Chol → Heel eBMD (Morris 2019, n = 426,824):",
    "  β = −0.051,  p = 1.4×10⁻⁴  (280 SNPs)",
    "  WM: −0.026;  WMode: −0.031;  PRESSO: −0.047",
    "",
    "HMGCR drug-target → Heel eBMD:",
    "  β = −0.115,  p < 10⁻¹²  (37 cis-SNPs)",
    "  ~2× larger than genome-wide estimate",
    "",
    "PCSK9 → Heel eBMD:  β = −0.014,  p = 0.61  ✗ null",
    "  Mevalonate-specific fingerprint (same as Ca)",
    "",
    "Fractures (Dönertaş, n ≈ 484k):  β ≈ 0.0002,  p = 0.73",
    "  BMD effect without fracture effect = density, not quality",
    "",
    "Bayesian model weight: 97.3% for BMD route",
],
    Inches(0.3), Inches(3.56), Inches(4.9), Inches(2.44), sz=11,
    color=RGBColor(0x00, 0x44, 0x00))

fig(sl, FIGURES["bmd_forest"], Inches(5.6), Inches(1.2), Inches(4.0), Inches(5.85))
fig(sl, FIGURES["bayes_bar"],  Inches(9.65), Inches(1.2), Inches(3.5), Inches(5.85))


# ─── SLIDE 7: MEDIATION IS BROKEN ────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Bone Loss and Calcium Rise Are Parallel — Not Serial — Effects",
       "MVMR mediation + mediator screen: BMD → Ca is null at every skeletal site")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

# Refuted hypothesis
txt(sl, "Original Hypothesis (REFUTED):", Inches(0.4), Inches(1.35), Inches(5), Inches(0.4),
    sz=14, bold=True, color=RGBColor(0xAA, 0x00, 0x00))
rect(sl, Inches(0.4), Inches(1.8), Inches(12.1), Inches(1.0),
     fill=RGBColor(0xFF, 0xEE, 0xEE), line=RGBColor(0xCC, 0x88, 0x88))
txt(sl, "↑ Cholesterol  →  ↓ BMD  →  ↑ Serum Ca   (bone as the source of elevated calcium)",
    Inches(0.6), Inches(1.92), Inches(11.8), Inches(0.45),
    sz=15, color=RGBColor(0x88, 0x00, 0x00))
txt(sl, "✗  Mediation link 2 (BMD → Ca) is NULL at every skeletal site",
    Inches(0.6), Inches(2.35), Inches(11.8), Inches(0.35),
    sz=13, italic=True, color=RGBColor(0x88, 0x00, 0x00))

# Current model
txt(sl, "Current Model (SUPPORTED):", Inches(0.4), Inches(3.0), Inches(5), Inches(0.4),
    sz=14, bold=True, color=RGBColor(0x00, 0x66, 0x00))
rect(sl, Inches(0.4), Inches(3.45), Inches(12.1), Inches(1.1),
     fill=RGBColor(0xEE, 0xFF, 0xEE), line=RGBColor(0x88, 0xCC, 0x88))
txt(sl, "↑ HMGCR activity  →  ↓ BMD  (parallel)\n"
       "↑ HMGCR activity  →  ↑ Serum Ca  (parallel; calcium-specific mechanism — probably renal/tubular)",
    Inches(0.6), Inches(3.55), Inches(11.8), Inches(0.88),
    sz=14, color=RGBColor(0x00, 0x55, 0x00))

# Evidence table
txt(sl, "Mediation Evidence:", Inches(0.4), Inches(4.7), Inches(4.5), Inches(0.35),
    sz=14, bold=True, color=BLUE)

rows2 = [
    ("Heel eBMD → Ca",           "β = −0.009,  p = 0.45",   "301 SNPs, F = 197"),
    ("Femoral neck BMD → Ca",    "β = −0.001,  p = 0.97",   "—"),
    ("Total body BMD → Ca",      "β = −0.012,  p = 0.42",   "—"),
    ("HMGCR via BMD (mediation)","~0.5% of total effect",    "Negligible"),
    ("HMGCR → Phosphate (PO₄)",  "β = 0.027,  p = 0.15  (FLAT)", "Rules out PTH, VitD, bone resorption"),
    ("HMGCR → Ca | albumin",     "112% of total effect",     "Not an albumin artefact"),
]
y2 = Inches(5.1)
for r2 in rows2:
    rect(sl, Inches(0.4), y2, Inches(3.8), Inches(0.26), fill=RGBColor(0xEE, 0xF4, 0xFF))
    txt(sl, r2[0], Inches(0.45), y2 + Inches(0.03), Inches(3.7), Inches(0.2),
        sz=11, bold=True, color=BLUE)
    txt(sl, r2[1], Inches(4.3),  y2 + Inches(0.03), Inches(4.3), Inches(0.2), sz=11, color=DGRAY)
    txt(sl, r2[2], Inches(8.65), y2 + Inches(0.03), Inches(4.6), Inches(0.2),
        sz=10, italic=True, color=RGBColor(0x55, 0x55, 0x55))
    y2 += Inches(0.29)


# ─── SLIDE 8: ROBUST MR ──────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Robust MR: Cholesterol → BMD Survives Overlap & Pleiotropy Corrections",
       "5 estimators × 2 exposure arms (overlapping GLGC 2021 + overlap-free GLGC 2013)")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

heads = ["Method", "@ 5×10⁻⁸", "@ 5×10⁻⁵", "Key feature"]
cw    = [2.2, 1.9, 1.9, 6.5]
cx    = [0.25, 2.5, 4.45, 6.4]
rect(sl, Inches(0.25), Inches(1.3), Inches(12.6), Inches(0.35), fill=BLUE)
for i, h in enumerate(heads):
    txt(sl, h, Inches(cx[i] + 0.05), Inches(1.32), Inches(cw[i] - 0.1), Inches(0.3),
        sz=12, bold=True, color=WHITE)

table_rows = [
    ("IVW (random effects)",        "−0.067", "—",     "Baseline; SE ~4× too narrow"),
    ("MR-RAPS (ℓ₂)",               "−0.084", "−0.056","Weak-IV corrected; balanced pleiotropy"),
    ("MR-RAPS (tukey robust)",      "−0.062", "−0.043","Robust loss; outlier-deletion pushes toward null"),
    ("MRBEE",                       "−0.072", "−0.051","Overlap (Rxy) + estimation-error bias correction"),
    ("MR-APSS",                     "—",      "−0.056","LDSC C-matrix + winner's-curse; overlap-free: −0.037"),
    ("CAUSE (γ)",                   "—",      "−0.040","Shared-factor model; model comparison inconclusive"),
    ("", "", "", ""),
    ("Sample overlap (MR-APSS)",    "0.004 (0.17 SE)", "—", "IMMATERIAL — measured, not assumed"),
    ("Pleiotropy SD (τ̂)",          "0.005 (~3× SE)",  "—", "Heavy but diffuse → random-effects treatment"),
]
for ri, row in enumerate(table_rows):
    bg = LGRAY if ri % 2 == 0 else WHITE
    if "overlap" in row[0].lower() or "Pleiotrop" in row[0]:
        bg = RGBColor(0xFF, 0xFF, 0xDD)
    rect(sl, Inches(0.25), Inches(1.65 + ri * 0.37), Inches(12.6), Inches(0.36), fill=bg)
    for ci2, val in enumerate(row):
        txt(sl, val, Inches(cx[ci2] + 0.05), Inches(1.65 + ri * 0.37 + 0.04),
            Inches(cw[ci2] - 0.1), Inches(0.28), sz=11.5, color=DGRAY)

txt(sl,
    "All 5 methods negative; honest interval ≈ −0.09 to −0.01.  "
    "CAUSE directionally consistent; model comparison inconclusive (low power).  "
    "Threshold drives more variation between methods than model assumptions do.",
    Inches(0.3), Inches(5.08), Inches(12.7), Inches(0.5),
    sz=12.5, bold=True, color=BLUE, align=PP_ALIGN.CENTER)

fig(sl, FIGURES["mechanisms"],  Inches(0.2),  Inches(5.65), Inches(6.2), Inches(1.5))
fig(sl, FIGURES["mech_primary"],Inches(6.5),  Inches(5.65), Inches(6.7), Inches(1.5))


# ─── SLIDE 9: CLINICAL — STATIN INDEPENDENCE ─────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Clinical Evidence: Statins Protect Bone Independent of LDL-C",
       "Michigan Genomics Initiative EHR; Fine-Gray competing risks (n ≈ 5,784; 595 osteoporosis events)")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(5.5), Inches(5.85), fill=WHITE)
txt(sl, "Primary Results (Osteoporosis Endpoint)", Inches(0.25), Inches(1.3), Inches(5.3), Inches(0.38),
    sz=14, bold=True, color=BLUE)

key_results = [
    ("LDL-C (per 10 mg/dL)", "Sub-HR = 0.985  (0.962–1.009)", "p = 0.21 — NULL", False),
    ("Ever statin use",       "Sub-HR = 0.650  (0.530–0.787)", "p < 0.001 — PROTECTIVE", True),
    ("Female sex",            "Sub-HR = 6.035  (4.941–7.372)", "p < 0.001", False),
    ("Age (per year)",        "Sub-HR = 1.055  (1.048–1.061)", "p < 0.001", False),
]
yr = Inches(1.72)
for label, val, p, hi in key_results:
    bg = RGBColor(0xEE, 0xFF, 0xEE) if hi else RGBColor(0xF5, 0xF8, 0xFF)
    rect(sl, Inches(0.25), yr, Inches(5.2), Inches(0.68), fill=bg)
    txt(sl, label, Inches(0.35), yr + Inches(0.04), Inches(5.0), Inches(0.27),
        sz=12.5, bold=True, color=BLUE)
    txt(sl, val + "  |  " + p, Inches(0.35), yr + Inches(0.33), Inches(5.0), Inches(0.28),
        sz=12, color=DGRAY)
    yr += Inches(0.74)

txt(sl, "Robust across 9 model specifications:", Inches(0.25), Inches(4.73), Inches(5.2), Inches(0.3),
    sz=13, bold=True, color=BLUE)
bullets(sl, [
    "Fine-Gray 3-state (MACE + death competing), baseline, LOCF",
    "Standard Cox, time-varying, statin-naive subgroup, all patients",
    "→ LDL-C null in all 9;  statin protective in primary model",
],
    Inches(0.25), Inches(5.05), Inches(5.2), Inches(0.9), sz=12, color=DGRAY)

rect(sl, Inches(0.25), Inches(6.1), Inches(5.2), Inches(0.75), fill=MAIZE)
txt(sl,
    "Statin bone-protection (Sub-HR 0.65) persists after full LDL-C adjustment —\n"
    "mevalonate pathway, not cholesterol lowering",
    Inches(0.3), Inches(6.2), Inches(5.1), Inches(0.55),
    sz=11.5, bold=True, color=BLUE, align=PP_ALIGN.CENTER)

fig(sl, FIGURES["statin_ldlc"], Inches(5.85), Inches(1.2), Inches(7.3), Inches(5.85))


# ─── SLIDE 10: PROPOSED MECHANISM ────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Proposed Mechanism: GGPP Depletion in Osteoclasts",
       "Mevalonate → GGPP → Rho GTPase prenylation → osteoclast-mediated bone demineralization")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

# Biochemistry framework figure (left 58%)
fig(sl, FIGURES["biochem"], Inches(0.3), Inches(1.3), Inches(7.4), Inches(4.95))

# Right-side evidence summary
txt(sl, "MR/EHR evidence summary", Inches(7.9), Inches(1.35), Inches(5.2), Inches(0.38),
    sz=13, bold=True, color=BLUE)

side = [
    ("HMGCR → Ca:  β = 0.185  (p < 10⁻²⁰)",              RGBColor(0x00, 0x66, 0x00)),
    ("HMGCR → BMD:  β = −0.115  (p < 10⁻¹²)",             RGBColor(0x00, 0x66, 0x00)),
    ("PCSK9 → Ca:  β = 0.004,  p = 0.87  ✗  NULL",        RGBColor(0xAA, 0x00, 0x00)),
    ("PCSK9 → BMD:  β = −0.014,  p = 0.61  ✗  NULL",      RGBColor(0xAA, 0x00, 0x00)),
    ("HMGCR eQTL (blood) → BMD:  β=−0.044,  p=4.4×10⁻⁵", MBLU),
    ("Statin (EHR):  Sub-HR = 0.65,  LDL-C independent",   MBLU),
]
ys = Inches(1.8)
for s_txt, s_col in side:
    rect(sl, Inches(7.85), ys, Inches(5.25), Inches(0.42), fill=RGBColor(0xF5, 0xF8, 0xFF))
    txt(sl, s_txt, Inches(7.95), ys + Inches(0.06), Inches(5.1), Inches(0.32),
        sz=11.5, color=s_col)
    ys += Inches(0.48)

txt(sl, "Key implication:", Inches(7.85), Inches(4.75), Inches(5.2), Inches(0.32),
    sz=12, bold=True, color=BLUE)
rect(sl, Inches(7.85), Inches(5.1), Inches(5.25), Inches(1.0),
     fill=RGBColor(0xFF, 0xFB, 0xEE), line=RGBColor(0xCC, 0xAA, 0x00))
txt(sl,
    "Statins: block HMGCR → deplete GGPP → impair\n"
    "Rho GTPase prenylation in osteoclasts.\n"
    "Ezetimibe / PCSK9i: bypass mevalonate → no bone benefit.\n"
    "N-bisphosphonates: block FPP synthase (same pathway).\n"
    "→  Aim 2 will test GGPP with Ggps1 KO + GGTI-298.",
    Inches(7.95), Inches(5.15), Inches(5.1), Inches(0.9),
    sz=11, color=DGRAY)


# ─── SLIDE 11: PROPOSED AIMS ─────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Proposed R01 Aims",
       "Hypothesis: mevalonate flux promotes bone demineralization via GGPP depletion in osteoclasts")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

rect(sl, Inches(0.25), Inches(1.3), Inches(12.6), Inches(0.4), fill=BLUE)
txt(sl, "Aim 1:  Define Whether Statin-Mediated Bone Protection Is Mevalonate-Dependent",
    Inches(0.35), Inches(1.33), Inches(12.4), Inches(0.34), sz=14, bold=True, color=WHITE)

aim1 = [
    ("1A", "WT mice (both sexes, NCD + HFHS): pravastatin vs. evolocumab (PCSK9i) vs. ezetimibe "
           "at equivalent LDL-C reductions.  Outcomes: bone density, architecture (μCT), strength, turnover."),
    ("1B", "Osteoclast-specific Hmgcr knockout mice:  defines cell-type specificity of the statin bone effect."),
    ("1C", "MVMR (human):  formally test HMGCR → BMD → Ca causal chain vs. parallel effects in UKB/MGI."),
    ("1D", "EHR (All of Us):  statin vs. ezetimibe initiators → incident fracture and BMD (active-comparator design)."),
]
ya = Inches(1.75)
for code, t in aim1:
    rect(sl, Inches(0.35), ya, Inches(0.5), Inches(0.4), fill=MBLU)
    txt(sl, code, Inches(0.38), ya + Inches(0.07), Inches(0.44), Inches(0.26),
        sz=12, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    txt(sl, t, Inches(0.9), ya + Inches(0.06), Inches(12.0), Inches(0.4), sz=12, color=DGRAY)
    ya += Inches(0.47)

rect(sl, Inches(0.25), Inches(3.85), Inches(12.6), Inches(0.4), fill=RGBColor(0x77, 0x44, 0x00))
txt(sl, "Aim 2:  Establish That Geranylgeranylation Mediates Statin-Mediated Bone Protection",
    Inches(0.35), Inches(3.88), Inches(12.4), Inches(0.34), sz=14, bold=True, color=WHITE)

aim2 = [
    ("2A", "Primary osteoclasts from Aim 1 mice: GGPP supplementation rescue test — "
           "does GGPP restore statin/bisphosphonate-inhibited osteoclast differentiation?  RhoA/Rac1 prenylation readout."),
    ("2B", "Osteoclast-specific Ggps1 KO (GGPP synthase):  in vivo ablation of GGPP in the osteoclast lineage."),
    ("2C", "Pharmacological validation:  GGPP synthase inhibitor (GGTI-298) in HFHS diet mice.  "
           "Primary outcomes: bone density, strength, architecture (μCT)."),
]
ya = Inches(4.3)
for code, t in aim2:
    rect(sl, Inches(0.35), ya, Inches(0.5), Inches(0.4), fill=RGBColor(0xCC, 0x77, 0x00))
    txt(sl, code, Inches(0.38), ya + Inches(0.07), Inches(0.44), Inches(0.26),
        sz=12, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    txt(sl, t, Inches(0.9), ya + Inches(0.06), Inches(12.0), Inches(0.4), sz=12, color=DGRAY)
    ya += Inches(0.47)

txt(sl,
    "Expected outcome: statins (not PCSK9i, not ezetimibe) are bone-protective via mevalonate/GGPP; "
    "osteoclast-specific Ggps1 KO phenocopies statin treatment on bone.",
    Inches(0.35), Inches(6.62), Inches(12.6), Inches(0.38),
    sz=12.5, bold=True, color=BLUE, align=PP_ALIGN.CENTER)


# ─── SLIDE 12: OPEN QUESTIONS / COLLABORATION ────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
header(sl, "Open Questions and Where We'd Value Your Expertise",
       "Areas where in vivo / cell bone biology expertise is critical")
footer(sl)

rect(sl, Inches(0.15), Inches(1.2), Inches(12.9), Inches(5.85), fill=WHITE)

rect(sl, Inches(0.25), Inches(1.3), Inches(5.8), Inches(5.65), fill=LGRAY)
txt(sl, "What We Have", Inches(0.35), Inches(1.35), Inches(5.6), Inches(0.38),
    sz=15, bold=True, color=BLUE)
bullets(sl, [
    "Strong human MR + epidemiological evidence base",
    "HMGCR/mevalonate pathway specificity (vs PCSK9, NPC1L1)",
    "EHR clinical cohort (MGI, ~20,000 patients)",
    "Bioinformatics / statistical genetics expertise",
    "",
    "Planned analyses (near-term):",
    "  Update meta-analysis (new studies identified)",
    "  MRBEE/CAUSE bidirectional MR in DO mice",
    "    (test Chol↔Ca direction in the mouse data)",
    "  Renal test: HMGCR → urinary Ca excretion",
    "  PTH via deCODE pQTL (Ferkingstad 2021, n≈35k)",
    "  All of Us: statin vs ezetimibe fracture analysis",
],
    Inches(0.35), Inches(1.78), Inches(5.6), Inches(5.0), sz=12.5, color=DGRAY)

rect(sl, Inches(6.3), Inches(1.3), Inches(6.6), Inches(5.65), fill=RGBColor(0xFF, 0xFB, 0xEE))
txt(sl, "Where We'd Value Your Input", Inches(6.4), Inches(1.35), Inches(6.4), Inches(0.38),
    sz=15, bold=True, color=RGBColor(0x77, 0x44, 0x00))
bullets(sl, [
    "Osteoclast model systems:",
    "  Best primary osteoclast protocols for statin/GGPP studies?",
    "  Cell lines vs. primary for differentiation assays?",
    "  RhoA/Rac1 prenylation readout approaches?",
    "",
    "In vivo bone phenotyping:",
    "  Bone strength testing (3-pt bending vs. compression)?",
    "  Turnover markers: plasma CTX / P1NP / osteocalcin?",
    "  μCT parameters most sensitive to osteoclast activity?",
    "",
    "Conditional KO models:",
    "  Preferred Cre driver: Ctsk-Cre vs. LysM-Cre?",
    "  Known osteoblast 'leak' issues to be aware of?",
    "",
    "Visualizations still needed (flagged):",
    "  Statin vs. ezetimibe in vivo bone phenotype schematic",
    "  Clinical trial fracture aggregate figure",
],
    Inches(6.4), Inches(1.78), Inches(6.4), Inches(5.0), sz=12, color=DGRAY)


# ─── SLIDE 13: SUMMARY ───────────────────────────────────────────────────────
sl = prs.slides.add_slide(blank)
rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=BLUE)
rect(sl, 0, Inches(1.65), SLIDE_W, Inches(0.06), fill=MAIZE)
rect(sl, 0, Inches(6.8), SLIDE_W, Inches(0.7), fill=RGBColor(0x00, 0x17, 0x30))

txt(sl, "Summary", Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.7),
    sz=34, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
txt(sl, "Seven converging lines of evidence for a causal cholesterol → bone demineralization pathway",
    Inches(0.5), Inches(0.95), Inches(12.3), Inches(0.55),
    sz=15, italic=True, color=MAIZE, align=PP_ALIGN.CENTER)

summary = [
    "1.  DO mice + human meta-analysis: consistent Ca–Chol correlation (r = 0.176, p = 4.2×10⁻⁸; Spearman ρ ≈ 0.4–0.5 in mice)",
    "2.  Bidirectional MR: Chol → Ca (p < 0.001); Ca → Chol = NULL.  96.5% Bayesian posterior for unidirectional direction",
    "3.  Pathway specificity: HMGCR β ≈ 0.185 (3–5× genome-wide); PCSK9 β = 0.004 (F > 2000, null) → mevalonate, not LDL-C",
    "4.  Bone as intermediate: HMGCR → BMD β = −0.115 (p < 10⁻¹²); PCSK9 → BMD = null.  Identical mevalonate fingerprint",
    "5.  Serum Ca is homeostatically defended: BMD → Ca null at all 3 skeletal sites.  Bone loss and Ca rise are parallel effects",
    "6.  Robust to pleiotropy/overlap: 5 methods × 2 exposure arms all negative.  Honest CI: β ≈ −0.09 to −0.01",
    "7.  Clinical: statins protect bone (Sub-HR 0.65) independent of LDL-C in MGI EHR — mevalonate mechanism confirmed",
]
ys = Inches(1.8)
for i, s in enumerate(summary):
    rect(sl, Inches(0.3), ys, Inches(0.42), Inches(0.42), fill=MAIZE)
    txt(sl, str(i + 1) + ".", Inches(0.32), ys + Inches(0.05), Inches(0.38), Inches(0.32),
        sz=13, bold=True, color=BLUE, align=PP_ALIGN.CENTER)
    txt(sl, s[3:], Inches(0.8), ys + Inches(0.03), Inches(12.3), Inches(0.37),
        sz=12.5, color=WHITE)
    ys += Inches(0.52)

txt(sl,
    "Hypothesis: mevalonate flux → GGPP depletion in osteoclasts → bone demineralization + Ca rise.  "
    "Statins (not PCSK9i, not ezetimibe) block this.  Bisphosphonates act downstream at FPP synthase.",
    Inches(0.5), Inches(6.1), Inches(12.3), Inches(0.52),
    sz=13, bold=True, color=MAIZE, align=PP_ALIGN.CENTER)

txt(sl, "Dave Bridges Lab  ·  davebrid@umich.edu  ·  University of Michigan",
    Inches(0.5), Inches(6.85), Inches(12.3), Inches(0.28),
    sz=10, color=RGBColor(0xAA, 0xAA, 0xAA), align=PP_ALIGN.CENTER)


# ─── SAVE ────────────────────────────────────────────────────────────────────
out = "/Users/davebrid/Documents/GitHub/PrecisionNutrition/Presentations/2026-11-13 Surface Lab Update.pptx"
prs.save(out)
print(f"Saved: {out}  ({len(prs.slides)} slides)")
