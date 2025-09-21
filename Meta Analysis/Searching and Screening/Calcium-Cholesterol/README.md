## Initial Screening of Manuscripts

Following the search strategy described in the [PROSPERO documentation](../../PROSPERO%20Registration/Calcium-Cholesterol/search_strategy.md).  Searches were done on PubMed, EMBASE, Cochrane (CENTRAL, and CINAHL), and Web of Science with bib files saved. These were imported into Zotero (3765 total references) and de-duplicated to 2195 references.

### Raw Files

These are all found in this folder

- PubMed: pubmed-serumcalci-set.nbib (583 records)
- EMBASE: records.ris (1209 records)
- Cochrane CENTRAL: citation-export.txt (218 records)
- Cochrane CINAHL: ae563724-ae52-4a0f-b528-4bc02b1205ad.bib (734 records).  Had to manually clean double quotes and extra curly brackets in the abstracts.
- Web of Science: savedrecs.ris and savedrecs (1).ris (1021 records, had to download in two parts) https://www.webofscience.com/wos/woscc/summary/55f88581-4004-40e7-9bdb-1cfb392fbddf-0179f0d203/relevance/1

### Initial GPT-Based Screening

We used a scripted GPT-5 prompt (Supplementary Table S1) to classify and exclude non-eligible records based on Title and Abstract. Prompt text and code are available in our repository.

#### Prompt

You are acting as an expert systematic-reviewer following PRISMA guidelines.  
You will classify each article in a deduplicated list by reading its Title and Abstract (fields: `Title`, `Abstract Note`).

Flag any article as:
- "Case report/series" if the title or abstract clearly states “case report”, “case series”, or “a case of …”.
- "Review/Systematic review/Meta-analysis" if the title or abstract clearly states “systematic review”, “meta-analysis”, “review article”, or “narrative review”.
- "Non-human study" only if it is clear the study population is not human and no human data are included.
- "Unsure" if you cannot confidently classify it into one of the above categories.

Be conservative: only flag an article if it clearly matches the criteria.  
Do not exclude any article with both human and non-human data.

Output a CSV file with two columns:  
`Key` (article key from input file) and `Reason` (one of the four categories above).

#### GPT Screening Validation

## Notes:

- DB indicates Dave Bridges under screener
- [] around a paper indicates it was translated by PubMed from a different intial language
- If a paper is excluded for multiple reasons, they are separated by a ";"

## Deep Screening of Manuscripts