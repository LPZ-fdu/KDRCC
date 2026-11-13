# Chinese methods and simulation appendix

`Appendix_KDRCC_zh.docx` is the submission document. `appendix_zh.json` stores its ordered sections, paragraphs, equations, tables and references. `equation_index.json` records the display-equation numbering.

## Organization

Sections S.1–S.9 define KDR-CC for a general fixed-dimensional exposure vector and finite-dimensional auxiliary models, and provide identification, estimation, covariance, effect-output and variable-testing theory. Sections S.10–S.11 specify the data-generating mechanism, auxiliary designs, A–E simulation modules and evaluation measures. Section S.12 records the numerical analysis configuration and its implementation mapping. The A035 record is defined there by its complete parameter values.

## Build

The builder uses Python 3, `python-docx`, `lxml` and Microsoft Office's `MML2OMML.XSL` stylesheet. The stylesheet is supplied by the installed Office distribution. The existing appendix is the formatting template. From this directory on a machine with Office installed:

```powershell
python build_appendix.py --xsl "C:/Program Files/Microsoft Office/root/Office16/MML2OMML.XSL"
```

The builder reads `appendix_zh.json`, uses the existing `Appendix_KDRCC_zh.docx` as its template and writes the document and equation index. Optional `--template`, `--source`, `--output` and `--index` arguments specify alternate paths.

Equations are editable Office Math objects. Each display occupies a one-row, two-column table with a centered equation and right-aligned section-local number. Ordinary paragraphs use a two-character first-line indent. Data tables use 1.0-point outer horizontal rules and a 0.75-point header separator. The template supplies fonts, paragraph styles, page geometry, headers and footers.

## Verification

Open the document in Microsoft Word, update page fields and inspect its rendered pages. The release was exported to PDF through Word and all 32 pages were visually inspected. The stage archive inventories record the document digest; `../docs/FINAL_REVIEW.md` records the methods review. The implementation definitions and simulation-specific settings are in `../docs/METHOD_SPECIFICATION.md`.
