"""Supplementary-table parser used in v0.13.

The runtime's spreadsheet RPC failed while importing the supplied XLSX files.
To avoid changing spreadsheet libraries or silently converting content, this
release reads the XLSX Open Packaging Convention directly:

- sharedStrings.xml
- workbook.xml / workbook relationships
- worksheet sheetData cells

No formulas are evaluated and no workbook is modified.

Outputs are derivative CSV/JSON crosswalks only.
"""
