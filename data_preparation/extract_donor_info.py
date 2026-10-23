"""Extract donor characteristics into a spreadsheet.

Usage: python extract_donor_info.py [pdf_dir] [output.xlsx]
"""

import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import pandas as pd
import pdfplumber

PDF_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/donor_reports")
OUT_XLSX = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("output/donors_extracted.xlsx")


def pdf_to_text(pdf_path: Path) -> str:
    with pdfplumber.open(str(pdf_path)) as pdf:
        return "\n".join(page.extract_text() or "" for page in pdf.pages)


def first_match(text: str, pattern: str) -> Optional[str]:
    m = re.search(pattern, text, re.IGNORECASE)
    # value is always the last capture group (some patterns capture a label first)
    return m.group(m.lastindex).strip() if m else None


def to_float(x):
    try:
        return float(x.replace(",", "")) if x else None
    except ValueError:
        return None


def to_int(x):
    try:
        return int(float(x)) if x else None
    except ValueError:
        return None


def yn_to_10(x):
    if not x:
        return None
    x = x.lower()
    if x in ("yes", "y", "positive", "pos"):
        return 1
    if x in ("no", "n", "negative", "neg"):
        return 0
    return None


def norm_sex(x):
    if not x:
        return None
    x = x.upper()
    if x in ("M", "MALE"):
        return "M"
    if x in ("F", "FEMALE"):
        return "F"
    return x


def donor_id_from_filename(pdf_path: Path):
    nums = re.findall(r"\d+", pdf_path.stem)
    if not nums:
        return None
    return to_int(max(nums, key=len))


def extract_fields(text: str) -> Dict[str, Any]:
    return {
        "donor_id": None,
        "donor_age": to_int(first_match(text, r"\bAge\b\s*:?\s*(\d+)")),
        "sex_FM": norm_sex(first_match(text, r"\bSex\b\s*:?\s*(Male|Female|M|F)")),
        "donor_type": first_match(text, r"Donor\s*Type\s*:?\s*([A-Za-z0-9\/\- ]+)"),
        "donor_cod": first_match(text, r"Cause\s+of\s+Death\s*:?\s*(.+)"),
        "donor_weight": to_float(first_match(text, r"Weight\s*:?\s*(\d+\.?\d*)")),
        "donor_height": to_float(first_match(text, r"Height\s*:?\s*(\d+\.?\d*)")),
        "donor_bmi": to_float(first_match(text, r"BMI\s*:?\s*(\d+\.?\d*)")),
        "donor_los_days": to_int(first_match(text, r"Length\s+of\s+stay\s*:?\s*(\d+)")),
        "donor_ooha": yn_to_10(first_match(text, r"(OOHA|Out\s*of\s*hospital\s*arrest)\s*:?\s*(Yes|No)")),
        "donor_downtime": to_float(first_match(text, r"Downtime\s*:?\s*(\d+\.?\d*)")),
        "donor_htn": yn_to_10(first_match(text, r"Hypertension\s*:?\s*(Yes|No)")),
        "donor_cvd": yn_to_10(first_match(text, r"(CVD|Cardiovascular\s+disease)\s*:?\s*(Yes|No)")),
        "donor_diabetes": yn_to_10(first_match(text, r"Diabetes\s*:?\s*(Yes|No)")),
        "donor_smoking": yn_to_10(first_match(text, r"(Smoking|Smoker)\s*:?\s*(Yes|No)")),
        "donor_cmv": yn_to_10(first_match(text, r"CMV\s*:?\s*(Positive|Negative|Yes|No)")),
        "donor_creat_adm": to_float(first_match(text, r"Admission\s+creatinine\s*:?\s*(\d+\.?\d*)")),
        "donor_egfr_adm": to_float(first_match(text, r"Admission\s+eGFR\s*:?\s*(\d+\.?\d*)")),
        "donor_creat_last": to_float(first_match(text, r"Last\s+creatinine\s*:?\s*(\d+\.?\d*)")),
        "donor_egfr_last": to_float(first_match(text, r"(Last\s+eGFR|GFR)\s*:?\s*(\d+\.?\d*)")),
        "donor_abo": first_match(text, r"(ABO|Blood\s+group)\s*:?\s*(AB|A|B|O)"),
        "donor_rh": first_match(text, r"(Rh|Rhesus)\s*:?\s*([+\-]|Positive|Negative)"),
    }


def main():
    if not PDF_DIR.exists():
        sys.exit(f"PDF folder not found: {PDF_DIR}")

    rows = []
    for pdf in sorted(PDF_DIR.glob("*.pdf")):
        print("Extracting:", pdf.name)
        data = extract_fields(pdf_to_text(pdf))
        data["donor_id"] = donor_id_from_filename(pdf)
        rows.append(data)

    df = pd.DataFrame(rows).sort_values("donor_id").drop_duplicates("donor_id", keep="last")
    OUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
    df.to_excel(OUT_XLSX, index=False)
    print("Saved:", OUT_XLSX)


if __name__ == "__main__":
    main()
