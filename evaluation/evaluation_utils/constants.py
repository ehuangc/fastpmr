from pathlib import Path

EVALUATION_DIR = Path(__file__).resolve().parent.parent
FASTPMR_BIN = "fastpmr"
AADR_DIR = EVALUATION_DIR / "aadr"
AADR_DATA_PREFIX = AADR_DIR / "data" / "v66.1240K.aadr.PUB"
AADR_EXTS = (".anno", ".ind", ".snp", ".geno")
AADR_RUNS = 1
AADR_NPZ_PATH = AADR_DIR / "results" / "fastpmr" / "fastpmr_results.npz"
AADR_METADATA_PATH = AADR_DIR / "data" / "v66.1240K.aadr.PUB.anno"

CHICHEN_ITZA_DIR = EVALUATION_DIR / "chichen_itza"
CHICHEN_ITZA_DATA_DIR = CHICHEN_ITZA_DIR / "data"
# ENA accession with the submitted BAMs
CHICHEN_ITZA_ACCESSION = "PRJEB73567"
CHICHEN_ITZA_FILEREPORT = CHICHEN_ITZA_DATA_DIR / f"filereport_read_run_{CHICHEN_ITZA_ACCESSION}.json"
# ENA portal API query fetching the file report as JSON; submitted_ftp gives the BAM/BAI
# URLs and library_construction_protocol encodes single- vs. double-stranded (see run_eager.py)
CHICHEN_ITZA_FILEREPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    f"?accession={CHICHEN_ITZA_ACCESSION}&result=read_run"
    "&fields=submitted_ftp,library_construction_protocol&format=json"
)
CHICHEN_ITZA_BAM_DIR = CHICHEN_ITZA_DIR / "bams"
CHICHEN_ITZA_RENAMED_BAM_DIR = CHICHEN_ITZA_BAM_DIR / "renamed"
CHICHEN_ITZA_TSV = CHICHEN_ITZA_DIR / "config_eager_ych.tsv"
CHICHEN_ITZA_RESULTS_DIR = CHICHEN_ITZA_DIR / "results"
# Directory nf-core/eager's pileupcaller genotyping writes into (the per-strandedness and
# merged EIGENSTRAT prefixes are defined in run_fastpmr.py)
CHICHEN_ITZA_GENO_DIR = CHICHEN_ITZA_RESULTS_DIR / "genotyping"

# AADR metadata field names
LAT_FIELD = "Latitude"
LON_FIELD = "Longitude"
LOCALITY_FIELD = "Locality"
INDIVIDUAL_ID_FIELD = "Individual ID"
GROUP_ID_FIELD = "Group ID"
DATE_MEAN_BP_FIELD = (
    "Date mean in BP in years before 1950 CE "
    "[OxCal mu for a direct radiocarbon date, and average of range for a contextual date]"
)
FULL_DATE_FIELD = (
    "Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age "
    "(Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990+-40 BP, Ua-35016). "
    "(Format 2) Archaeological context range, e.g. 2500-1700 BCE"
)
PUBLICATION_FIELD = "Publication abbreviation"
SKELETAL_CODE_FIELD = "Skeletal code"
POLITICAL_ENTITY_FIELD = "Political Entity"

# fmt: off
COUNTRY_TO_REGION: dict[str, str] = {
    # ── Africa ──
    "AO": "africa", "BF": "africa", "BI": "africa", "BJ": "africa", "BW": "africa", "CD": "africa", "CF": "africa",
    "CG": "africa", "CI": "africa", "CM": "africa", "CV": "africa", "DJ": "africa", "DZ": "africa", "EG": "africa",
    "EH": "africa", "ER": "africa", "ET": "africa", "GA": "africa", "GH": "africa", "GM": "africa", "GN": "africa",
    "GQ": "africa", "GW": "africa", "KE": "africa", "KM": "africa", "LR": "africa", "LS": "africa", "LY": "africa",
    "MA": "africa", "MG": "africa", "ML": "africa", "MR": "africa", "MU": "africa", "MW": "africa", "MZ": "africa",
    "NA": "africa", "NE": "africa", "NG": "africa", "RE": "africa", "RW": "africa", "SC": "africa", "SD": "africa",
    "SH": "africa", "SL": "africa", "SN": "africa", "SO": "africa", "SS": "africa", "ST": "africa", "SZ": "africa",
    "TD": "africa", "TG": "africa", "TN": "africa", "TZ": "africa", "UG": "africa", "YT": "africa", "ZA": "africa",
    "ZM": "africa", "ZW": "africa",
    # ── Asia ──
    "AE": "asia", "AF": "asia", "AM": "asia", "AZ": "asia", "BD": "asia", "BH": "asia", "BN": "asia",
    "BT": "asia", "CN": "asia", "CY": "asia", "GE": "asia", "HK": "asia", "ID": "asia", "IL": "asia",
    "IN": "asia", "IQ": "asia", "IR": "asia", "JO": "asia", "JP": "asia", "KG": "asia", "KH": "asia",
    "KP": "asia", "KR": "asia", "KW": "asia", "KZ": "asia", "LA": "asia", "LB": "asia", "LK": "asia",
    "MM": "asia", "MN": "asia", "MO": "asia", "MV": "asia", "MY": "asia", "NP": "asia", "OM": "asia",
    "PH": "asia", "PK": "asia", "PS": "asia", "QA": "asia", "RU": "asia", "SA": "asia", "SG": "asia",
    "SY": "asia", "TH": "asia", "TJ": "asia", "TL": "asia", "TM": "asia", "TR": "asia", "TW": "asia",
    "UZ": "asia", "VN": "asia", "YE": "asia",
    # ── Europe ──
    "AD": "europe", "AL": "europe", "AT": "europe", "AX": "europe", "BA": "europe", "BE": "europe",
    "BG": "europe", "BY": "europe", "CH": "europe", "CZ": "europe", "DE": "europe", "DK": "europe",
    "EE": "europe", "ES": "europe", "FI": "europe", "FO": "europe", "FR": "europe", "GB": "europe",
    "GG": "europe", "GI": "europe", "GR": "europe", "HR": "europe", "HU": "europe", "IE": "europe",
    "IM": "europe", "IS": "europe", "IT": "europe", "JE": "europe", "LI": "europe", "LT": "europe",
    "LU": "europe", "LV": "europe", "MC": "europe", "MD": "europe", "ME": "europe", "MK": "europe",
    "MT": "europe", "NL": "europe", "NO": "europe", "PL": "europe", "PT": "europe", "RO": "europe",
    "RS": "europe", "SE": "europe", "SI": "europe", "SJ": "europe", "SK": "europe", "SM": "europe",
    "UA": "europe", "VA": "europe", "XK": "europe",
    # ── Oceania ──
    "AS": "oceania", "AU": "oceania", "CK": "oceania", "FJ": "oceania", "FM": "oceania", "GU": "oceania",
    "KI": "oceania", "MH": "oceania", "MP": "oceania", "NC": "oceania", "NF": "oceania", "NR": "oceania",
    "NU": "oceania", "NZ": "oceania", "PF": "oceania", "PG": "oceania", "PN": "oceania", "PW": "oceania",
    "SB": "oceania", "TK": "oceania", "TO": "oceania", "TV": "oceania", "VU": "oceania", "WF": "oceania",
    "WS": "oceania",
    # ── Americas ──
    "AG": "americas", "AI": "americas", "AR": "americas", "AW": "americas", "BB": "americas", "BL": "americas",
    "BM": "americas", "BO": "americas", "BQ": "americas", "BR": "americas", "BS": "americas", "BZ": "americas",
    "CA": "americas", "CL": "americas", "CO": "americas", "CR": "americas", "CU": "americas", "CW": "americas",
    "DM": "americas", "DO": "americas", "EC": "americas", "FK": "americas", "GD": "americas", "GF": "americas",
    "GL": "americas", "GP": "americas", "GT": "americas", "GY": "americas", "HN": "americas", "HT": "americas",
    "JM": "americas", "KN": "americas", "KY": "americas", "LC": "americas", "MF": "americas", "MQ": "americas",
    "MS": "americas", "MX": "americas", "NI": "americas", "PA": "americas", "PE": "americas", "PM": "americas",
    "PR": "americas", "PY": "americas", "SR": "americas", "SV": "americas", "SX": "americas", "TC": "americas",
    "TT": "americas", "US": "americas", "UY": "americas", "VC": "americas", "VE": "americas", "VG": "americas",
    "VI": "americas",
    # ── Miscellaneous territories ──
    "CC": "oceania", "CX": "oceania",  # Cocos/Christmas Islands (AU territory)
    "HM": "oceania",  # Heard/McDonald Islands (AU territory)
    "IO": "asia",  # British Indian Ocean Territory (nearest path via Cairo→India)
    "TF": "africa",  # French Southern Territories (least-wrong among five regions)
    "UM": "oceania",  # US Minor Outlying Islands
}
# fmt: on
EURASIA_REGIONS = {"asia", "europe"}
# US Pacific territories that reverse_geocoder reports as "US" but should classify as Oceania
US_OCEANIA_ADMIN1 = {"Guam", "Northern Mariana Islands"}

COMPARISON_DIR = EVALUATION_DIR / "comparison"
COMPARISON_DATA_PREFIX = COMPARISON_DIR / "data" / "maravall_lopez_et_al"

PERFORMANCE_DIR = EVALUATION_DIR / "performance"
PERFORMANCE_DATA_PREFIX = PERFORMANCE_DIR / "data" / "lazaridis_et_al"
PERFORMANCE_SAMPLE_SET_DIR = PERFORMANCE_DIR / "data" / "lazaridis_et_al_sample_sets"
PERFORMANCE_SAMPLE_SET_SIZES = (
    128,
    256,
    384,
    512,
    640,
    768,
    896,
    1024,
    1152,
    1280,
)
PERFORMANCE_RUNS = 3

EIGENSTRAT_EXTS = (".ind", ".snp", ".geno")
PLINK_EXTS = (".bed", ".bim", ".fam")
CSV_FIELDS = [
    "label",
    "mean_s",
    "stddev_s",
    "min_s",
    "max_s",
    "mean_bytes",
    "stddev_bytes",
    "min_bytes",
    "max_bytes",
]
