"""
Tests for classify_cells region classification.
"""

import pytest

from aadr.plot_pmr_maps import classify_cells


@pytest.mark.parametrize(
    "lat, lon, expected_region",
    [
        # ── Original cases ──
        (40, 80, "asia"),  # Xinjiang
        (62, 129, "asia"),  # Yakutsk (RU → asia)
        (15.21, 145.72, "oceania"),  # Saipan (US territory → oceania via admin1)
        (6.84, 158.34, "oceania"),  # Pohnpei (FM)
        (41, 29, "asia"),  # Istanbul area (TR → asia)
        (30, 31, "africa"),  # Cairo (EG)
        (61, -150, "americas"),  # Anchorage (US)
        (-33.87, 151.21, "oceania"),  # Sydney (AU)
        # ── Europe boundary / island cases ──
        (64.15, -21.95, "europe"),  # Reykjavik (IS)
        (35.90, 14.51, "europe"),  # Malta (MT)
        (62.01, 6.77, "europe"),  # Faroe Islands (FO)
        (78.22, 15.63, "europe"),  # Svalbard (SJ)
        (36.14, -5.35, "europe"),  # Gibraltar (GI)
        (37.97, 23.73, "europe"),  # Athens (GR) – near TR border
        (42.44, 19.26, "europe"),  # Montenegro (ME)
        # ── Africa edge cases ──
        (-21.11, 55.53, "africa"),  # Réunion (RE)
        (-12.78, 45.23, "africa"),  # Mayotte (YT)
        (-15.94, -5.72, "africa"),  # Saint Helena (SH)
        (36.75, 3.04, "africa"),  # Algiers (DZ) – North Africa near Europe
        (33.89, -6.92, "africa"),  # Rabat (MA) – near Strait of Gibraltar
        (-25.75, 28.19, "africa"),  # Pretoria (ZA)
        # ── Asia: Near East / Caucasus borders ──
        (33.89, 35.50, "asia"),  # Beirut (LB)
        (31.95, 35.93, "asia"),  # Amman (JO)
        (41.69, 44.80, "asia"),  # Tbilisi (GE) – Caucasus
        (40.18, 44.51, "asia"),  # Yerevan (AM)
        (40.41, 49.87, "asia"),  # Baku (AZ)
        (35.17, 33.36, "asia"),  # Nicosia (CY)
        # ── Asia: extremes ──
        (35.68, 139.69, "asia"),  # Tokyo (JP)
        (1.35, 103.82, "asia"),  # Singapore (SG)
        (47.91, 106.91, "asia"),  # Ulaanbaatar (MN)
        (37.57, 126.98, "asia"),  # Seoul (KR)
        (13.76, 100.50, "asia"),  # Bangkok (TH)
        # ── Americas: far-flung territories ──
        (18.47, -66.10, "americas"),  # San Juan (PR)
        (64.17, -51.74, "americas"),  # Nuuk, Greenland (GL)
        (18.04, -63.05, "americas"),  # Saint Martin (MF)
        (46.77, -56.18, "americas"),  # Saint Pierre and Miquelon (PM)
        (-51.80, -59.00, "americas"),  # Falkland Islands (FK)
        (4.71, -74.07, "americas"),  # Bogotá (CO)
        (-34.60, -58.38, "americas"),  # Buenos Aires (AR)
        # ── Oceania: Pacific islands ──
        (-6.31, 155.96, "oceania"),  # Bougainville, Papua New Guinea (PG)
        (-17.78, 177.96, "oceania"),  # Suva, Fiji (FJ)
        (-13.83, -171.76, "oceania"),  # Apia, Samoa (WS)
        (-22.28, 166.46, "oceania"),  # Nouméa, New Caledonia (NC)
        (13.44, 144.79, "oceania"),  # Guam (GU) – US territory in oceania
        (7.50, 134.62, "oceania"),  # Palau (PW)
        # ── US Pacific: lon/lat fallback path (lon > 144, lat < 25) ──
        (13.38, 144.95, "oceania"),  # Guam alternate coords
        (14.17, 145.25, "oceania"),  # Saipan alternate coords
        # ── Miscellaneous remote territories ──
        (-12.17, 96.83, "oceania"),  # Cocos Islands (CC → oceania)
        (-10.49, 105.63, "oceania"),  # Christmas Island (CX → oceania)
        (-53.08, 73.50, "africa"),  # French Southern Territories (TF → africa)
    ],
)
def test_classify_cells(lat: float, lon: float, expected_region: str) -> None:
    [region] = classify_cells([(lat, lon)])
    assert region == expected_region, f"({lat}, {lon}): expected {expected_region}, got {region}"


def test_classify_cells_batch_order() -> None:
    """Verify batch results match the input order."""
    cells = [(30, 31), (61, -150), (-33.87, 151.21)]
    regions = classify_cells(cells)
    assert regions == ["africa", "americas", "oceania"]
