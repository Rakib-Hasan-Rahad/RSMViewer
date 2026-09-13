"""Centralized, intelligent RNA motif name normalization.

Different annotation sources spell the same structural-motif family in many
ways: abbreviations (``SR``, ``KT``, ``CL``, ``EL``), separators (``k-turn``,
``k_turn``, ``kink turn``), and family-index suffixes (``sarcin-ricin-1``,
``k-turn-2``). This module resolves every such variant to a single canonical
family name so a query for ``SR`` matches Atlas ``Sarcin-Ricin``, Rfam
``sarcin-ricin-1``/``sarcin-ricin-2``, and RNAMotifScanX ``sarcin-ricin`` rows
alike, while keeping genuinely distinct families (e.g. ``K-TURN`` vs
``REVERSE-K-TURN``) apart.

The definitions are enriched from the four sources RSMViewer integrates:
RNA 3D Motif Atlas (BGSU), Rfam, RNAMotifScanX, and FR3D.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional


def _alnum(text: object) -> str:
    """Uppercase, alphanumeric-only normalization used for alias lookup."""
    return re.sub(r"[^A-Z0-9]", "", str(text).upper())


# Canonical display name -> accepted free-form spellings/abbreviations.
# Values are matched after alphanumeric normalization, so separators and case
# do not matter here. Add new families or aliases in one place only.
_FAMILY_DEFINITIONS = {
    # ── Core tertiary motifs (cross-source) ──────────────────────────────
    "SARCIN-RICIN": [
        "SR", "SARCIN", "SARCINRICIN", "SARCINRICINLOOP", "SRL",
        "SARCINRICIN1", "SARCINRICIN2", "SARCINRICINMOTIF",
    ],
    "K-TURN": [
        "KT", "KTURN", "KINKTURN", "KINK",
        "KTURN1", "KTURN2", "KINKTURN1", "KINKTURN2",
    ],
    "REVERSE-K-TURN": [
        "REVERSEKTURN", "REVERSEKINKTURN", "REVKTURN", "REVERSEKT", "RKT",
        "REVERSEKINK",
    ],
    "PSEUDO-KINK-TURN": ["PK-TURN", "PKTURN", "PKT", "PSEUDOKINKTURN"],
    "C-LOOP": ["CL", "CLOOP"],
    "E-LOOP": ["EL", "ELOOP"],
    "T-LOOP": ["TL", "TLOOP"],
    "U-TURN": ["UT", "UTURN"],
    # ── Tetraloops ───────────────────────────────────────────────────────
    "GNRA": ["GNRA", "GNRATETRALOOP"],
    "UNCG": ["UNCG", "UNCGTETRALOOP"],
    "CUYG": ["CUYG"],
    # ── Secondary-structure context (Atlas generic keys) ─────────────────
    "HL": ["HL", "HAIRPIN", "HAIRPINLOOP", "HAIRPINLOOPHL", "OTHERHL"],
    "IL": ["IL", "INTERNAL", "INTERNALLOOP", "INTERNALLOOPIL", "OTHERIL"],
    "J3": ["J3", "3WJ", "THREEWAYJUNCTION", "3WAYJUNCTIONJ3", "3WAYJUNCTION"],
    "J4": ["J4", "4WJ", "FOURWAYJUNCTION", "4WAYJUNCTIONJ4", "4WAYJUNCTION"],
    "J5": ["J5", "5WJ", "FIVEWAYJUNCTION", "5WAYJUNCTIONJ5", "5WAYJUNCTION"],
    "J6": ["J6", "6WJ", "SIXWAYJUNCTION", "6WAYJUNCTIONJ6", "6WAYJUNCTION"],
    "J7": ["J7", "7WJ", "SEVENWAYJUNCTION", "7WAYJUNCTIONJ7", "7WAYJUNCTION"],
    "PSEUDOKNOT": ["PSEUDOKNOT", "PK"],
    "BULGED": ["BULGED", "BULGE", "BULGEDLOOP"],
    # ── Rfam / other named families ──────────────────────────────────────
    "TANDEM-GA": ["TANDEMGA", "TANDEMGASHEARED"],
    "RIGHT-ANGLE": ["RIGHTANGLE", "RIGHTANGLE2", "RIGHTANGLE3"],
    "DOCKING-ELBOW": ["DOCKINGELBOW"],
    "TWIST-UP": ["TWISTUP"],
    "UAA-GAN": ["UAAGAN", "UAAGANMOTIF"],
    "DOMAIN-V": ["DOMAINV"],
    "SRP-S-DOMAIN": ["SRPSDOMAIN", "SRPS"],
    "ANYA": ["ANYA"],
    "UMAC": ["UMAC"],
    "TRIT": ["TRIT"],
    # ── Atlas descriptive families (as annotated by BGSU) ────────────────
    "ANTICODON LOOP RELATED": ["ANTICODONLOOPRELATED"],
    "LSU P LOOP": ["LSUPLOOP"],
    "LSU A LOOP": ["LSUALOOP"],
    "RIBOSOMAL LSU H95": ["RIBOSOMALLSUH95", "RIBSOMALLSUH95", "LSUH95"],
    "5S RRNA JUNCTION": ["5SRRNAJUNCTION"],
    "TRIPLE SHEARED": ["TRIPLESHEARED"],
    "DOUBLE SHEARED": ["DOUBLESHEARED"],
    "DOUBLE SHEARED WITH NON-CANONICAL CWW": ["DOUBLESHEAREDWITHNONCANONICALCWW"],
    "DOUBLE SHEARED; A IN SYN": ["DOUBLESHEAREDAINSYN"],
    "SINGLE STACK BEND": ["SINGLESTACKBEND"],
    "MAJOR GROOVE PLATFORM": ["MAJORGROOVEPLATFORM"],
    "MINOR GROOVE PLATFORM": ["MINORGROOVEPLATFORM"],
    "MAJOR GROOVE INTERCALATION": ["MAJORGROOVEINTERCALATION"],
    "TANDEM NON-CANONICAL CWW PAIRS": ["TANDEMNONCANONICALCWWPAIRS"],
    "TRIPLE NON-CANONICAL CWW PAIRS": ["TRIPLENONCANONICALCWWPAIRS"],
    "ISOLATED NON-CANONICAL CWW PAIR": ["ISOLATEDNONCANONICALCWWPAIR"],
    "ISOLATED NON-CANONICAL CWW WITH BULGES": ["ISOLATEDNONCANONICALCWWWITHBULGES"],
    "ISOLATED NON-CANONICAL CWW CONTACT": ["ISOLATEDNONCANONICALCWWCONTACT"],
    "NON CANONICAL CWW AND NEAR PAIR": ["NONCANONICALCWWANDNEARPAIR"],
    "ISOLATED CWS BASEPAIR": ["ISOLATEDCWSBASEPAIR"],
    "ISOLATED CWH BASEPAIR": ["ISOLATEDCWHBASEPAIR"],
    "ISOLATED THS BASEPAIR WITH BULGES": ["ISOLATEDTHSBASEPAIRWITHBULGES"],
    "ISOLATED TWW TURN": ["ISOLATEDTWWTURN"],
    "INTERCALATED TWH": ["INTERCALATEDTWH"],
    "AAA CROSS-STRAND STACK": ["AAACROSSSTRANDSTACK"],
    "THW-THW CROSS-STRAND STACK": ["THWTHWCROSSSTRANDSTACK"],
    "AG THS OUTSIDE LOOP": ["AGTHSOUTSIDELOOP"],
    "THS DOUBLE PLATFORM": ["THSDOUBLEPLATFORM"],
    "EXTERNALLY STRUCTURED": ["EXTERNALLYSTRUCTURED"],
    "180 DEGREE TURN": ["180DEGREETURN"],
    "TSH-TWH-THS": ["TSHTWHTHS"],
    "TSH-THH-THS": ["TSHTHHTHS"],
    "TSH-THS-THW": ["TSHTHSTHW"],
    "TSH-TSH-THH-THS": ["TSHTSHTHHTHS"],
    "TSH-THW-THH-THS": ["TSHTHWTHHTHS"],
    "TWH-TWH-THW-THW": ["TWHTWHTHWTHW"],
}

# Reverse index: normalized alias -> canonical display name.
_ALIAS_TO_CANONICAL = {}


def _register_alias(alias: str, canonical: str) -> None:
    """Register an alias, raising if it already maps to a different family.

    Silent overwrites would let a future edit point one spelling at two
    families; failing loudly at import time keeps the table unambiguous.
    """
    key = _alnum(alias)
    if not key:
        return
    existing = _ALIAS_TO_CANONICAL.get(key)
    if existing is not None and existing != canonical:
        raise ValueError(
            f"Alias collision: {alias!r} maps to both {existing!r} and {canonical!r}"
        )
    _ALIAS_TO_CANONICAL[key] = canonical


for _canonical, _aliases in _FAMILY_DEFINITIONS.items():
    _register_alias(_canonical, _canonical)
    for _alias in _aliases:
        _register_alias(_alias, _canonical)


# ── Free-text keyword scan (used for FR3D query names) ───────────────────────
# FR3D query names embed the family loosely, e.g. "geometric_5_sarcin_ricin" or
# "geometric_3_sarcin3geometric". These distinctive keywords are scanned as
# substrings of the alnum-normalized text. Ordered most-specific first so that
# "reverse kink" resolves to REVERSE-K-TURN before "kink" can match K-TURN.
# Short/ambiguous abbreviations (SR, KT, HL, IL, CL, EL) are deliberately
# excluded here because they would match inside unrelated words.
_TEXT_KEYWORDS = [
    ("REVERSEKINKTURN", "REVERSE-K-TURN"),
    ("REVERSEKTURN", "REVERSE-K-TURN"),
    ("REVERSEKINK", "REVERSE-K-TURN"),
    ("SARCINRICIN", "SARCIN-RICIN"),
    ("SARCIN", "SARCIN-RICIN"),
    ("KINKTURN", "K-TURN"),
    ("KTURN", "K-TURN"),
    ("CLOOP", "C-LOOP"),
    ("ELOOP", "E-LOOP"),
    ("TLOOP", "T-LOOP"),
    ("UTURN", "U-TURN"),
    ("GNRA", "GNRA"),
    ("UNCG", "UNCG"),
    ("CUYG", "CUYG"),
    ("TANDEMGA", "TANDEM-GA"),
    ("RIGHTANGLE", "RIGHT-ANGLE"),
    ("DOCKINGELBOW", "DOCKING-ELBOW"),
    ("TWISTUP", "TWIST-UP"),
    ("UAAGAN", "UAA-GAN"),
    ("PSEUDOKNOT", "PSEUDOKNOT"),
]


def canonical_motif(value: str) -> str:
    """Return the canonical family display name for any spelling/abbreviation.

    Unknown families fall back to their upper-cased form with a trailing
    ``-<n>`` / ``_<n>`` family index removed, so novel motifs still normalize
    consistently and match themselves.
    """
    raw = str(value).strip()
    if not raw:
        return raw
    key = _alnum(raw)
    if key in _ALIAS_TO_CANONICAL:
        return _ALIAS_TO_CANONICAL[key]
    stripped = re.sub(r"\d+$", "", key)
    if stripped and stripped in _ALIAS_TO_CANONICAL:
        return _ALIAS_TO_CANONICAL[stripped]
    return re.sub(r"[-_]\d+$", "", raw.upper())


def is_known_family(value: str) -> bool:
    """True when *value* resolves to a family in the alias table."""
    key = _alnum(value)
    if key in _ALIAS_TO_CANONICAL:
        return True
    stripped = re.sub(r"\d+$", "", key)
    return bool(stripped) and stripped in _ALIAS_TO_CANONICAL


def family_short_code(value: str) -> str:
    """Return a short, group-name-friendly code for a family (e.g. ``SR``).

    Prefers the shortest defined abbreviation (``SR`` for Sarcin-Ricin, ``KT``
    for K-Turn, ``CL``/``EL`` for the loops). Families without a short alias
    fall back to their word initials, then to a sanitized short form, so the
    suggested group names stay compact (``group_SR`` rather than
    ``group_Sarcin_Ricin``).
    """
    canonical = canonical_motif(value)
    aliases = _FAMILY_DEFINITIONS.get(canonical, [])
    if aliases:
        short = min(aliases, key=lambda alias: (len(alias), alias))
        if len(short) <= 4:
            return short.upper()
    words = [w for w in re.split(r"[^A-Za-z0-9]+", canonical) if w]
    initials = "".join(w[0] for w in words)
    if 2 <= len(initials) <= 5:
        return initials.upper()
    return (re.sub(r"[^A-Za-z0-9]+", "_", canonical).strip("_") or "motif").upper()


def family_from_text(text: str) -> Optional[str]:
    """Infer a canonical family from loose free text such as an FR3D query name.

    Scans distinctive family keywords as substrings, most-specific first, so
    ``geometric_5_sarcin_ricin`` and ``geometric_3_sarcin3geometric`` both map
    to ``SARCIN-RICIN`` while ``reverse_kturn_*`` maps to ``REVERSE-K-TURN`` and
    never to ``K-TURN``. Returns None when no keyword is found.
    """
    if not text:
        return None
    norm = _alnum(text)
    if not norm:
        return None
    for keyword, canonical in _TEXT_KEYWORDS:
        if keyword in norm:
            return canonical
    return None


def labels_match_motif(
    query_motif: str,
    labels: Iterable[str],
    free_text_labels: Iterable[str] = (),
) -> bool:
    """True when any source label denotes the same family as *query_motif*.

    Matching is exact on canonical family, which keeps distinct families such
    as ``K-TURN`` and ``REVERSE-K-TURN`` separate. A lenient substring fallback
    is applied only to labels that are not themselves a *different* known
    family, so compound or novel labels (e.g. ``sarcin-ricin core region``)
    still match without letting abbreviations bleed across families.

    ``free_text_labels`` are matched with :func:`family_from_text` instead of
    exact canonicalization. Pass FR3D query names here so a query for ``SR``
    matches an FR3D motif named ``geometric_5_sarcin_ricin`` by keyword, while
    other sources keep strict family matching.
    """
    target = canonical_motif(query_motif)
    target_norm = _alnum(target)
    for label in labels:
        if not label:
            continue
        label_canonical = canonical_motif(label)
        if label_canonical == target:
            return True
        # Substring fallback for compound/novel labels only (e.g.
        # "sarcin-ricin core region"). Restricted to canonical names of at
        # least 5 characters so short codes (HL, IL, J3, SR, KT) cannot bleed
        # into unrelated words such as "THLXYZ" or "SILENCER".
        if (
            len(target_norm) >= 5
            and target_norm in _alnum(label)
            and not is_known_family(label_canonical)
        ):
            return True
    for text in free_text_labels:
        if family_from_text(text) == target:
            return True
    return False
