"""Canonical public source identities and parsing for RSMViewer."""

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class SourceInfo:
    """Public source metadata with an internal provider key."""

    name: str
    provider_id: str
    source_type: str
    description: str


class SourceRegistry:
    """Registry for the four supported, publicly visible annotation sources."""

    SOURCE_MAP: Dict[str, SourceInfo] = {
        "rna3dmotifatlas": SourceInfo(
            name="RNA3DMotifAtlas",
            provider_id="bgsu_api",
            source_type="online",
            description="RNA 3D Motif Atlas annotations",
        ),
        "rfam": SourceInfo(
            name="Rfam",
            provider_id="rfam_api",
            source_type="online",
            description="Rfam annotations",
        ),
        "fr3d": SourceInfo(
            name="FR3D",
            provider_id="fr3d",
            source_type="external",
            description="FR3D annotations",
        ),
        "rnamotifscanx": SourceInfo(
            name="RNAMotifScanX",
            provider_id="rmsx",
            source_type="external",
            description="RNAMotifScanX annotations",
        ),
    }

    def __init__(self) -> None:
        self.sources = dict(self.SOURCE_MAP)

    @staticmethod
    def _key(name: str) -> str:
        return name.strip().casefold()

    def get_source(self, name: str) -> Optional[SourceInfo]:
        """Return source metadata for a case-insensitive canonical name."""
        return self.sources.get(self._key(name))

    def get_all_sources(self) -> Dict[str, SourceInfo]:
        """Return all sources keyed by normalized canonical name."""
        return dict(self.sources)

    def get_source_by_name(self, name: str) -> Optional[SourceInfo]:
        """Compatibility alias for name-based lookup."""
        return self.get_source(name)

    def get_public_name_for_provider(self, provider_id: str) -> Optional[str]:
        """Return the canonical public name for an internal provider ID."""
        provider = str(provider_id).strip()
        for source in self.sources.values():
            if source.provider_id == provider:
                return source.name
        return None

    def parse_names(self, value: str) -> List[str]:
        """Parse a comma-separated source list into canonical display names."""
        raw_names = value.split(",")
        if not value.strip() or any(not item.strip() for item in raw_names):
            raise ValueError(self._usage_error("Source names cannot be empty."))

        canonical_names: List[str] = []
        seen = set()
        for raw_name in raw_names:
            source = self.get_source(raw_name)
            if source is None:
                raise ValueError(
                    self._usage_error(f"Unknown source '{raw_name.strip()}'.")
                )
            key = self._key(source.name)
            if key in seen:
                raise ValueError(
                    self._usage_error(f"Duplicate source '{source.name}'.")
                )
            seen.add(key)
            canonical_names.append(source.name)
        return canonical_names

    def validate_source_names(self, names: Sequence[str]) -> Tuple[bool, str]:
        """Validate already separated names without changing their order."""
        try:
            self.parse_names(",".join(names))
        except ValueError as exc:
            return False, str(exc)
        return True, "OK"

    def get_provider_ids(self, names: Sequence[str]) -> List[str]:
        """Map canonical public names to internal provider IDs."""
        sources = [self.get_source(name) for name in self.get_source_names(names)]
        return [source.provider_id for source in sources if source is not None]

    def get_source_names(self, names: Sequence[str]) -> List[str]:
        """Normalize source names to canonical display form."""
        return self.parse_names(",".join(names))

    def get_source_descriptions(self, names: Sequence[str]) -> List[str]:
        """Return descriptions in the supplied source order."""
        sources = [self.get_source(name) for name in self.get_source_names(names)]
        return [source.description for source in sources if source is not None]

    def _usage_error(self, message: str) -> str:
        valid = ", ".join(source.name for source in self.sources.values())
        return f"{message} Valid sources: {valid}. Example: rmv_db RNA3DMotifAtlas,Rfam"


_source_registry = SourceRegistry()


def get_source_registry() -> SourceRegistry:
    """Return the process-wide source registry."""
    return _source_registry
