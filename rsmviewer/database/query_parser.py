"""Parser for the annotation-selection command grammar."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import FrozenSet, List, Sequence, Tuple, Union

from .source_registry import SourceRegistry, get_source_registry


class QuerySyntaxError(ValueError):
    """A query error with a source-expression location."""

    def __init__(self, message: str, expression: str, position: int = 0) -> None:
        self.position = max(0, min(position, len(expression)))
        pointer = " " * self.position + "^"
        super().__init__(f"{message}\n{expression}\n{pointer}")


@dataclass(frozen=True)
class SourcePredicate:
    operator: str
    value: str = ""
    left: "SourcePredicate | None" = None
    right: "SourcePredicate | None" = None

    def matches(self, available_sources: Sequence[str]) -> bool:
        available = {source.casefold() for source in available_sources}
        if self.operator == "name":
            return self.value.casefold() in available
        if self.operator == "not":
            return not self.left.matches(available_sources)
        if self.operator == "and":
            return self.left.matches(available_sources) and self.right.matches(available_sources)
        if self.operator == "or":
            return self.left.matches(available_sources) or self.right.matches(available_sources)
        raise RuntimeError(f"Unknown source predicate operator: {self.operator}")


@dataclass(frozen=True)
class QueryExpression:
    motif: str
    structures: Union[str, Tuple[str, ...]]
    sources: SourcePredicate
    group: str
    text: str

    def matches_structure(self, structure_id: str) -> bool:
        if self.structures == "all":
            return True
        return structure_id.upper() in self.structures


_GROUP_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def canonical_motif_name(value: str) -> str:
    """Normalize documented short motif aliases without imposing biology."""
    name = value.strip()
    if not name:
        raise ValueError("Motif expression cannot be empty")
    from .motif_aliases import canonical_motif
    return canonical_motif(name)


def parse_selection_query(
    text: str, registry: SourceRegistry | None = None
) -> QueryExpression:
    """Parse ``motif, structures, sources, as group``."""
    query = text.strip()
    clauses = [clause.strip() for clause in query.split(",")]
    if len(clauses) != 4:
        raise QuerySyntaxError(
            "Expected four comma-separated clauses: motif, structures, sources, as group.",
            query,
            query.find(",") if "," in query else len(query),
        )

    motif = canonical_motif_name(clauses[0])
    structures = _parse_structure_expression(clauses[1], query)
    source_text, group = _split_source_and_group(clauses[2], clauses[3], query)
    source_parser = _SourceExpressionParser(source_text, registry or get_source_registry())
    sources = source_parser.parse()
    return QueryExpression(motif, structures, sources, group, query)


def _split_source_and_group(
    source_clause: str, group_clause: str, full_query: str
) -> Tuple[str, str]:
    match = re.match(r"^(.+?)\s+as\s+([A-Za-z_][A-Za-z0-9_]*)$", group_clause, re.IGNORECASE)
    if match:
        source_clause = f"{source_clause} {match.group(1)}"
        group = match.group(2)
    else:
        match = re.match(r"^as\s+([A-Za-z_][A-Za-z0-9_]*)$", group_clause, re.IGNORECASE)
        if not match:
            raise QuerySyntaxError(
                "Expected a group name using 'as <group>'.",
                full_query,
                full_query.lower().rfind("as") if "as" in full_query.lower() else len(full_query),
            )
        group = match.group(1)

    if not _GROUP_RE.fullmatch(group):
        raise QuerySyntaxError("Invalid group name.", full_query, full_query.rfind(group))
    return source_clause.strip(), group


def _parse_structure_expression(value: str, full_query: str) -> Union[str, Tuple[str, ...]]:
    expression = value.strip()
    if expression.casefold() == "all":
        return "all"
    tokens = re.split(r"\s+(and|or)\s+", expression, flags=re.IGNORECASE)
    identifiers = [tokens[index].strip().upper() for index in range(0, len(tokens), 2)]
    operators = [tokens[index].casefold() for index in range(1, len(tokens), 2)]
    if not identifiers or any(not re.fullmatch(r"[A-Z0-9][A-Z0-9_-]*", identifier) for identifier in identifiers):
        raise QuerySyntaxError("Invalid structure expression.", full_query, full_query.find(value))
    if len(identifiers) > 1 and len(operators) != len(identifiers) - 1:
        raise QuerySyntaxError("Invalid structure expression operators.", full_query, full_query.find(value))
    return tuple(dict.fromkeys(identifiers))


class _SourceExpressionParser:
    def __init__(self, expression: str, registry: SourceRegistry) -> None:
        self.expression = expression.strip()
        self.registry = registry
        self.tokens = re.findall(r"[A-Za-z0-9]+|[^\sA-Za-z0-9]", self.expression)
        self.index = 0

    def parse(self) -> SourcePredicate:
        if not self.tokens:
            raise QuerySyntaxError("Source expression cannot be empty.", self.expression)
        result = self._parse_or()
        if self.index != len(self.tokens):
            raise QuerySyntaxError("Unexpected source-expression token.", self.expression, self._position())
        return result

    def _parse_or(self) -> SourcePredicate:
        result = self._parse_and()
        while self._accept("or"):
            result = SourcePredicate("or", left=result, right=self._parse_and())
        return result

    def _parse_and(self) -> SourcePredicate:
        result = self._parse_not()
        while self._accept("and"):
            result = SourcePredicate("and", left=result, right=self._parse_not())
        return result

    def _parse_not(self) -> SourcePredicate:
        if self._accept("not"):
            return SourcePredicate("not", left=self._parse_not())
        return self._parse_primary()

    def _parse_primary(self) -> SourcePredicate:
        if self._accept("("):
            result = self._parse_or()
            if not self._accept(")"):
                raise QuerySyntaxError("Missing ')' in source expression.", self.expression, self._position())
            return result
        if self.index >= len(self.tokens):
            raise QuerySyntaxError("Expected a source name.", self.expression, self._position())
        token = self.tokens[self.index]
        self.index += 1
        source = self.registry.get_source(token)
        if source is None:
            raise QuerySyntaxError(f"Unknown source '{token}'.", self.expression, self._position() - len(token))
        return SourcePredicate("name", value=source.name)

    def _accept(self, expected: str) -> bool:
        if self.index < len(self.tokens) and self.tokens[self.index].casefold() == expected.casefold():
            self.index += 1
            return True
        return False

    def _position(self) -> int:
        return len(" ".join(self.tokens[:self.index]))